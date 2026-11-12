###
#
###

import pandas as pd
import csv
import itertools

configfile: "./config.yaml"

NCBI_KEY = os.environ.get("NCBI_API_KEY")

def sniff_delimiter(file_path, sample_size=2048):
    with open(file_path, 'r') as fp:
        sample = fp.read(sample_size)
        return csv.Sniffer().sniff(sample).delimiter

SAMPLES = []
SYMBOLS = []
max_rows = 5 if config.get("test") else None
sample_delim = sniff_delimiter(config['sample_metadata'])

with open(config['sample_metadata'], mode='r', newline='') as fp:
    reader = csv.DictReader(fp, delimiter=sample_delim)
    for row in itertools.islice(reader, max_rows):
        SAMPLES.append(row['run'])

OUTDIR = [ config.get('output_directory') if config.get('output_directory') is not None else './workflow-results' ]

dbdir = config.get('database_directory')

# Dictionary for correct slurm batch allocations
HICPUS_JOBS = {1: ['low', 1], 2: ['low', 1], 3: ['low', 25], 4: ['high', 25], 5: ['high', 50]}
BIGMEM_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bml', 25], 4: ['bmh', 25], 5: ['bmh', 50]}

wildcard_constraints:
    sample = "(SRR|ERR|DRR)\d+",

metadata_df = pd.read_csv(config['protein_metadata'], header=None, names=['protein_id'])
row_count = len(metadata_df)

# Conditionally define the output filename based on row count
if row_count == 1:
    # Extract the value from the single row using the assigned column name
    protein_name = metadata_df.iloc[0]['protein_id']
    protein_filename = protein_name.replace(" ", "_")
    target_output = f"{{o}}/final_results_{protein_filename}.csv"
    target_protein = f"{{o}}/protein_seqs_{protein_filename}.faa.gz"
else:
    # Use a default filename if there is more than one row
    target_output = "{o}/final_results_multi.csv"
    target_protein = "{o}/protein_seqs_multi.faa.gz"

rule all:
    input:
        #expand('{o}/cds_info/cds_idents.txt', o=OUTDIR),
        #expand('{o}/cds_info/cds_idents_data.zip', o=OUTDIR),
        #expand("{o}/cds_output/combined_cds.fna", o=OUTDIR),
        #expand("{o}/cds_output/combined_prot.faa", o=OUTDIR),
        #expand("{o}/blast_dbs/cds_db.pdb", o=OUTDIR),
        #expand("{o}/blast_results/blastp_output.csv", o=OUTDIR),
        #expand("{o}/cds_output/filtered_cds.fna", o=OUTDIR),
        #expand("{o}/mapped_reads/stats/{sample}_mapping_stats.tsv", o=OUTDIR, sample=SAMPLES),
        #expand("{o}/mapped_reads/stats/{sample}_coverage_stats.tsv", o=OUTDIR, sample=SAMPLES),
        #expand("{o}/mapped_reads/stats/{sample}_average_coverage_stats.tsv", o=OUTDIR, sample=SAMPLES),
        expand(target_protein, o=OUTDIR),
        #expand("{o}/results/filtered_{sample}_merged_stats.csv", o=OUTDIR, sample=SAMPLES),
        expand(target_output, o=OUTDIR),


rule download_sra:
    output:
        sra_file = temp("{o}/sra/{sample}.sra"),
        checksum = "{o}/validate/{sample}-checksum.log",
        stats = '{o}/stats/{sample}-counts.xml',
    threads:
        10
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition = lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    conda:
        "envs/sra.yaml"
    shell:
       """
           (
               aws s3 cp --no-sign-request \
                   s3://sra-pub-run-odp/sra/{wildcards.sample}/{wildcards.sample} \
                   {output.sra_file} \
               || prefetch {wildcards.sample} -o {output.sra_file}
           )
           # Validate the sra file against the checksums provided
           vdb-validate {output.sra_file} -I no &> {output.checksum}

           if grep -q 'err' {output.checksum};
             then
               echo 'Validation of {wildcards.sample} failed. Removing {output.sra_file}...'
               rm {output.sra_file}
               echo 'Validation failed, exiting.'
               exit 1
           else
               echo 'Sample {wildcards.sample} validated without error'
           fi

           # Stats for the sequence contained in the sra file
           sra-stat --alignment off --quick -x {output.sra_file} &> {output.stats}
       """

rule dump_fastq:
    input:
        xml_script = 'scripts/process_xml_file.py',
        seq_script = 'scripts/process_fastq_files.sh',
        sra_file = "{o}/sra/{sample}.sra",
        stats = '{o}/stats/{sample}-counts.xml',
    output:
        r1 = temp("{o}/fastq/{sample}_1.fastq.gz"),
        r2 = temp("{o}/fastq/{sample}_2.fastq.gz"),
    conda:
        "envs/sra.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition = lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    threads:
        16
    params:
        sra_counts = 'all-sra-counts.csv',
        fastq_counts = 'all-fastq-counts.csv',
    shell:
        """
        # job ID returning 0 in log files when running on cluster (use SLURM_JOB_ID if available, otherwise use Snakemake jobid)
        job=${{SLURM_JOB_ID:-{jobid}}}
        # make a directory specific to user and job
        export MYTMP="/scratch/${{USER}}/slurm_{rule}.${{job}}"
        mkdir -v -p "$MYTMP"/scripts

        # force clean it up after job script ends
        function cleanup() {{
            echo "Cleaning up $MYTMP"
            cd /tmp || cd /
            rm -rf "$MYTMP" && echo "$MYTMP removed"
        }}
        trap cleanup EXIT

        # change to it!
        cp {input.xml_script} "$MYTMP"/scripts
        cp {input.seq_script} "$MYTMP"/scripts
        cd "$MYTMP"

        # run your code, telling it to use that dir:
        echo "Changing to $(pwd)"

        echo "Running fasterq-dump for {wildcards.sample}"
        echo ''

        # Get fastqc file
        # To verify any errors due to threads -> https://github.com/sourmash-bio/sourmash/issues/2941
        # Quality check, if failed check use Worts conservative approach
        # https://github.com/sourmash-bio/wort/blob/09d6ee552c85360fdc60e9d59e2fc757c2c6dd99/wort/blueprints/compute/tasks.py#L38-L61
        fasterq-dump {input.sra_file} --skip-technical --split-files --progress \
                     --threads {threads} --bufsize 1000MB --curcache 10000MB --mem {resources.mem_mb}

        echo ''
        echo "Fastq files generated by fasterq-dump:"
        ls -h "$MYTMP"
        echo ''

        # For checking the status of fastq files because the md5checksums are difficult to use with this file type
        {input.xml_script} {input.stats} -o {params.sra_counts}
        {input.xml_script} {input.stats} -o {wildcards.o}/stats/{params.sra_counts}

        # capture the awk variables in the 3rd and 4th column of the counts file for each sample
        read sra_read_count sra_base_count < <(awk -F',' -v sample="{wildcards.sample}" '
            $0 ~ sample {{
                sra_read_count = $3
                sra_base_count = $4
                print sra_read_count, sra_base_count
                exit
            }}
        ' {params.sra_counts})

        sra_base_count=$(echo "$sra_base_count" | tr -d '\r')

        R1=$MYTMP/{wildcards.sample}_1.fastq
        R2=$MYTMP/{wildcards.sample}_2.fastq
        R1_T=$MYTMP/{wildcards.sample}_1_trim.fastq
        R2_T=$MYTMP/{wildcards.sample}_2_trim.fastq

        {input.seq_script} {threads} ./{params.fastq_counts} "$R1" "$R2"
        # Only add the output if it is not in the output...
        if ! grep -q "{wildcards.sample}_1" ./{params.fastq_counts} && ! grep -q "{wildcards.sample}_2" ./{params.fastq_counts}; then
            cat ./{params.fastq_counts} >> {wildcards.o}/stats/{params.fastq_counts}
        fi

        # Create the count variables from the file
        read fq1_read_count fq1_base_count < <(awk -F',' -v sample="{wildcards.sample}_1" '
            $0 ~ sample {{
                fq1_read_count = $2
                fq1_base_count = $3
                print fq1_read_count, fq1_base_count
                exit
            }}
        ' {params.fastq_counts})
        read fq2_read_count fq2_base_count < <(awk -F',' -v sample="{wildcards.sample}_2" '
            $0 ~ sample {{
                fq2_read_count = $2
                fq2_base_count = $3
                print fq2_read_count, fq2_base_count
                exit
            }}
        ' {params.fastq_counts})

        # if there is a comma, remove the comma
        fq1_read_count=$(echo "$fq1_read_count" | tr -d ',')
        fq2_read_count=$(echo "$fq2_read_count" | tr -d ',')
        fq1_base_count=$(echo "$fq1_base_count" | tr -d ',')
        fq2_base_count=$(echo "$fq2_base_count" | tr -d ',')

        echo "$fq1_read_count reads in {wildcards.sample}_1"
        echo "$fq2_read_count reads in {wildcards.sample}_2"
        echo "$sra_read_count reads in sra file"
        echo "$(( fq1_base_count + fq2_base_count )) bases in {wildcards.sample}_1 and {wildcards.sample}_2"
        echo "$sra_base_count in sra file"

        # Compare the SRA file counts and the FASTQ file counts
        if [[ $fq1_read_count == $sra_read_count && $fq2_read_count == $sra_read_count && $(( fq1_base_count + fq2_base_count )) == $sra_base_count ]]; then
            echo "Read and base counts match!!!"

       # Remove and use the unthreaded, tried-and-true approach
        else
            echo "Read or base counts do not match!!!"
            rm "$R1" "$R2"
            fastq-dump --disable-multithreading --skip-technical --split-files --clip {input.sra_file} -O "$MYTMP"
        fi

        echo ''
        echo "Processing R1..."
        seqtk seq -C "$R1" | \
            perl -ne 's/\.([12])$/\/$1/; print $_' | \
            gzip -9 -c > {output.r1} &

        echo "Processing R2..."
        seqtk seq -C "$R2" | \
            perl -ne 's/\.([12])$/\/$1/; print $_' | \
            gzip -9 -c > {output.r2} &

        wait
        echo "Finished downloading raw reads"
        """


rule fastp:
  input:
    fq1 = "{o}/fastq/{sample}_1.fastq.gz",
    fq2 = "{o}/fastq/{sample}_2.fastq.gz"
  output: 
    trim1 = temporary("{o}/fastp_output/trim/{sample}_1_trim.fastq.gz"),
    trim2 = temporary("{o}/fastp_output/trim/{sample}_2_trim.fastq.gz"),
    json = "{o}/fastp_output/reports/{sample}_trim.json",  
    html = "{o}/fastp_output/reports/{sample}_trim.html",
    log = temp("{o}/fastp_output/log/{sample}_stderr.log"),
  resources:
    mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    time = lambda wildcards, attempt: 12 * 60 * attempt,
    runtime = lambda wildcards, attempt: 12 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
        echo "Running fastp normally..."
        fastp \
          --in1 {input.fq1} \
          -I {input.fq2} \
          --out1 {output.trim1} \
          -O {output.trim2} \
          --json {output.json} \
          --html {output.html} 2> {output.log} || fallback=1

        if grep -q "igzip: unexpected eof" {output.log} ; then
            echo "igzip error detected — retrying with manual decompression..."
            fastp \
              --in1 <(gunzip -c {input.fq1}) \ #decompress and read in
              -I <(gunzip -c {input.fq2}) \ #decompress and read in
              --out1 {output.trim1} \
              -O {output.trim2} \
              --json {output.json} \
              --html {output.html}
        fi
    """

rule fastqc:
  input:
    trim1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq.gz",
    trim2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq.gz",
  output:
    html1 = "{o}/fastqc_reports/{sample}_1_fastqc.html",
    html2 = "{o}/fastqc_reports/{sample}_2_fastqc.html",
    zip1 = "{o}/fastqc_reports/{sample}_1_fastqc.zip",
    zip2 = "{o}/fastqc_reports/{sample}_2_fastqc.zip",
  resources:
    mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    time = lambda wildcards, attempt: 15 * 60 * attempt,
    runtime = lambda wildcards, attempt: 15 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    fastqc \
      {input.trim1} \
      {input.trim2} \
      --outdir {FASTQC_DIR}
    """

rule sketch:
    input: 
        trim1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq.gz",
        trim2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq.gz",
    output:
        trim1 = "{o}/sigs/{sample}_1_trim.zip",
        trim2 = "{o}/sigs/{sample}_2_trim.zip",
        sig = "{o}/sigs/{sample}_trim.k31.zip",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 15 * 60 * attempt,
        runtime = lambda wildcards, attempt: 15 * 60 * attempt,
        allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
        partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
    conda: "envs/branchwater.yaml"
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in config["k_sizes"]]),
        scale = config["scale"],
    shell:"""
        sourmash scripts singlesketch --input-moltype dna \
             -p scaled={params.scale},{params.k_list},abund \
             {input.trim1} -o {output.trim1}
        sourmash scripts singlesketch --input-moltype dna \
             -p scaled={params.scale},{params.k_list},abund \
             {input.trim2} -o {output.trim2}
        sourmash sig merge {output.trim1} {output.trim2} \
           -o {output.sig} -k 31 --name 'rawreads:{wildcards.sample}'
    """

rule pathlist:
    input:
        sig = expand("{o}/sigs/{sample}_trim.k31.zip", o=OUTDIR, sample=SAMPLES),
    output:
        "{o}/summary-manifest.csv"
    conda: 'envs/branchwater.yaml'
    shell:"""
        sourmash sig collect {input.sig} -o {output} -F csv --abspath
    """

rule fastmultigather:
    input:
        sig = "{o}/summary-manifest.csv",
	db = f"{dbdir}/gtdb-rs226-k31.rocksdb",
    output:
        gather = '{o}/gather/gather-k31.csv',
    conda:
        'envs/branchwater.yaml'
    threads:
        16
    resources:
        disk_mb = lambda wildcards, attempt: 96 * 1024 * attempt,
        mem_mb = lambda wildcards, attempt: 96 * 1024 * attempt,
        time = lambda wildcards, attempt: 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 24 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    params:
        scale = config["scale"],
    shell: """
        sourmash scripts fastmultigather {input.sig} {input.db} -c {threads} \
            -k 31 -s {params.scale} -t 0 -m DNA \
            -o {output.gather}
    """

rule collect_cds:
    input:
        gat = '{o}/gather/gather-k31.csv'
    output:
        cds_idents = '{o}/cds_info/cds_idents.txt',
    shell:"""
        awk -F, '{{
            gsub(/"/, "", $10)

            if (match($10, /GC[^ \t"]*/)) {{
                print substr($10, RSTART, RLENGTH)
            }}
        }}' {input} | sort -u > {output.cds_idents}
    """    

rule get_cds:
    input:
        cds_idents = '{o}/cds_info/cds_idents.txt',
    output:
        data = '{o}/cds_info/cds_idents_data.zip'
    conda: 'envs/ncbi.yaml'
    params:
        apikey=NCBI_KEY,
        meta = '{o}/cds_info/cds_idents_metadata.zip',
        data = '{o}/cds_info/cds_idents_data',
    threads:
        32
    resources:
        disk_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
        mem_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
        time = lambda wildcards, attempt: 6 * 60 * attempt,
        runtime = lambda wildcards, attempt: 6 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    shell:"""
        echo "Gathering metadata..."

        # gather the cds nucl and protein aa metadata
        datasets download genome accession \
            --inputfile "{input.cds_idents}" \
            --include cds,protein \
            --api-key "{params.apikey}" \
            --filename "{params.meta}" \
            --dehydrated

        echo "Metadata complete. Gathering data..."

        unzip -o {params.meta} -d {params.data}

        datasets rehydrate --directory {params.data} \
           --api-key "{params.apikey}" \
           --max-workers 30

        cd {wildcards.o}/cds_info && \
            zip -r -0 -q {output.data} "cds_idents_data" && \
            cd -

        echo "Combined zip created {output.data}"
        """

rule combine_cds:
    input:
        cds = '{o}/cds_info/cds_idents_data.zip',
    output:
        combined_cds = "{o}/cds_output/combined_cds.fna",
        combined_prot = "{o}/cds_output/combined_prot.faa",
    conda: "envs/sebastian_env.yml"
    resources:
        mem_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
        time = lambda wildcards, attempt: 2 * 60 * attempt,
        runtime = lambda wildcards, attempt: 2 * 60 * attempt,
        allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
        partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
    shell:"""
        unzip -Z1 {input.cds} | grep 'cds_from_genomic\.fna$' | while IFS= read -r file; do
            base=$(basename "$(dirname "$file")")  # e.g., GCF_964240095.1
            unzip -p {input.cds} "$file" | \
                awk -v prefix="$base" '/^>/ {{print ">" prefix "_" substr($0, 2)}} !/^>/' 
        done > 	{output.combined_cds}

        unzip -Z1 {input.cds} | grep 'protein\.faa$' | while IFS= read -r file; do
            base=$(basename "$(dirname "$file")")  # e.g., GCF_964240095.1
            unzip -p {input.cds} "$file" | \
                awk -v prefix="$base" '/^>/ {{print ">" prefix "_" substr($0, 2)}} !/^>/' 
        done > 	{output.combined_prot}
    """


rule blast_cds_database:
    input:
        cds = "{o}/cds_output/combined_prot.faa"
    output:
        pjs = "{o}/blast_dbs/cds_db.pjs",
        pdb = "{o}/blast_dbs/cds_db.pdb"
    params:
        out_prefix = "{o}/blast_dbs/cds_db"
    resources:
        mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
        time = lambda wildcards, attempt: 8 * 60 * attempt,
        runtime = lambda wildcards, attempt: 8 * 60 * attempt,
        allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
        partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
    conda: "envs/sebastian_env.yml"
    shell:
      """
      makeblastdb -in {input.cds} -dbtype prot -out {params.out_prefix}
      """

def build_ncbi_query(wildcards, input):
    with open(input[0], "r") as f:
        proteins = [line.strip() for line in f if line.strip()]
    
    or_terms = [f'"{p}"[Protein Name]' for p in proteins]
    combined_queries = " OR ".join(or_terms)
    
    return f"({combined_queries}) AND (Bacteria[Organism]) AND (3[Protein Count] : 1000000[Protein Count])"


rule get_genes:
    input:
        csv = config['protein_metadata'],
    output:
        pgenes = target_protein
    conda: 'envs/ncbi.yaml'
    resources:
        time = lambda wildcards, attempt: 18 * 60 * attempt,
        runtime = lambda wildcards, attempt: 18 * 60 * attempt,
    params:
        apikey=NCBI_KEY,
        query_string = build_ncbi_query
    shell:
        """
        echo "Sending batched query to NCBI..."
        
        # ipg: identical protein group or prot: protein
        # ipg are identical sequences found in the prot db
        esearch -db ipg -query '{params.query_string}' | efetch -format fasta | gzip -c > {output}
        """

checkpoint split_seqs:
    input:
        pgenes = target_protein
    output:
        split_dir = temporary(directory("{o}/protein_seqs_split_files"))
    shell:
        """
        python scripts/split-fasta.py {input.pgenes} -o {wildcards.o} -n 1000
        """

rule blastn:
  input:
    genes = "{o}/protein_seqs_split_files/protein_seqs_{chunk}.fasta",
    db = "{o}/blast_dbs/cds_db.pdb"
  output:
     temp = temporary("{o}/blast_results/chunks/blast_{chunk}_processed.csv"),
  params:
     db = "{o}/blast_dbs/cds_db"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 18 * 60 * attempt,
    runtime = lambda wildcards, attempt: 18 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  threads: 32
  conda: "envs/sebastian_env.yml"
  shell:
    """
    blastp \
      -query {input.genes} \
      -db {params.db} \
      -evalue 1e-20 \
      -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
      -out {output.temp} \
      -num_threads {threads}

    """

def aggregate_chunks(wildcards):
    # Wait for the checkpoint to complete and get its output directory
    checkpoint_output = checkpoints.split_seqs.get(o=wildcards.o).output.split_dir
    
    # Identify all the generated {chunk} numbers (1, 2, 3, etc.)
    file_pattern = os.path.join(checkpoint_output, "protein_seqs_{chunk}.fasta")
    chunks = glob_wildcards(file_pattern).chunk
    
    # Tell Snakemake to expect a processed output file for every chunk found
    return expand("{o}/blast_results/chunks/blast_{chunk}_processed.csv", 
                  o=wildcards.o, 
                  chunk=chunks)

rule gather_blastp:
    input:
        temp = aggregate_chunks
    output:
        csv = "{o}/blast_results/blastp_output.csv",
    shell:
        """
        cat <(printf "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,e_value,bit_score,qcovs\\n") \
            {input.temp} > {output.csv}
        """

rule make_contig_list: 
  input:
    blast_hits = "{o}/blast_results/blastp_output.csv",
    cds = "{o}/cds_output/combined_cds.fna",
    script = "scripts/extract_cds_contigs.py",
  output:
    contigs = "{o}/cds_output/filtered_cds.fna",
  priority: 1
  resources:
    mem_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    time = lambda wildcards, attempt: 3 * 60 * attempt,
    runtime = lambda wildcards, attempt: 3 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    {input.script} {input.blast_hits} {input.cds} {output.contigs}
    """

rule minimap2:
  input:
    contigs = "{o}/cds_output/filtered_cds.fna",
    fq1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq.gz", #"{o}/fastq/{sample}_1.fastq.gz",
    fq2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq.gz", #"{o}/fastq/{sample}_2.fastq.gz",
  output:
    sorted_bam = temporary("{o}/mapped_reads/{sample}.sorted.bam"),
  threads: 16
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 18 * 60 * attempt,
    runtime = lambda wildcards, attempt: 18 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:"""
    minimap2 -ax sr -t {threads} \
      {input.contigs} \
      {input.fq1} {input.fq2} \
    | samtools sort -@ {threads} -o {output.sorted_bam}
    """
   
rule index_bam:
  input:
    sorted_bam = "{o}/mapped_reads/{sample}.sorted.bam"
  output:
    sorted_bai = temporary("{o}/mapped_reads/{sample}.sorted.bam.bai"),
  conda: "envs/sebastian_env.yml"
  resources:
    mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    time = lambda wildcards, attempt: 4 * 60 * attempt,
    runtime = lambda wildcards, attempt: 4 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  shell:
    """
    samtools index {input.sorted_bam}
    """

rule reads_mapped:
  input:
    sorted_bam = "{o}/mapped_reads/{sample}.sorted.bam.bai",
    bam = "{o}/mapped_reads/{sample}.sorted.bam",
    script = "scripts/aggregating_depth.py",
  output:
    map = "{o}/mapped_reads/stats/{sample}_mapping_stats.tsv",
    coverage = "{o}/mapped_reads/stats/{sample}_coverage_stats.tsv",
    ave_cov = "{o}/mapped_reads/stats/{sample}_average_coverage_stats.tsv",
  resources:
    mem_mb = lambda wildcards, attempt: 12 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 14 * 1024 * attempt,
    time = lambda wildcards, attempt: 4 * 60 * attempt,
    runtime = lambda wildcards, attempt: 4 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  priority: 1
  shell:
    """
    (echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
    samtools idxstats {input.bam} | \
    awk -F'\t' '{{print $0 "\t" ($3 - $4)}}') > {output.map}

    samtools depth -aa {input.bam} > {output.coverage}

    {input.script} {output.coverage} {output.ave_cov}
    """

rule merge_stats: 
  input:
    map = "{o}/mapped_reads/stats/{sample}_mapping_stats.tsv",
    ave_cov = "{o}/mapped_reads/stats/{sample}_average_coverage_stats.tsv",
    script = "scripts/merge_stats_3.py",
  output:
    merged = "{o}/results/filtered_{sample}_merged_stats.csv",
  params:
    contig_list_csv = "{o}/contiglist.csv",
    fna_dir = "{o}/prodigal_output/", 
    fna_ext = "fna"
  resources:
    mem_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    time = lambda wildcards, attempt: 1 * 60 * attempt,
    runtime = lambda wildcards, attempt: 1 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    {input.script} \
      {wildcards.sample} \
      -d $(dirname {input.map}) \
      -o {output.merged} \
      --filter
    """

rule merge_stats_wholistic: 
  input:
    meta = config['sample_metadata'],
    merged = expand("{o}/results/filtered_{sample}_merged_stats.csv", o=OUTDIR, sample = SAMPLES, symbol = SYMBOLS),
    script = "scripts/merge_stats_4.py",
  output:
    target_output
  params:
    whole = "{o}/all-mapped-coverage-stats.csv",
    target = f"-t '{protein_name}'" if protein_name else "",
  resources:
    mem_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 4 * 1024 * attempt,
    time = lambda wildcards, attempt: 1 * 60 * attempt,
    runtime = lambda wildcards, attempt: 1 * 60 * attempt,
    allowed_jobs = lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    {input.script} \
      -m {input.meta} \
      {params.target} \
      -i $(dirname {input.merged[0]}) \
      -o {output} 
    """
