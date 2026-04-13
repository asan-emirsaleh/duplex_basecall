import os

# ---------------------------------------------------------------------------
# Path helpers — keeps rules concise and paths consistent throughout
# ---------------------------------------------------------------------------
def out(*args):
    return os.path.join(OUTDIR, *args)

def log(*args):
    return os.path.join(LOGDIR, *args)


# ---------------------------------------------------------------------------
# Format conversion
# ---------------------------------------------------------------------------
rule convert_2_pod5:
    input:  "_data/{sample}/fast5/{pore}/{run}"
    output: directory(out("{sample}/pod5/{pore}/{run}"))
    conda:  "_envs/pod5-env.yaml"
    log:    log("convert_2_pod5_{sample}_{pore}_{run}.log")
    message: "Rule {rule} processing {input}."
    shell:  "pod5 convert fast5 {input} -t 12 -r --output {output} 2> {log}"


# ---------------------------------------------------------------------------
# Simplex basecalling — produces BAM with move tables for duplex pairing
# ---------------------------------------------------------------------------
rule dorado_basecall:
    input:  out("{sample}/pod5/{pore}/{run}")
    output: directory(out("{sample}/basecalled-simplex/{pore}/{run}"))
    log:    log("dorado_basecall_{sample}_{pore}_{run}.log")
    message: "Rule {rule} processing {input}."
    shell:
        "bash scripts/run_dorado_basecaller.sh -i {input} -o {output} --pore {wildcards.pore} 2> {log}"


# ---------------------------------------------------------------------------
# Duplex pair identification — produces pair lists and split POD5s
# ---------------------------------------------------------------------------
rule find_duplexes:
    input:
        bam_dir = out("{sample}/basecalled-simplex/{pore}/{run}"),
        pod5s   = out("{sample}/pod5/{pore}/{run}")
    output:
        out_dir  = directory(out("{sample}/duplex-data/{pore}/{run}")),
        pair_ids = out("{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt")
    conda:   "_envs/duplex-tools-env.yaml"
    params:  threads = 4
    log:     log("find_duplexes_{sample}_{pore}_{run}.log")
    message: "Rule {rule} processing {input.bam_dir}."
    shell:
        """
        python scripts/find_duplex_candidates.py \
            --out-dir {output.out_dir} \
            --in-dir  {input.bam_dir} \
            --pod5s   {input.pod5s} \
            --threads {params.threads} 2> {log}
        """


# ---------------------------------------------------------------------------
# Duplex basecalling — two iterations, outputs kept separate by type
#   distant.fastq  duplex reads from different molecules
#   split.fastq    duplex reads from self-ligated molecules
# ---------------------------------------------------------------------------
rule duplex_basecall:
    input:
        pod5s       = out("{sample}/pod5/{pore}/{run}"),
        duplex_data = out("{sample}/duplex-data/{pore}/{run}"),
        pair_ids    = out("{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt"),
        simplex_dir = out("{sample}/basecalled-simplex/{pore}/{run}")
    output:
        duplex_dir    = directory(out("{sample}/basecalled-duplex/{pore}/{run}")),
        distant_fastq = out("{sample}/basecalled-duplex/{pore}/{run}/distant.fastq"),
        split_fastq   = out("{sample}/basecalled-duplex/{pore}/{run}/split.fastq")
    params:  threads = 4
    log:     log("duplex_basecall_{sample}_{pore}_{run}.log")
    message: "Rule {rule} processing {input.duplex_data}."
    shell:
        """
        python scripts/run_duplex_basecall.py \
            --pod5_dir    {input.pod5s} \
            --duplex_data {input.duplex_data} \
            --simplex_dir {input.simplex_dir} \
            -o            {output.duplex_dir} \
            --pore        {wildcards.pore} \
            --threads     {params.threads} 2> {log}
        """


# ---------------------------------------------------------------------------
# Filter simplex — remove ancestors of both duplex types, convert BAM → FASTQ
# ---------------------------------------------------------------------------
rule filter_simplex:
    input:
        simplex_dir = out("{sample}/basecalled-simplex/{pore}/{run}"),
        duplex_data = out("{sample}/duplex-data/{pore}/{run}")
    output:
        directory(out("{sample}/basecalled-simplex-filtered/{pore}/{run}"))
    conda:   "_envs/seqkit-env.yaml"
    log:     log("filter_simplex_{sample}_{pore}_{run}.log")
    message: "Rule {rule} removing duplex ancestors from simplex set."
    shell:
        """
        mkdir -p {output}
        ids=$(mktemp)

        # Collect ancestor IDs from both pair files
        awk '{{print $1; print $2}}' {input.duplex_data}/pair_ids_filtered.txt > "$ids"
        awk '{{print $1}}'           {input.duplex_data}/split_duplex_pair_ids.txt >> "$ids"
        sort -u -o "$ids" "$ids"

        for bam in {input.simplex_dir}/*.bam; do
            base=$(basename "${{bam%.bam}}")
            samtools fastq "$bam" \
                | seqkit grep -v -f "$ids" \
                > {output}/${{base}}.fastq
        done 2>> {log}

        rm "$ids"
        """


# ---------------------------------------------------------------------------
# Adapter trimming — one rule handles simplex and both duplex types via wildcard
#   read_type ∈ {simplex, duplex/distant, duplex/split}
# ---------------------------------------------------------------------------
rule trim_reads:
    input:
        out("{sample}/basecalled-{read_type}-filtered/{pore}/{run}")
            if "{read_type}" == "simplex"
            else out("{sample}/basecalled-duplex/{pore}/{run}/{read_type}.fastq")
    output:
        reads    = directory(out("{sample}/basecalled-{read_type}-trimmed/{pore}/{run}")),
        adapters = log("adapters_{read_type}_{sample}_{pore}_{run}.txt")
    params:  threads = 8
    conda:   "_envs/porechop-env.yml"
    log:     log("trim_{read_type}_{sample}_{pore}_{run}.log")
    message: "Rule {rule} trimming {wildcards.read_type} reads."
    shell:
        """
        bash scripts/trim_adapters.sh \
            -i {input} \
            -o {output.reads} \
            -s {output.adapters} \
            --threads {params.threads} 2> {log}
        """


# ---------------------------------------------------------------------------
# Trimming assessment — same wildcard pattern as trim_reads
# ---------------------------------------------------------------------------
rule assess_trimming:
    input:
        untrimmed = out("{sample}/basecalled-{read_type}-filtered/{pore}/{run}")
                    if "{read_type}" == "simplex"
                    else out("{sample}/basecalled-duplex/{pore}/{run}/{read_type}.fastq"),
        trimmed   = out("{sample}/basecalled-{read_type}-trimmed/{pore}/{run}")
    output:
        directory(out("{sample}/trimming-assessment/{read_type}/{pore}/{run}"))
    conda:   "_envs/bioawk-env.yaml"
    log:     log("assess_trimming_{read_type}_{sample}_{pore}_{run}.log")
    message: "Rule {rule} assessing {wildcards.read_type} trimming."
    shell:
        """
        bash scripts/assess_trimming.sh \
            --untrimmed {input.untrimmed} \
            --trimmed   {input.trimmed} \
            --output    {output} 2> {log}
        """


# ---------------------------------------------------------------------------
# Calibration strand removal — three typed inputs, three clean outputs
# config["calib_ref"]  path to calibration FASTA (Lambda / PhiX)
# config["keep_calib"] set True to retain calibration reads in clean output
# ---------------------------------------------------------------------------
rule remove_calibration_strands:
    input:
        simplex        = out("{sample}/basecalled-simplex-trimmed/{pore}/{run}"),
        duplex_distant = out("{sample}/basecalled-duplex/distant-trimmed/{pore}/{run}"),
        duplex_split   = out("{sample}/basecalled-duplex/split-trimmed/{pore}/{run}")
    output:
        calib_dir      = directory(out("{sample}/calibration-strands/{pore}/{run}")),
        calib_report   = out("{sample}/calibration-strands/{pore}/{run}/quality_report.tsv"),
        simplex_clean  = directory(out("{sample}/basecalled-simplex-clean/{pore}/{run}")),
        distant_clean  = directory(out("{sample}/basecalled-duplex/distant-clean/{pore}/{run}")),
        split_clean    = directory(out("{sample}/basecalled-duplex/split-clean/{pore}/{run}"))
    params:
        threads   = 8,
        calib_ref = config["calib_ref"],
        keep_flag = lambda wc: "--keep" if config.get("keep_calib", False) else ""
    conda:   "_envs/minimap2-env.yaml"
    log:     log("remove_calibration_strands_{sample}_{pore}_{run}.log")
    message: "Rule {rule} removing calibration strands for {wildcards.sample}."
    shell:
        """
        bash scripts/remove_calibration_strands.sh \
            --simplex              {input.simplex} \
            --duplex-distant       {input.duplex_distant} \
            --duplex-split         {input.duplex_split} \
            --calib-out            {output.calib_dir} \
            --simplex-clean        {output.simplex_clean} \
            --duplex-distant-clean {output.distant_clean} \
            --duplex-split-clean   {output.split_clean} \
            --report               {output.calib_report} \
            --ref                  {params.calib_ref} \
            --threads              {params.threads} \
            {params.keep_flag} 2> {log}
        """


# ---------------------------------------------------------------------------
# Read classification — assemble the four canonical subsets
#   simplex.fastq              all simplex reads (ancestors already removed)
#   duplex_distant.fastq       duplex from different molecules
#   duplex_split.fastq         duplex from self-ligated molecules
#   simplex_no_ancestors.fastq explicit copy of simplex for downstream clarity
# ---------------------------------------------------------------------------
rule classify_reads:
    input:
        simplex  = out("{sample}/basecalled-simplex-clean/{pore}/{run}"),
        distant  = out("{sample}/basecalled-duplex/distant-clean/{pore}/{run}"),
        split    = out("{sample}/basecalled-duplex/split-clean/{pore}/{run}")
    output:
        directory(out("{sample}/read-subsets/{pore}/{run}"))
    conda:   "_envs/bioawk-env.yaml"
    log:     log("classify_reads_{sample}_{pore}_{run}.log")
    message: "Rule {rule} assembling read subsets for {wildcards.sample}."
    shell:
        """
        mkdir -p {output}
        cat {input.simplex}/*.fastq > {output}/simplex.fastq
        cat {input.distant}/*.fastq > {output}/duplex_distant.fastq
        cat {input.split}/*.fastq   > {output}/duplex_split.fastq
        cp  {output}/simplex.fastq    {output}/simplex_no_ancestors.fastq

        # Log read counts per subset
        for f in {output}/*.fastq; do
            printf "%s\\t%d\\n" \
                "$(basename $f)" \
                "$(bioawk -c fastx 'END{{print NR}}' $f)"
        done > {log}
        """


# ---------------------------------------------------------------------------
# Merge and classify — per-run merge then cross-run total
# ---------------------------------------------------------------------------
rule merge_and_classify:
    input:
        out("{sample}/read-subsets/{pore}/{run}")
    output:
        merged_dir       = directory(out("{sample}/merged/{pore}/{run}")),
        total_merged_dir = directory(out("{sample}/total-merged"))
    conda:   "_envs/bioawk-env.yaml"
    log:     log("merge_and_classify_{sample}_{pore}_{run}.log")
    message: "Rule {rule} merging read subsets for {wildcards.sample}."
    shell:
        """
        python scripts/merge_and_classify.py \
            --subsets-dir {input} \
            --merged-dir  {output.merged_dir} \
            --total       {output.total_merged_dir} 2> {log}
        """
