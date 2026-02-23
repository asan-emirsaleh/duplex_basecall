import os


rule convert_2_pod5:
    input:
        "_data/{sample}/fast5/{pore}/{run}"
    output:
        directory(os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}"))
    conda:
        "_envs/pod5-env.yaml"
    log:
        os.path.join(LOGDIR, "convert_2_pod5_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input}"
    shell:
        "pod5 convert fast5 {input} -t 12 -r --output {output} 2> {log}"


rule dorado_basecall:
    input:
        os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}")
    output:
        directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"))
    message:
        "Rule {rule} started processing {input}."
    log:
        os.path.join(LOGDIR, "dorado_basecall_{sample}_{pore}_{run}.log")
    shell:
        "bash scripts/run_dorado_basecaller.sh -i {input} -o {output} --pore {wildcards.pore} 2> {log}"


rule find_duplexes:
    input:
        bam_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"),
        pod5s   = os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}")
    output:
        out_dir  = directory(os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}")),
        pair_ids = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt")
    conda:
        "_envs/duplex-tools-env.yaml"
    params:
        threads = 4
    log:
        os.path.join(LOGDIR, "find_duplexes_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input.bam_dir}."
    shell:
        """
        python scripts/find_duplex_candidates.py \
            --out-dir {output.out_dir} \
            --in-dir {input.bam_dir} \
            --pod5s {input.pod5s} \
            --threads {params.threads} 2> {log}
        """


rule duplex_basecall:
    input:
        pod5s       = os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}"),
        duplex_data = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}"),
        pair_ids    = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt"),
        simplex_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}")
    output:
        duplex_dir = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}")),
        merged_dir = directory(os.path.join(OUTDIR, "{sample}/merged/{pore}/{run}"))
    params:
        threads = 4
    message:
        "Rule {rule} started processing {input.duplex_data}."
    log:
        os.path.join(LOGDIR, "duplex_basecall_{sample}_{pore}_{run}.log")
    shell:
        """
        python scripts/run_duplex_basecall.py \
            --pod5_dir {input.pod5s} \
            --duplex_data {input.duplex_data} \
            --simplex_dir {input.simplex_dir} \
            -o {output.duplex_dir} \
            --merged_dir {output.merged_dir} \
            --pore {wildcards.pore} \
            --threads {params.threads} 2> {log}
        """


rule trim_single_reads:
    input:
        directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"))
    output:
        out_dir  = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}")),
        adapters = os.path.join(LOGDIR, "adapters_simplex_{sample}_{pore}_{run}.txt")
    params:
        threads = 8
    conda:
        "_envs/porechop-env.yml"
    message:
        "Rule {rule} started processing {input}."
    log:
        os.path.join(LOGDIR, "trim_single_reads_{sample}_{pore}_{run}.log")
    shell:
        """
        mkdir -p {output.out_dir}
        bash scripts/trim_adapters.sh \
            -i {output.out_dir} \
            -o {output.out_dir} \
            -s {output.adapters} \
            --threads {params.threads} 2> {log}
        """


rule assess_single_trimming:
    input:
        untrimmed = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}")),
        trimmed   = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}"))
    output:
        directory(os.path.join(OUTDIR, "{sample}/trimming-assessment/simplex/{pore}/{run}"))
    conda:
        "_envs/bioawk-env.yaml"
    message:
        "Rule {rule} comparing trimmed vs untrimmed simplex reads for {wildcards.sample}."
    log:
        os.path.join(LOGDIR, "assess_single_trimming_{sample}_{pore}_{run}.log")
    shell:
        "bash scripts/assess_trimming.sh --untrimmed {input.untrimmed} --trimmed {input.trimmed} --output {output} 2> {log}"


rule trim_duplex_reads:
    input:
        directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}"))
    output:
        out_dir  = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}")),
        adapters = os.path.join(LOGDIR, "adapters_duplex_{sample}_{pore}_{run}.txt")
    params:
        threads = 8
    conda:
        "_envs/porechop-env.yml"
    message:
        "Rule {rule} started processing {input}."
    log:
        os.path.join(LOGDIR, "trim_duplex_reads_{sample}_{pore}_{run}.log")
    shell:
        """
        mkdir -p {output.out_dir}
        for bam in {input}/*.bam; do
            base=$(basename "$bam" .bam)
            samtools fastq "$bam" > {output.out_dir}/${{base}}.fastq
        done
        bash scripts/trim_adapters.sh \
            -i {output.out_dir} \
            -o {output.out_dir} \
            -s {output.adapters} \
            --threads {params.threads} 2> {log}
        """


rule assess_duplex_trimming:
    input:
        untrimmed = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}")),
        trimmed   = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}"))
    output:
        directory(os.path.join(OUTDIR, "{sample}/trimming-assessment/duplex/{pore}/{run}"))
    conda:
        "_envs/bioawk-env.yaml"
    message:
        "Rule {rule} comparing trimmed vs untrimmed duplex reads for {wildcards.sample}."
    log:
        os.path.join(LOGDIR, "assess_duplex_trimming_{sample}_{pore}_{run}.log")
    shell:
        "bash scripts/assess_trimming.sh --untrimmed {input.untrimmed} --trimmed {input.trimmed} --output {output} 2> {log}"


rule remove_calibration_strands:
    # Calibration strands (Lambda phage or PhiX) are spiked-in sequencing controls.
    # Reads are aligned to the calibration reference; hits are moved to a dedicated
    # directory and per-run quality metrics are computed on them.
    # Requires config["calib_ref"] pointing to the calibration FASTA.
    input:
        simplex_dir = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}")),
        duplex_dir  = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}"))
    output:
        calib_dir     = directory(os.path.join(OUTDIR, "{sample}/calibration-strands/{pore}/{run}")),
        calib_report  = os.path.join(OUTDIR, "{sample}/calibration-strands/{pore}/{run}/quality_report.tsv"),
        simplex_clean = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-clean/{pore}/{run}")),
        duplex_clean  = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}"))
    params:
        threads   = 8,
        calib_ref = config["calib_ref"]
    conda:
        "_envs/minimap2-env.yaml"
    message:
        "Rule {rule} finding and removing calibration strands for {wildcards.sample}."
    log:
        os.path.join(LOGDIR, "remove_calibration_strands_{sample}_{pore}_{run}.log")
    shell:
        "bash scripts/remove_calibration_strands.sh"
        " --simplex {input.simplex_dir}"
        " --duplex {input.duplex_dir}"
        " --calib-out {output.calib_dir}"
        " --simplex-clean {output.simplex_clean}"
        " --duplex-clean {output.duplex_clean}"
        " --report {output.calib_report}"
        " --ref {params.calib_ref}"
        " --threads {params.threads}"
        " 2> {log}"


rule merge_and_classify:
    input:
        simplex_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex-clean/{pore}/{run}"),
        duplex_dir  = os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}")
    output:
        merged_dir       = directory(os.path.join(OUTDIR, "{sample}/merged/{pore}/{run}")),
        total_merged_dir = directory(os.path.join(OUTDIR, "{sample}/total-merged"))
    conda:
        "_envs/bioawk-env.yaml"
    message:
        "Rule {rule} started processing basecalled reads."
    log:
        os.path.join(LOGDIR, "merge_and_classify_{sample}_{pore}_{run}.log")
    shell:
        """
        python scripts/merge_and_classify.py \
            --simplex_dir {input.simplex_dir} \
            --duplex_dir {input.duplex_dir} \
            --merged_dir {output.merged_dir} \
            --total {output.total_merged_dir} 2> {log}
        """

