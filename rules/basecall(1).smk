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
        "Rule {rule} started processing {input}."
    shell:
        "pod5 convert fast5 {input} -t 12 -r --output {output} 2> {log}"


rule dorado_basecall:
    # Produces a single BAM with move tables; --emit-moves is required for duplex pairing.
    input:
        os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}")
    output:
        directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"))
    log:
        os.path.join(LOGDIR, "dorado_basecall_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input}."
    shell:
        "bash scripts/run_dorado_basecaller.sh -i {input} -o {output} --pore {wildcards.pore} 2> {log}"


rule find_duplexes:
    # Identifies distant pairs (pair_ids_filtered.txt) and creates split POD5 files.
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
    # Two basecalling iterations, outputs kept separate:
    #   distant.fastq — duplex reads from different molecules
    #   split.fastq   — duplex reads from self-ligated molecules
    # These are NOT merged here; merging and subset assembly happens in classify_reads.
    input:
        pod5s       = os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}"),
        duplex_data = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}"),
        pair_ids    = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt"),
        simplex_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}")
    output:
        duplex_dir     = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}")),
        distant_fastq  = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/distant.fastq"),
        split_fastq    = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/split.fastq")
    params:
        threads = 4
    log:
        os.path.join(LOGDIR, "duplex_basecall_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input.duplex_data}."
    shell:
        """
        python scripts/run_duplex_basecall.py \
            --pod5_dir {input.pod5s} \
            --duplex_data {input.duplex_data} \
            --simplex_dir {input.simplex_dir} \
            -o {output.duplex_dir} \
            --pore {wildcards.pore} \
            --threads {params.threads} 2> {log}
        """


rule filter_split_reads:
    # Remove from the simplex BAM any reads that became ancestors of duplex pairs
    # (both distant and split), to prevent double-counting in downstream analysis.
    # BAM → FASTQ conversion happens here; output is FASTQ for all downstream rules.
    input:
        simplex_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"),
        duplex_data = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}")
    output:
        out_dir        = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-filtered/{pore}/{run}")),
        id_split_orig  = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/ids/id_split_orig.txt"),
        id_split_left  = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/ids/id_split_left.txt"),
        id_split_right = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/ids/id_split_right.txt")
    conda:
        "_envs/seqkit-env.yaml"
    log:
        os.path.join(LOGDIR, "filter_split_reads_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} removing duplex-ancestor reads from simplex set for {wildcards.sample}."
    shell:
        """
        split_pairs="{input.duplex_data}/split_duplex_pair_ids.txt"
        mkdir -p {output.out_dir} $(dirname {output.id_split_orig})

        awk '{{print $1}}' "$split_pairs" > {output.id_split_orig}
        awk '{{print $2}}' "$split_pairs" > {output.id_split_left}
        awk '{{print $3}}' "$split_pairs" > {output.id_split_right}

        # Collect all ancestor IDs: distant pair IDs + split pair original IDs
        distant_pairs="{input.duplex_data}/pair_ids_filtered.txt"
        cat <(awk '{{print $1}}' "$distant_pairs") \
            <(awk '{{print $2}}' "$distant_pairs") \
            {output.id_split_orig} \
            | sort -u > {output.out_dir}/../_all_ancestor_ids.txt

        for bam in {input.simplex_dir}/*.bam; do
            base=$(basename "$bam")
            samtools fastq "$bam" \
                | seqkit grep -v -f {output.out_dir}/../_all_ancestor_ids.txt \
                > {output.out_dir}/${{base%.bam}}.fastq
        done 2>> {log}
        """


rule trim_single_reads:
    input:
        directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-filtered/{pore}/{run}"))
    output:
        out_dir  = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}")),
        adapters = os.path.join(LOGDIR, "adapters_simplex_{sample}_{pore}_{run}.txt")
    params:
        threads = 8
    conda:
        "_envs/porechop-env.yml"
    log:
        os.path.join(LOGDIR, "trim_single_reads_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input}."
    shell:
        """
        bash scripts/trim_adapters.sh \
            -i {input} \
            -o {output.out_dir} \
            -s {output.adapters} \
            --threads {params.threads} 2> {log}
        """


rule assess_single_trimming:
    input:
        untrimmed = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-filtered/{pore}/{run}")),
        trimmed   = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}"))
    output:
        directory(os.path.join(OUTDIR, "{sample}/trimming-assessment/simplex/{pore}/{run}"))
    conda:
        "_envs/bioawk-env.yaml"
    log:
        os.path.join(LOGDIR, "assess_single_trimming_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} comparing trimmed vs untrimmed simplex reads for {wildcards.sample}."
    shell:
        "bash scripts/assess_trimming.sh --untrimmed {input.untrimmed} --trimmed {input.trimmed} --output {output} 2> {log}"


rule trim_duplex_reads:
    # Trims distant.fastq and split.fastq separately to preserve their identity.
    input:
        distant = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/distant.fastq"),
        split   = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/split.fastq")
    output:
        distant_dir      = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/distant")),
        split_dir        = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/split")),
        adapters_distant = os.path.join(LOGDIR, "adapters_duplex_distant_{sample}_{pore}_{run}.txt"),
        adapters_split   = os.path.join(LOGDIR, "adapters_duplex_split_{sample}_{pore}_{run}.txt")
    params:
        threads = 8
    conda:
        "_envs/porechop-env.yml"
    log:
        os.path.join(LOGDIR, "trim_duplex_reads_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input.distant} and {input.split}."
    shell:
        """
        bash scripts/trim_adapters.sh \
            -i {input.distant} \
            -o {output.distant_dir} \
            -s {output.adapters_distant} \
            --threads {params.threads} 2> {log}

        bash scripts/trim_adapters.sh \
            -i {input.split} \
            -o {output.split_dir} \
            -s {output.adapters_split} \
            --threads {params.threads} 2>> {log}
        """


rule assess_duplex_trimming:
    input:
        untrimmed_distant = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/distant.fastq"),
        untrimmed_split   = os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}/split.fastq"),
        trimmed_distant   = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/distant")),
        trimmed_split     = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/split"))
    output:
        distant_assessment = directory(os.path.join(OUTDIR, "{sample}/trimming-assessment/duplex-distant/{pore}/{run}")),
        split_assessment   = directory(os.path.join(OUTDIR, "{sample}/trimming-assessment/duplex-split/{pore}/{run}"))
    conda:
        "_envs/bioawk-env.yaml"
    log:
        os.path.join(LOGDIR, "assess_duplex_trimming_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} comparing trimmed vs untrimmed duplex reads for {wildcards.sample}."
    shell:
        """
        bash scripts/assess_trimming.sh \
            --untrimmed {input.untrimmed_distant} \
            --trimmed {input.trimmed_distant} \
            --output {output.distant_assessment} 2> {log}

        bash scripts/assess_trimming.sh \
            --untrimmed {input.untrimmed_split} \
            --trimmed {input.trimmed_split} \
            --output {output.split_assessment} 2>> {log}
        """


rule remove_calibration_strands:
    # Align reads to calibration reference (Lambda/PhiX).
    # Produces per-type alignment stats (flagstat + samtools stats), separate
    # calibration FASTQs per type, and a combined report TSV.
    # Use --keep in params if the reference organism shares sequences with the
    # calibration spike-in (e.g. rRNA locus comparison) and mapped reads should
    # be retained in the clean output.
    # Requires config["calib_ref"] pointing to the calibration FASTA.
    # Requires config["keep_calib"] = True/False (optional, default False).
    input:
        simplex         = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-trimmed/{pore}/{run}")),
        duplex_distant  = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/distant")),
        duplex_split    = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-trimmed/{pore}/{run}/split"))
    output:
        calib_dir              = directory(os.path.join(OUTDIR, "{sample}/calibration-strands/{pore}/{run}")),
        calib_report           = os.path.join(OUTDIR, "{sample}/calibration-strands/{pore}/{run}/quality_report.tsv"),
        simplex_clean          = directory(os.path.join(OUTDIR, "{sample}/basecalled-simplex-clean/{pore}/{run}")),
        duplex_distant_clean   = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}/distant")),
        duplex_split_clean     = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}/split"))
    params:
        threads    = 8,
        calib_ref  = config["calib_ref"],
        keep_flag  = lambda wc: "--keep" if config.get("keep_calib", False) else ""
    conda:
        "_envs/minimap2-env.yaml"
    log:
        os.path.join(LOGDIR, "remove_calibration_strands_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} finding and removing calibration strands for {wildcards.sample}."
    shell:
        "bash scripts/remove_calibration_strands.sh"
        " --simplex {input.simplex}"
        " --duplex-distant {input.duplex_distant}"
        " --duplex-split {input.duplex_split}"
        " --calib-out {output.calib_dir}"
        " --simplex-clean {output.simplex_clean}"
        " --duplex-distant-clean {output.duplex_distant_clean}"
        " --duplex-split-clean {output.duplex_split_clean}"
        " --report {output.calib_report}"
        " --ref {params.calib_ref}"
        " --threads {params.threads}"
        " {params.keep_flag}"
        " 2> {log}"


rule classify_reads:
    # Assemble the four canonical read subsets from clean data:
    #   simplex.fastq          — all simplex reads (ancestors removed)
    #   duplex_distant.fastq   — duplex reads from different molecules
    #   duplex_split.fastq     — duplex reads from self-ligated molecules
    #   simplex_no_ancestors.fastq — simplex reads with all duplex ancestors removed
    #                                (subset 4; same as simplex here since filter_split_reads
    #                                already removed both distant and split ancestors)
    # The first three subsets are primary; subset 4 is retained as an explicit file
    # for clarity and downstream QC comparison.
    input:
        simplex         = os.path.join(OUTDIR, "{sample}/basecalled-simplex-clean/{pore}/{run}"),
        duplex_distant  = os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}/distant"),
        duplex_split    = os.path.join(OUTDIR, "{sample}/basecalled-duplex-clean/{pore}/{run}/split")
    output:
        out_dir             = directory(os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}")),
        simplex_fq          = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/simplex.fastq"),
        duplex_distant_fq   = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/duplex_distant.fastq"),
        duplex_split_fq     = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/duplex_split.fastq"),
        simplex_no_anc_fq   = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/simplex_no_ancestors.fastq")
    conda:
        "_envs/bioawk-env.yaml"
    log:
        os.path.join(LOGDIR, "classify_reads_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} assembling read subsets for {wildcards.sample}."
    shell:
        """
        mkdir -p {output.out_dir}

        cat {input.simplex}/*.fastq        > {output.simplex_fq}
        cat {input.duplex_distant}/*.fastq > {output.duplex_distant_fq}
        cat {input.duplex_split}/*.fastq   > {output.duplex_split_fq}

        # Subset 4: simplex reads with both distant and split ancestors removed.
        # filter_split_reads already excluded all ancestors from the simplex set,
        # so this is a copy — retained as an explicit named file for downstream use.
        cp {output.simplex_fq} {output.simplex_no_anc_fq}

        echo "Read counts:" 2> {log}
        for f in {output.simplex_fq} {output.duplex_distant_fq} \
                 {output.duplex_split_fq} {output.simplex_no_anc_fq}; do
            n=$(bioawk -c fastx 'END{{print NR}}' "$f")
            echo "  $(basename $f): $n reads" >> {log}
        done
        """


rule merge_and_classify:
    # Per-run merge of the four clean subsets, then appended to total-merged/.
    input:
        simplex_fq        = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/simplex.fastq"),
        duplex_distant_fq = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/duplex_distant.fastq"),
        duplex_split_fq   = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/duplex_split.fastq"),
        simplex_no_anc_fq = os.path.join(OUTDIR, "{sample}/read-subsets/{pore}/{run}/simplex_no_ancestors.fastq")
    output:
        merged_dir       = directory(os.path.join(OUTDIR, "{sample}/merged/{pore}/{run}")),
        total_merged_dir = directory(os.path.join(OUTDIR, "{sample}/total-merged"))
    conda:
        "_envs/bioawk-env.yaml"
    log:
        os.path.join(LOGDIR, "merge_and_classify_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} merging read subsets for {wildcards.sample}."
    shell:
        """
        python scripts/merge_and_classify.py \
            --simplex {input.simplex_fq} \
            --duplex-distant {input.duplex_distant_fq} \
            --duplex-split {input.duplex_split_fq} \
            --simplex-no-ancestors {input.simplex_no_anc_fq} \
            --merged-dir {output.merged_dir} \
            --total {output.total_merged_dir} 2> {log}
        """
