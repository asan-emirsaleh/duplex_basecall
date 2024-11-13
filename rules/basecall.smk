import os


rule convert_2_pod5:
    input:
        "_data/{sample}/fast5/{pore}/{run}"
    output:
        directory(os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}"))
    conda:
        "envs/pod5-env.yaml"
    log:
        os.path.join(LOGDIR, "convert_2_pod5_{sample}_{pore}_{run}.log")
    message:
        "Rule {rule} started processing {input}."
    shell:
        "pod5 convert fast5 {input} -t 12 -r --output {output}"


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
        "bash scripts/run_dorado_basecaller.sh -i {input} -o {output} --pore {wildcards.pore}"


rule find_dupleces:
    input:
        bam_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}"),
        pod5s = os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}")
    output:
        out_dir = directory(os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}")),
        pair_ids = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt")
    conda:
        "envs/duplex-tools-env.yaml"
    params:
        threads = 4
    shell:
        """
        python scripts/process_duplex.py \
            --out-dir {output.out_dir} \
            --in-dir {input.bam_dir} \
            --pod5s {input.pod5s} \
            --threads {params.threads}
        """


rule duplex_basecall:
    input:
        pod5s = os.path.join(OUTDIR, "{sample}/pod5/{pore}/{run}"),
        duplex_data = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}")
        pair_ids = os.path.join(OUTDIR, "{sample}/duplex-data/{pore}/{run}/split_duplex_pair_ids.txt")
        simplex_dir = os.path.join(OUTDIR, "{sample}/basecalled-simplex/{pore}/{run}")
    output:
        duplex_dir = directory(os.path.join(OUTDIR, "{sample}/basecalled-duplex/{pore}/{run}"))
        merged_dir = directory(os.path.join(OUTDIR, "{sample}/merged/{pore}/{run}"))
    params:
        threads = 4
    message:
        "Rule {rule} started processing {input.duplex_data}."
    log:
        os.path.join(LOGDIR, "dorado_basecall_{sample}_{pore}_{run}.log")
    shell:
        """
        bash scripts/do_duplex_basecall.py \
            --pod5_dir {input.pod5s} \
            --duplex_data {input.duplex_data} \
            --simplex_dir {input.simplex_dir} \
            -o {output.duplex_dir} \
            --merged_dir {output.merged_dir} \
            --pore {wildcards.pore} \
            --threads {params.threads}
        """