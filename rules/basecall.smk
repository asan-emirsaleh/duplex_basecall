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

