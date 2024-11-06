rule prepare_samples:
    input:
        "_data"
    output:
        "samples.yaml"
    log:
        os.path.join(LOGDIR, "prepare_samples.log")
    message:
        "Rule {rule} started processing {input}."
    shell:
        "[[ -d {input} ]] && python scripts/prepare_sampleslist.py"


rule do_checks:
    input:
        os.path.join(OUTDIR, "{sample}/basecalled/{pore}/{run}")
    output:
        "checkpoints/{sample}/{pore}/{run}"
    log:
        os.path.join(LOGDIR, "do_check_{sample}_{pore}_{run}.log")
    shell:
        "[[ -d {input} ]] && mkdir -p {output} && touch {output}/bc.success"