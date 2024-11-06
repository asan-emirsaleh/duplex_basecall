import os
import gzip
import sys
import yaml

configfile: "config.yaml"

# Load the samples information
def get_samples():
    if os.path.exists("samples.yaml"):
        with open("samples.yaml") as f:
            return yaml.safe_load(f)
    return {}  # Return an empty dict if samples.yaml doesn't exist yet

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

#import scripts.prepare_sampleslist as prep
#prep.prepare_samples()

samples = get_samples()

# Flatten it into tuples of (sample, pore, run)
pores_runs = [(sample, pore, run) 
              for sample in samples
              for pore in samples[sample]['pores']
              for run in samples[sample]['pores'][pore]]


sample=[s for s, p, r in pores_runs]
pore=[p for s, p, r in pores_runs]
run=[r for s, p, r in pores_runs]


print(f'samples: {sample}\npores: {pore}\nruns: {run}')


rule all:
    input:
        "samples.yaml",
        expand(os.path.join(OUTDIR, "{sample}/basecalled/{pore}/{run}"), 
                sample=[s for s, p, r in pores_runs],
                pore=[p for s, p, r in pores_runs],
                run=[r for s, p, r in pores_runs])


include: "rules/basecall.smk"
