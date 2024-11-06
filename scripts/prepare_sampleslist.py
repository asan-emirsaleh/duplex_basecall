import os
import yaml
from rich.console import Console


DATA_DIR = "_data"
OUTPUT_FILE = "samples.yaml"

console = Console()

def prepare_samples():
    # Initialize an empty dictionary to store sample data
    samples = {}

    # Walk through the `_data` directory to gather samples, pores, and runs
    # assuming that directory structure is "{sample}/{pore}/{run}/["fast5_fail", "fast5_pass"]"
    for root, dirs, files in os.walk(DATA_DIR):
        # Look for fast5 directories in the path
        if "fast5" in root:
            # Split the path to get sample, pore, and run information
            path_parts = root.split(os.sep)
            if len(path_parts) > 4:
                sample = path_parts[1]   # Extract sample
                pore = path_parts[3]     # Extract pore
                run = path_parts[4]      # Extract run

                # Initialize sample if not added
                if sample not in samples:
                    samples[sample] = {"pores": {}}

                if pore not in samples[sample]["pores"]:
                    samples[sample]["pores"] = {pore: []}

                # Append the run to the sample's run list
                if run not in samples[sample]["pores"][pore]:
                    samples[sample]["pores"][pore].append(run)

    # Save the output to a YAML file
    with open(OUTPUT_FILE, 'w') as outfile:
        yaml.dump(samples, outfile, default_flow_style=False)

    console.print(f"[bold green]Samples information saved to {OUTPUT_FILE}[/bold green]")

prepare_samples()
