import os
import click
import rich.console as console
import subprocess


## Pseudocode block
#
# if pore == r10.4...
# do dorado duplex basecall workflow, check commands in github
#
# if pore == 9.4.1
# do guppy_basecaller_duplex
#
# ITERATION 1
#
# IN="${R}/pod5/"
# OUT="${R}/basecalled_duplex/distant"
# PAIRS="${R}/basecalled_duplex/pair_ids_filtered.txt"
#
# guppy_basecaller_duplex \
#    -i ${IN} \
#    -s ${OUT} \
#    -x 'cuda:0' \
#    -c dna_r9.4.1_450bps_sup.cfg \
#    --duplex_pairing_mode from_pair_list  \
#    --duplex_pairing_file ${PAIRS}
#
# do output merge
#
# ITERATION 2
# 
# IN="${DIR}/basecalled_duplex/pod5s_splitduplex/"
# OUT="${DIR}/basecalled_duplex/split"
# PAIRS="${R}/basecalled_duplex/split_duplex_pair_ids.txt"
#  
#  guppy_basecaller_duplex \
#    -i ${IN} \
#    -s ${OUT} \
#    -x 'cuda:0' \
#    -c dna_r9.4.1_450bps_sup.cfg \
#    --duplex_pairing_mode from_pair_list  \
#    --duplex_pairing_file ${PAIRS}
#
# do output merge

import os
import click
import subprocess
from rich.console import Console
from pathlib import Path
from enum import Enum

console = Console()

model_map = {
    "R9.4.1": "dna_r9.4.1_450bps_sup.cfg",
    "R10.0": "",
    "R10.3": "",
    "R10.4.0": "",
    "R10.4.1": "dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
}

@click.command()
@click.option('--root-dir', 
              '-r', required=True, 
              type=click.Path(exists=True), 
              help='Root directory containing input files',
              )
@click.option('--pore', 
              type=click.Choice(model_map.keys(), case_sensitive=False),
              required=True, 
              help='Pore type used for sequencing',
              )


def run_command(cmd, echo=True):
    """Run a shell command and handle errors."""
    if echo:
        console.print(f"[dim]Running command: {' '.join(cmd)}[/dim]")
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        return result
    except subprocess.CalledProcessError as e:
        console.print(f"[bold red]Error running command: {e}[/bold red]")
        console.print(f"[red]Error output: {e.stderr}[/red]")
        raise click.Abort()


def run_dorado_duplex(input_dir, output_dir, pairs_file, model):
    """Run Dorado duplex basecalling workflow."""
    console.print("[bold blue]Running Dorado duplex basecalling...[/bold blue]")
    
    cmd = [
        "dorado", "duplex",
        "--pairs-file", pairs_file,
        "--model", model,
        "--device", "cuda:0",
        input_dir,
        output_dir
    ]
    
    run_command(cmd)

def run_guppy_duplex(input_dir, output_dir, pairs_file, model):
    """Run Guppy duplex basecalling workflow."""
    console.print("[bold blue]Running Guppy duplex basecalling...[/bold blue]")
    
    cmd = [
        "guppy_basecaller_duplex",
        "-i", input_dir,
        "-s", output_dir,
        "-x", "cuda:0",
        "-c", model,
        "--duplex_pairing_mode", "from_pair_list",
        "--duplex_pairing_file", pairs_file
    ]
    
    run_command(cmd)



def do_duplex_basecall(root_dir, pore):
    """Duplex basecalling workflow for different pore types."""
    
    root_dir = Path(root_dir)
    if not pore:
        pore = "R9.4.1"

    model = model_map[pore]
    
    # Common paths
    pod5_dir = root_dir / "pod5"
    basecalled_dir = root_dir / "basecalled_duplex"
    distant_out = basecalled_dir / "distant"
    pairs_filtered = basecalled_dir / "pair_ids_filtered.txt"

    if pore == "R10.4.1":
        # run dorado duplex pipeline
        run_dorado_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered)
        )

    else:
        # run legacy guppy duplex pipeline
        #
        # Iteration 1
        # processing standalone duplex reads
        console.print("\n[bold green]Starting Iteration 1: Distant pairs[/bold green]")
        console.print("processing standalone distand reads...")

        run_guppy_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered)
        )

        # Iteration 1
        # splitting and processing concatenated duplex reads
        console.print("\n[bold green]Starting Iteration 2: Split pairs[/bold green]")
        console.print("splitting and processing concatenated reads...")

        run_guppy_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs)
        )

    console.print("[bold green]Duplex reads were basecalled successfully![/bold green]")

    # Ensure directories exist
    os.makedirs(basecalled_dir, exist_ok=True)
    
    # Iteration 1: Distant pairs
    console.print("\n[bold green]Starting Iteration 1: Distant pairs[/bold green]")
    

    
    if not pairs_filtered.exists():
        raise click.FileError(str(pairs_filtered), "Pairs file not found")

    else:  # R9.4.1
        
    
    merge_outputs(distant_out)
    
    # Iteration 2: Split pairs
    console.print("\n[bold green]Starting Iteration 2: Split pairs[/bold green]")
    
    split_input = basecalled_dir / "pod5s_splitduplex"
    split_out = basecalled_dir / "split"
    split_pairs = basecalled_dir / "split_duplex_pair_ids.txt"
    
    if not split_pairs.exists():
        raise click.FileError(str(split_pairs), "Split pairs file not found")
    
    if pore_type == PoreType.R10_4:
        run_dorado_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs)
        )
    else:  # R9.4.1
        run_guppy_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs)
        )
    
    merge_outputs(split_out)
    
    console.print("[bold green]Basecalling workflow completed successfully![/bold green]")

def merge_outputs(output_dir):
    """Merge basecalling outputs."""
    console.print("[yellow]Merging outputs...[/yellow]")
    # Add your output merging logic here
    pass

if __name__ == "__main__":
    main()