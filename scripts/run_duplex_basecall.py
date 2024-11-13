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
import glob

console = Console()

model_map = {
    "R9.4.1": "dna_r9.4.1_450bps_sup.cfg",
    "R10.0": "",
    "R10.3": "",
    "R10.4.0": "",
    "R10.4.1": "dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
}

@click.command()
@click.option('--pod5_dir', 
              required=True, 
              type=click.Path(exists=True), 
              help='Directory containing pod5 files',
              )
@click.option('--duplex_data', 
              required=True, 
              type=click.Path(exists=True), 
              help='Directory containing duplex info',
              )
@click.option('--simplex_dir', 
              required=True, 
              type=click.Path(exists=True), 
              help='Directory for basecalled simplex read files',
              )
@click.option('--duplex_dir', 
              '-o', required=True, 
              type=click.Path(exists=True), 
              help='Root directory for output files',
              )
@click.option('--merged_dir', 
              required=True, 
              type=click.Path(exists=True), 
              help='Directory for merged data outputting',
              )
@click.option('--pore', 
              type=click.Choice(model_map.keys(), case_sensitive=False),
              required=True, 
              help='Pore type used for sequencing',
              )
@click.option('--threads', 
              '-t', required=False, 
              type=click.INT, 
              help='CPU threads number',
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


def run_dorado_duplex(input_dir, output_dir, pairs_file, model, threads):
    """Run Dorado duplex basecalling workflow."""
    console.print("[bold blue]Running Dorado duplex basecalling...[/bold blue]")
    
    cmd = [
        "dorado", "duplex",
        "--pairs-file", pairs_file,
        "--model", model,
        "--device", "cuda:0",
        "--threads", threads,
        input_dir,
        output_dir,
    ]
    
    run_command(cmd)

def run_guppy_duplex(input_dir, output_dir, pairs_file, model, threads):
    """Run Guppy duplex basecalling workflow."""
    console.print("[bold blue]Running Guppy duplex basecalling...[/bold blue]")
    
    cmd = [
        "guppy_basecaller_duplex",
        "-i", input_dir,
        "-s", output_dir,
        "-x", "cuda:0",
        "-c", model,
        "--threads", threads,
        "--duplex_pairing_mode", "from_pair_list",
        "--duplex_pairing_file", pairs_file
    ]
    
    run_command(cmd)

def remove_previous(in_file):
    """Remove fastq files that were merged before,
    in order to avoid duplicating via concatenation"""
    
    if Path.exists(in_file):
        os.remove(in_file)

    console.print(f"[yellow]Merged file {in_file} existed before was removed.[/yellow]")


def merge_outputs(source_dir, destination, destination_merged):
    """Merge basecalling outputs."""
    console.print("[yellow]Merging outputs...[/yellow]")
    
    remove_previous(destination)

    source_files = source_dir.glob('*/*.fastq')

    for source_file in source_files:
        console.print(f"concatenating {source_files}...")
        
        with open(destination, "a") as outfile:
            with open(source_file, "r") as infile:
                outfile.write(infile.read())

        with open(destination_merged, "a") as outfile:
            with open(source_file, "r") as infile:
                outfile.write(infile.read())

    console.print("[green]Outputs has been merged.[/green]")


def do_duplex_basecall(pod5_dir, duplex_data, simplex_dir, duplex_dir, merged_dir, pore, threads):
    """Duplex basecalling workflow for different pore types."""
    
    pod5_dir = Path(pod5_dir)
    duplex_data = Path(duplex_data)
    simplex_dir = Path(simplex_dir)
    duplex_dir = Path(duplex_dir)
    merged_dir = Path(merged_dir)
    
    if not pore:
        pore = "R9.4.1"
    else:
        pore = pore.capitalize()

    if not threads:
        threads = 4

    model = model_map[pore]
    
    # Common paths
    distant_out = duplex_dir / "distant"
    pairs_filtered = duplex_data / "pair_ids_filtered.txt"
    distant_merged = merged_dir / "distant_merged.fastq"
    all_merged = merged_dir / "all_merged.fastq"

    remove_previous(all_merged)

    if pore == "R10.4.1":
        # run dorado duplex pipeline
        run_dorado_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered),
            model=str(model),
        )

        # The code for fastq exporting from BAM should be added 
        # before merge_output() function
        #merge_outputs(
        #    source_dir='',
        #    destination='',
        #    destination_merged=''
        #)

    else:
        # run legacy guppy duplex pipeline
        #
        # Iteration 1
        # processing standalone duplex reads
        console.print("\n[bold green]Starting Iteration 1: Distant pairs[/bold green]")
        console.print("processing standalone distand reads...")

        if not pairs_filtered.exists():
            raise click.FileError(str(pairs_filtered), "Pairs file not found")

        run_guppy_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered),
            model=str(model),
            threads=threads,
        )

        merge_outputs(
            source_dir=distant_out,
            destination=distant_merged,
            destination_merged=all_merged
        )

        # Iteration 2
        # splitting and processing concatenated duplex reads
        console.print("\n[bold green]Starting Iteration 2: Split pairs[/bold green]")
        console.print("splitting and processing concatenated reads...")

        split_input = duplex_dir / "pod5s_splitduplex" # generated by duplex-tools
        split_out = duplex_dir / "split"
        split_pairs = duplex_data / "split_duplex_pair_ids.txt"
        split_merged = merged_dir / "split_merged.fastq"

        if not split_pairs.exists():
            raise click.FileError(str(split_pairs), "Split pairs file not found")

        run_guppy_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs),
            model=str(model),
            threads=threads,
        )

        merge_outputs(
            source_dir=distant_out,
            destination=split_merged,
            destination_merged=all_merged
        )

    console.print("[bold green]Duplex reads were basecalled successfully![/bold green]")
 
do_duplex_basecall(
    pod5_dir, 
    duplex_data, 
    simplex_dir, 
    duplex_dir, 
    merged_dir, 
    pore, 
    threads
)