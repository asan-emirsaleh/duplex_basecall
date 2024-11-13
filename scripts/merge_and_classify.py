import os
import click
import subprocess
from rich.console import Console
from pathlib import Path

console = Console()


@click.command()
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
@click.option('--total', 
              required=True, 
              type=click.Path(exists=True), 
              help='Directory for merged data outputting',
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


def remove_previous(in_file):
    """Remove fastq files that were merged before,
    in order to avoid duplicating via concatenation"""
    
    if Path.exists(in_file):
        os.remove(in_file)

    console.print(f"[yellow]Merged file {in_file} existed before was removed.[/yellow]")


def merge_outputs(source_dir: Path, destination: Path, destination_merged: Path, suffix: str):
    """Merge basecalling outputs."""
    console.print("[yellow]Merging outputs...[/yellow]")
    
    remove_previous(destination)

    source_files = source_dir.glob('*/*.fastq')

    for source_file in source_files:
        console.print(f"concatenating {source_file}...")
        
        with open(destination, "a") as outfile:
            with open(source_file, "r") as infile:
                outfile.write(infile.read())

        with open(destination_merged, "a") as outfile:
            with open(source_file, "r") as infile:
                outfile.write(infile.read())

    console.print("[green]Outputs has been merged.[/green]")


def parse_id(fastq):
    """Parse FASTQ id with awk."""

    cmd = f'awk '
    subprocess.call(cmd, shell=True)

def do_duplex_basecall(pod5_dir, duplex_data, simplex_dir, duplex_dir, merged_dir, pore, threads):
    """Duplex basecalling workflow for different pore types."""
    
    duplex_data = Path(duplex_data)
    simplex_dir = Path(simplex_dir)
    duplex_dir = Path(duplex_dir)
    merged_dir = Path(merged_dir)

    
    # Common paths
    distant_out = duplex_dir / "distant"
    pairs_filtered = duplex_data / "pair_ids_filtered.txt"
    distant_merged = merged_dir / "duplex-distant_merged.fastq"
    all_merged = merged_dir / "duplex-all_merged.fastq"

    remove_previous(all_merged)

    merge_outputs(
            source_dir=distant_out,
            destination=distant_merged,
            destination_merged=all_merged
        )

    # Iteration 2
    # splitting and processing concatenated duplex reads
    console.print("\n[bold green]Starting Iteration 2: Split pairs[/bold green]")
    console.print("splitting and processing concatenated reads...")

    split_input = duplex_data / "pod5s_splitduplex" # generated by duplex-tools.py
    # Note that in run_duplex-tools_dorado.sh, the destination category 
    # is `duplex_calls_dorado/pod5s_splitduplex`, but not `duplex_data/pod5s_splitduplex`

    split_out = duplex_dir / "split"
    split_pairs = duplex_data / "split_duplex_pair_ids.txt"
    split_merged = merged_dir / "duplex-split_merged.fastq"

    if not split_pairs.exists():
        raise click.FileError(str(split_pairs), "Split pairs file not found")


    merge_outputs(
            source_dir=split_out,
            destination=split_merged,
            destination_merged=all_merged
        )

    console.print("[bold green]Duplex reads were basecalled successfully![/bold green]")

if __name__ == "__main__":
    do_duplex_basecall()