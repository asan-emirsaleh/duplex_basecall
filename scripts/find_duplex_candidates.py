from rich.console import Console
import os
import subprocess
import click
import glob
from pathlib import Path
import shutil

console = Console()

def check_file_exists(file_path):
    """Check if a file exists and raise an error if it doesn't."""
    if not os.path.exists(file_path):
        raise click.FileError(file_path, hint='File does not exist')
    return file_path

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

@click.command()
@click.option('--out-dir', required=True, type=click.Path(), help='Output directory for results')
@click.option('--in-dir', required=True, type=click.Path(exists=True), help='Input directory containing BAM files')
@click.option('--pod5s', required=True, type=click.Path(exists=True), help='Directory containing POD5 files')
@click.option('--threads', '-t', default=1, help='Number of threads to use')


def process_duplex_reads(out_dir, in_dir, pod5s, threads):
    """Process duplex reads using duplex_tools."""
    
    console.print("[bold yellow]Searching for duplex reads...[/bold yellow]")
    
    # Create output directory
    os.makedirs(out_dir, exist_ok=True)
    
    # Find BAM files
    bam_files = list(glob.glob(os.path.join(in_dir, "*.bam")))
    if not bam_files:
        raise click.BadParameter(f"No BAM files found in {in_dir}")
    elif len(bam_files) > 1:
        raise click.BadParameter(f"Multiple BAM files found in {in_dir}. Please provide a directory with only one BAM file.")
    else:
        bam_file = bam_files[0]  # Take the first BAM file if multiple exist
    
    console.print(f"Processing {bam_file}...")
    
    # Looking for pairs
    console.print("Looking for pairs...")
    run_command([
        "duplex_tools", "pair",
        "--threads", str(threads),
        "--output_dir", out_dir,
        bam_file
    ])
    
    # Splitting pairs
    console.print("Splitting pairs...")
    split_output_dir = os.path.join(out_dir, "pod5s_splitduplex")
    
    # Remove previous split_output_dir if it exists
    if os.path.exists(split_output_dir):
        shutil.rmtree(split_output_dir)
    
    run_command([
        "duplex_tools", "split_pairs",
        "--threads", str(threads),
        bam_file,
        pod5s,
        split_output_dir + "/"
    ])
    
    # Concatenate pair IDs into new file
    pair_files = glob.glob(os.path.join(split_output_dir, "*_pair_ids.txt"))
    with open(os.path.join(out_dir, "split_duplex_pair_ids.txt"), 'w') as outfile:
        for pair_file in pair_files:
            with open(pair_file) as infile:
                outfile.write(infile.read())
    
    console.print("[bold green]Processing complete![/bold green]")


process_duplex_reads()