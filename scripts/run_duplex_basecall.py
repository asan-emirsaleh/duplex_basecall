import os
import click
import subprocess
from rich.console import Console
from pathlib import Path

console = Console()

model_map = {
    "R9.4.1":  "dna_r9.4.1_450bps_sup.cfg",
    "R10.0":   "",
    "R10.3":   "",
    "R10.4.0": "",
    "R10.4.1": "dna_r10.4.1_e8.2_400bps_sup@v4.2.0",
}

DORADO = "_apps/dorado/bin/dorado"


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


def remove_previous(in_file: Path):
    """Remove a merged file if it exists, to avoid duplicating on re-run."""
    if in_file.exists():
        os.remove(in_file)
        console.print(f"[yellow]Removed previous merged file: {in_file}[/yellow]")


def merge_outputs(source_dir: Path, destination: Path, destination_merged: Path):
    """Merge per-batch FASTQ outputs from a guppy output directory into
    a per-iteration file (destination) and the combined all-merged file."""
    console.print("[yellow]Merging outputs...[/yellow]")

    remove_previous(destination)

    source_files = list(source_dir.glob("*/*.fastq"))
    if not source_files:
        console.print(f"[red]No FASTQ files found in {source_dir}[/red]")
        raise click.Abort()

    for source_file in source_files:
        console.print(f"concatenating {source_file}...")
        with open(destination, "a") as outfile, open(source_file) as infile:
            outfile.write(infile.read())
        with open(destination_merged, "a") as outfile, open(source_file) as infile:
            outfile.write(infile.read())

    console.print("[green]Outputs merged.[/green]")


def run_dorado_duplex(input_dir, output_dir, pairs_file, model, threads):
    """Run Dorado duplex basecalling and export FASTQ from the resulting BAM."""
    console.print("[bold blue]Running Dorado duplex basecalling...[/bold blue]")

    bam_out = Path(output_dir) / "duplex.bam"
    cmd = [
        DORADO, "duplex",
        "--pairs", pairs_file,
        "--device", "cuda:all",
        "--threads", str(threads),
        model,
        input_dir,
    ]
    result = run_command(cmd)

    # Write BAM output (dorado writes to stdout by default)
    bam_out.parent.mkdir(parents=True, exist_ok=True)
    with open(bam_out, "wb") as f:
        f.write(result.stdout.encode() if isinstance(result.stdout, str) else result.stdout)

    # Convert BAM to FASTQ
    fastq_out = Path(output_dir) / "duplex.fastq"
    run_command(["samtools", "fastq", "-@", str(threads), str(bam_out), "-o", str(fastq_out)])

    return fastq_out


def run_guppy_duplex(input_dir, output_dir, pairs_file, model, threads):
    """Run Guppy duplex basecalling."""
    console.print("[bold blue]Running Guppy duplex basecalling...[/bold blue]")

    cmd = [
        "guppy_basecaller_duplex",
        "-i", input_dir,
        "-s", output_dir,
        "-x", "cuda:0",
        "-c", model,
        "--threads", str(threads),
        "--detect_adapter", "--detect_primer",
        "--detect_mid_strand_adapter",
        "--trim_strategy", "dna",
        "--trim_adapters", "--trim_primers",
        "--do_read_splitting",
        "--disable_pings",
        "--duplex_pairing_mode", "from_pair_list",
        "--duplex_pairing_file", pairs_file,
    ]
    run_command(cmd)


@click.command()
@click.option("--pod5_dir",
              required=True, type=click.Path(exists=True),
              help="Directory containing POD5 files")
@click.option("--duplex_data",
              required=True, type=click.Path(exists=True),
              help="Directory containing duplex pair info (duplex_tools output)")
@click.option("--simplex_dir",
              required=True, type=click.Path(exists=True),
              help="Directory with basecalled simplex BAM files")
@click.option("-o", "--duplex_dir",
              required=True, type=click.Path(),
              help="Output directory for duplex basecalling results")
@click.option("--pore",
              required=True, type=click.Choice(list(model_map.keys()), case_sensitive=False),
              help="Pore type")
@click.option("--threads", "-t",
              required=False, type=click.INT, default=4,
              help="CPU thread count")
def do_duplex_basecall(pod5_dir, duplex_data, simplex_dir, duplex_dir, pore, threads):
    """Duplex basecalling workflow: two iterations (distant + split), merged into a single FASTQ."""

    pod5_dir    = Path(pod5_dir)
    duplex_data = Path(duplex_data)
    simplex_dir = Path(simplex_dir)
    duplex_dir  = Path(duplex_dir)
    model = model_map[pore.upper()]
    if not model:
        raise click.BadParameter(f"No model configured for pore type {pore_key}. Please add it to model_map.")

    duplex_dir.mkdir(parents=True, exist_ok=True)

    # Output paths — all inside duplex_dir, matching SMK expectation of merged.fastq
    distant_out    = duplex_dir / "distant"
    split_out      = duplex_dir / "split"
    distant_merged = duplex_dir / "duplex-distant_merged.fastq"
    split_merged   = duplex_dir / "duplex-split_merged.fastq"
    all_merged     = duplex_dir / "merged.fastq"   # primary output consumed by SMK

    pairs_filtered = duplex_data / "pair_ids_filtered.txt"
    split_pairs    = duplex_data / "split_duplex_pair_ids.txt"
    split_input    = duplex_data / "pod5s_splitduplex"

    remove_previous(all_merged)

    if pore_key == "R10.4.1":
        # Dorado duplex pipeline
        # Iteration 1: distant pairs
        console.print("\n[bold green]Starting Iteration 1: Distant pairs (dorado)[/bold green]")
        if not pairs_filtered.exists():
            raise click.FileError(str(pairs_filtered), "Distant pairs file not found")

        distant_fastq = run_dorado_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered),
            model=model,
            threads=threads,
        )
        with open(all_merged, "a") as out, open(distant_fastq) as src:
            out.write(src.read())

        # Iteration 2: split pairs
        console.print("\n[bold green]Starting Iteration 2: Split pairs (dorado)[/bold green]")
        if not split_pairs.exists():
            raise click.FileError(str(split_pairs), "Split pairs file not found")

        split_fastq = run_dorado_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs),
            model=model,
            threads=threads,
        )
        with open(all_merged, "a") as out, open(split_fastq) as src:
            out.write(src.read())

    else:
        # Legacy Guppy duplex pipeline
        # Iteration 1: distant pairs
        console.print("\n[bold green]Starting Iteration 1: Distant pairs (guppy)[/bold green]")
        if not pairs_filtered.exists():
            raise click.FileError(str(pairs_filtered), "Distant pairs file not found")

        run_guppy_duplex(
            input_dir=str(pod5_dir),
            output_dir=str(distant_out),
            pairs_file=str(pairs_filtered),
            model=model,
            threads=threads,
        )
        merge_outputs(source_dir=distant_out, destination=distant_merged, destination_merged=all_merged)

        # Iteration 2: split pairs
        console.print("\n[bold green]Starting Iteration 2: Split pairs (guppy)[/bold green]")
        if not split_pairs.exists():
            raise click.FileError(str(split_pairs), "Split pairs file not found")

        run_guppy_duplex(
            input_dir=str(split_input),
            output_dir=str(split_out),
            pairs_file=str(split_pairs),
            model=model,
            threads=threads,
        )
        merge_outputs(source_dir=split_out, destination=split_merged, destination_merged=all_merged)

    console.print(f"[bold green]Duplex basecalling complete. Merged output: {all_merged}[/bold green]")


if __name__ == "__main__":
    do_duplex_basecall()