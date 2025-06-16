# main.py

from pathlib import Path
from typing import Annotated

import pandas as pd
import polars as pl
import typer

from . import core
from . import khipu

app = typer.Typer()


def load_dataframe(file_path: Path):
    """Loads a CSV or TSV file into a pandas or polars DataFrame."""
    if file_path.suffix == ".csv":
        try:
            return pl.read_csv(file_path)
        except Exception:
            return pd.read_csv(file_path)
    elif file_path.suffix == ".tsv":
        try:
            return pl.read_csv(file_path, separator="\\t")
        except Exception:
            return pd.read_csv(file_path, sep="\\t")
    else:
        raise ValueError("Unsupported file type. Please provide a .csv or .tsv file.")


@app.command()
def combine(
    ms1_path: Annotated[
        Path, typer.Argument(help="Path to the first mass spec file (.csv or .tsv)")
    ],
    ms2_path: Annotated[
        Path, typer.Argument(help="Path to the second mass spec file (.csv or .tsv)")
    ],
    output_path: Annotated[
        Path, typer.Option(help="Path to save the combined file.")
    ] = "combined.csv",
    mode: Annotated[str, typer.Option(help="Mode to run khipu in.")] = "pos",
    mz_tolerance_ppm: Annotated[float, typer.Option(help="M/Z tolerance in ppm.")] = 5,
    rt_tolerance: Annotated[float, typer.Option(help="RT tolerance in seconds.")] = 2,
    id_col: Annotated[int, typer.Option(help="Column index for the id column.")] = 0,
    rt_col: Annotated[int, typer.Option(help="Column index for the rt column.")] = 1,
    mz_col: Annotated[int, typer.Option(help="Column index for the mz column.")] = 2,
    int_cols: Annotated[
        tuple[int, int], typer.Option(help="Column indices for the intensity columns.")
    ] = (3, 4),
    delimiter: Annotated[str, typer.Option(help="Delimiter for the file.")] = ",",
):
    """
    Combines two mass spec files using the mass_combine logic.
    """
    print(f"Loading files: {ms1_path} and {ms2_path}")
    print("Running khipu...")

    ms1_df = khipu.khipu(
        ms1_path,
        id_col,
        rt_col,
        mz_col,
        int_cols,
        delimiter,
        mode,
        mz_tolerance_ppm,
        rt_tolerance,
    )
    ms2_df = khipu.khipu(
        ms2_path,
        id_col,
        rt_col,
        mz_col,
        int_cols,
        delimiter,
        mode,
        mz_tolerance_ppm,
        rt_tolerance,
    )

    print("Combining files...")
    combined_df = core.mass_combine(ms1_df, ms2_df)

    print(f"Saving combined file to: {output_path}")
    if isinstance(combined_df, pl.DataFrame):
        combined_df.write_csv(output_path)
    else:
        combined_df.to_csv(output_path, index=False)

    print("Done.")


if __name__ == "__main__":
    app()
