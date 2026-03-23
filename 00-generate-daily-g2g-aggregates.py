"""
Generate one daily aggregate NetCDF per day from hourly G2G source files.

This is intended as a one-time preprocessing step before publishing daily G2G
files to public web storage.
"""

import argparse
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd
import xarray as xr


VAR_NAME = "g2gflow"
NODATA_VALUE = -999.0


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Create daily mean G2G NetCDF files from hourly source files. "
            "Outputs one file per day."
        )
    )
    parser.add_argument(
        "years",
        nargs="*",
        help=(
            "Optional years to process (for example: 2019 2020). "
            "Defaults to all years found."
        ),
    )
    parser.add_argument(
        "--start-date",
        type=str,
        default=None,
        help="Optional start date (YYYY-MM-DD).",
    )
    parser.add_argument(
        "--end-date",
        type=str,
        default=None,
        help="Optional end date (YYYY-MM-DD).",
    )
    parser.add_argument(
        "--input-root",
        type=str,
        default="src-data/g2g-data/extracted_files",
        help="Root folder containing yearly source files.",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default="src-data/g2g-data/daily-aggregated",
        help="Root folder where daily aggregate files are written.",
    )
    return parser.parse_args()


def discover_years(input_root: Path, requested_years: list[str]) -> list[str]:
    """Resolve which year folders should be processed."""
    if requested_years:
        return sorted(requested_years)

    years = [d.name for d in input_root.iterdir() if d.is_dir() and d.name.isdigit()]
    return sorted(years)


def build_file_pattern() -> re.Pattern[str]:
    """Return a regex for matching source filenames and extracting date."""
    return re.compile(
        rf"^sidb2netcdf_{re.escape(VAR_NAME)}_(\d{{4}}-\d{{2}}-\d{{2}})\.nc$"
    )


def iter_source_files(
    year_dir: Path,
    filename_pattern: re.Pattern[str],
    start_ts: pd.Timestamp | None,
    end_ts: pd.Timestamp | None,
) -> Iterable[tuple[Path, str]]:
    """Yield matching source files and their YYYY-MM-DD date string."""
    for source_file in sorted(year_dir.glob("*.nc")):
        match = filename_pattern.match(source_file.name)
        if match is None:
            continue

        date_str = match.group(1)
        date_ts = pd.Timestamp(date_str)
        if start_ts is not None and date_ts < start_ts:
            continue
        if end_ts is not None and date_ts > end_ts:
            continue
        yield source_file, date_str


def aggregate_day(
    source_file: Path,
    output_file: Path,
    date_str: str,
) -> None:
    """Aggregate one source file to a daily mean and write output atomically."""
    output_tmp = output_file.with_suffix(".tmp.nc")
    if output_tmp.exists():
        output_tmp.unlink()

    with xr.open_dataset(source_file) as ds:
        if VAR_NAME not in ds:
            raise ValueError(f"Variable '{VAR_NAME}' not found in {source_file}")

        data = ds[VAR_NAME].where(ds[VAR_NAME] != NODATA_VALUE)
        daily = data.mean(dim="time", skipna=True, keep_attrs=True)
        daily = daily.expand_dims(time=[pd.Timestamp(date_str)])
        output_ds = daily.to_dataset(name=VAR_NAME)

        # Copy all source global metadata, then append provenance fields.
        output_ds.attrs = dict(ds.attrs)
        output_ds.attrs["source_file"] = source_file.name
        output_ds.attrs["aggregation"] = "Daily mean from hourly values"
        output_ds.attrs["generated_by"] = "Australian Institute of Marine Science (AIMS)"
        output_ds.attrs["references"] = "https://doi.org/10.26274/0K5S-5E13"
        output_ds.attrs["doi"] = "10.26274/0K5S-5E13"
        output_ds.attrs["processing_description"] = (
            "Processed by the Australian Institute of Marine Science (AIMS): "
            "created from hourly G2G flow by masking no-data (-999.0) and "
            "computing daily mean over the time dimension."
        )

        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        existing_history = ds.attrs.get("history", "")
        history_line = (
            f"{timestamp}: Aggregated hourly {VAR_NAME} to one daily mean "
            "slice for this date."
        )
        if existing_history:
            output_ds.attrs["history"] = f"{existing_history}\n{history_line}"
        else:
            output_ds.attrs["history"] = history_line

        # Explicitly preserve coordinate and variable metadata from source.
        for coord_name in output_ds.coords:
            if coord_name in ds.coords:
                output_ds[coord_name].attrs = dict(ds[coord_name].attrs)
        output_ds[VAR_NAME].attrs = dict(ds[VAR_NAME].attrs)

        # Update comment to reflect daily aggregation (not hourly).
        output_ds.attrs["comment"] = (
            "This file contains daily mean values. Each timestamp represents "
            "one calendar date. The time coordinate has no sub-daily component."
        )

        output_ds.to_netcdf(output_tmp)

    os.replace(output_tmp, output_file)


def main() -> None:
    """Entry point for generating daily aggregate G2G files."""
    args = parse_args()

    input_root = Path(args.input_root)
    output_root = Path(args.output_root)
    start_ts = pd.Timestamp(args.start_date) if args.start_date else None
    end_ts = pd.Timestamp(args.end_date) if args.end_date else None

    if start_ts is not None and end_ts is not None and start_ts > end_ts:
        raise ValueError("--start-date must be less than or equal to --end-date")

    if not input_root.exists():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")

    years = discover_years(input_root, args.years)
    if not years:
        raise ValueError(f"No year folders found under {input_root}")

    filename_pattern = build_file_pattern()
    total_written = 0
    total_skipped = 0

    for year in years:
        year_dir = input_root / year
        if not year_dir.exists():
            print(f"Skipping missing year folder: {year_dir}")
            continue

        year_output_dir = output_root / year
        year_output_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nProcessing year {year}...")

        processed_in_year = 0
        for source_file, date_str in iter_source_files(
            year_dir,
            filename_pattern,
            start_ts,
            end_ts,
        ):
            output_file = year_output_dir / f"BOM_eReefs-{VAR_NAME}_daily_{date_str}.nc"

            if output_file.exists():
                total_skipped += 1
                print(f"  Skipping {date_str}: output exists")
                continue

            print(f"  Aggregating {date_str}")
            aggregate_day(
                source_file=source_file,
                output_file=output_file,
                date_str=date_str,
            )
            total_written += 1
            processed_in_year += 1

        print(f"Year {year} complete: wrote {processed_in_year} file(s)")

    print("\nDone.")
    print(f"Wrote: {total_written} file(s)")
    print(f"Skipped existing: {total_skipped} file(s)")


if __name__ == "__main__":
    main()
