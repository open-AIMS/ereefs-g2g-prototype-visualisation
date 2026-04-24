"""
Download pre-generated daily G2G NetCDF files as a ZIP and extract them locally.

The download is sourced from a Nextcloud public share and requests a ZIP for
only the selected year sub-folder.
"""

import argparse
import shutil
import tempfile
import urllib.request
import zipfile
import re
from pathlib import Path


BASE_URL_TEMPLATE = (
    "https://nextcloud.eatlas.org.au/s/LiRXpzLFBCWPf4f/download?path=%2Fdaily%2Fg2gflow-data%2F{year}"
)
EXCLUDED_YEARS = {}

DAILY_FILE_PATTERN = re.compile(
    r"^BOM_eReefs-g2gflow_daily_\d{4}-\d{2}-\d{2}\.nc$"
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Download daily G2G NetCDF ZIP for a year and extract NetCDF files "
            "to src-data/g2g-data/daily-aggregated/<year>."
        )
    )
    parser.add_argument("year", type=str, help="Year to download, for example 2019.")
    parser.add_argument(
        "--destination-root",
        type=str,
        default="src-data/g2g-data/daily-aggregated",
        help="Root output folder for extracted daily files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing destination year folder if present.",
    )
    return parser.parse_args()


def download_zip(url: str, zip_path: Path) -> None:
    """Download a ZIP file from URL to local path."""
    print(f"Downloading {url}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "ereefs-g2g-prototype-visualisation/1.0"},
    )
    with urllib.request.urlopen(request) as response, open(zip_path, "wb") as f:
        shutil.copyfileobj(response, f)


def extract_daily_files(zip_path: Path, destination_year_dir: Path) -> int:
    """Extract daily G2G NetCDF files from ZIP into destination year folder."""
    with tempfile.TemporaryDirectory(prefix="g2g-daily-unzip-") as tmp_dir:
        tmp_path = Path(tmp_dir)

        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(tmp_path)

        matches = sorted(
            path
            for path in tmp_path.rglob("*.nc")
            if DAILY_FILE_PATTERN.match(path.name)
        )
        if not matches:
            raise ValueError("No daily G2G NetCDF files found in downloaded ZIP")

        destination_year_dir.mkdir(parents=True, exist_ok=True)
        for src_file in matches:
            dst_file = destination_year_dir / src_file.name
            shutil.move(str(src_file), str(dst_file))

    return len(matches)


def main() -> None:
    """Download and extract one year of daily G2G files."""
    args = parse_args()

    year = args.year
    if year in EXCLUDED_YEARS:
        print(f"Skipping excluded year from publication: {year}")
        return

    url = BASE_URL_TEMPLATE.format(year=year)
    destination_year_dir = Path(args.destination_root) / year

    if destination_year_dir.exists() and not args.force:
        print(
            f"Destination already exists: {destination_year_dir}. "
            "Use --force to replace it."
        )
        return

    if destination_year_dir.exists() and args.force:
        print(f"Removing existing destination: {destination_year_dir}")
        shutil.rmtree(destination_year_dir)

    destination_year_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="g2g-daily-download-") as tmp_dir:
        zip_path = Path(tmp_dir) / f"g2g_daily_{year}.zip"
        download_zip(url, zip_path)

        if not zipfile.is_zipfile(zip_path):
            raise ValueError(
                "Downloaded file is not a ZIP archive. "
                "Please verify the URL and permissions."
            )

        count = extract_daily_files(zip_path, destination_year_dir)

    print(f"Extracted {count} files to {destination_year_dir}")


if __name__ == "__main__":
    main()
