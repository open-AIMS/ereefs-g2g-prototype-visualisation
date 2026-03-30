#!/bin/bash
set -euo pipefail

# Required environment variables:
#   YEAR_START     (for example: 2011)
#   YEAR_END       (for example: 2023)
: "${YEAR_START:?Environment variable YEAR_START is required}"
: "${YEAR_END:?Environment variable YEAR_END is required}"

echo "Downloading base map data (Step 1)..."
python3 01-download-base-map-data.py

for YEAR in $(seq "$YEAR_START" "$YEAR_END"); do
    echo "Processing year $YEAR..."

    echo "Step 2: Downloading eReefs Hydro salinity for $YEAR"
    python3 02-get-daily-ereefs-hydro-data.py "$YEAR"

    echo "Step 3: Downloading daily G2G ZIP for $YEAR"
    python3 03-download-daily-g2g-data.py "$YEAR"

    echo "Step 4: Generating animations for $YEAR"
    python3 04-animate-G2G-and-salinity.py "$YEAR"

    echo "Cleaning downloaded data for $YEAR (keeping base map)"
    rm -f "src-data/eReefs-hydro/GBR4_H2p0_salt_crop_${YEAR}"*.nc
    rm -f "src-data/g2g-data/daily-aggregated/${YEAR}/"*.nc
    rmdir "src-data/g2g-data/daily-aggregated/${YEAR}" 2>/dev/null || true

    echo "Finished processing year $YEAR."
    echo
done

echo "All years processed."
