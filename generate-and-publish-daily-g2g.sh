#!/bin/bash
set -euo pipefail

# Required environment variables:
#   SOURCE_S3_PREFIX   (for example: s3://ereefs-g2g-data)
#   DEST_S3_PREFIX     (for example: s3://bucket/path)
#   SOURCE_AWS_PROFILE (AWS CLI profile with read access to source bucket)
#   DEST_AWS_PROFILE   (AWS CLI profile with write access to destination prefix)
: "${SOURCE_S3_PREFIX:?Environment variable SOURCE_S3_PREFIX is required}"
: "${DEST_S3_PREFIX:?Environment variable DEST_S3_PREFIX is required}"
: "${SOURCE_AWS_PROFILE:?Environment variable SOURCE_AWS_PROFILE is required}"
: "${DEST_AWS_PROFILE:?Environment variable DEST_AWS_PROFILE is required}"


# Normalize S3 prefixes to avoid accidental double slashes
SOURCE_S3_PREFIX="${SOURCE_S3_PREFIX%/}"
DEST_S3_PREFIX="${DEST_S3_PREFIX%/}"

TEMP_DIR="./src-data/g2g-data/temp_s3_files"
EXTRACTED_ROOT="./src-data/g2g-data/extracted_files"
DAILY_OUTPUT_ROOT="./src-data/g2g-data/daily-aggregated"

# If "true", keep extracted/daily files after upload.
KEEP_LOCAL="${KEEP_LOCAL:-false}"

if [[ "$#" -eq 0 ]]; then
  YEARS=($(seq 2011 2023))
else
  YEARS=("$@")
fi

log() {
  echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] $*"
}

mkdir -p "$TEMP_DIR" "$EXTRACTED_ROOT" "$DAILY_OUTPUT_ROOT"

for YEAR in "${YEARS[@]}"; do
  log "Processing year $YEAR"

  rm -f "$TEMP_DIR"/*.tar.gz || true

  log "Downloading source archives from $SOURCE_S3_PREFIX for year $YEAR"
  aws --profile "$SOURCE_AWS_PROFILE" s3 cp "$SOURCE_S3_PREFIX/" "$TEMP_DIR/" \
    --recursive --exclude "*" --include "*$YEAR*.tar.gz"

  shopt -s nullglob
  archives=("$TEMP_DIR"/*.tar.gz)
  shopt -u nullglob
  if [[ "${#archives[@]}" -eq 0 ]]; then
    log "No archives found for $YEAR. Skipping."
    continue
  fi

  log "Extracting ${#archives[@]} archive(s)"
  for file in "${archives[@]}"; do
    tar -xzf "$file" -C "$EXTRACTED_ROOT"
  done

  log "Generating daily NetCDF aggregates for $YEAR"
  python3 00-generate-daily-g2g-aggregates.py "$YEAR" \
    --input-root "$EXTRACTED_ROOT" \
    --output-root "$DAILY_OUTPUT_ROOT"

  if [[ ! -d "$DAILY_OUTPUT_ROOT/$YEAR" ]]; then
    log "Expected output folder missing: $DAILY_OUTPUT_ROOT/$YEAR"
    exit 1
  fi

  shopt -s nullglob
  daily_files=("$DAILY_OUTPUT_ROOT/$YEAR"/*.nc)
  shopt -u nullglob
  if [[ "${#daily_files[@]}" -eq 0 ]]; then
    log "No daily NetCDF files were generated for $YEAR. Aborting upload."
    exit 1
  fi

  log "Uploading daily files for $YEAR to $DEST_S3_PREFIX/$YEAR/"
  aws --profile "$DEST_AWS_PROFILE" s3 cp "$DAILY_OUTPUT_ROOT/$YEAR/" \
    "$DEST_S3_PREFIX/$YEAR/" --recursive

  if [[ "$KEEP_LOCAL" != "true" ]]; then
    log "Cleaning local artifacts for $YEAR"
    rm -rf "$EXTRACTED_ROOT/$YEAR" "$DAILY_OUTPUT_ROOT/$YEAR"
    rm -f "$TEMP_DIR"/*.tar.gz || true
  fi

  log "Finished year $YEAR"
  echo

done

log "All requested years processed."
