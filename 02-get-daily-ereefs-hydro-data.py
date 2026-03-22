"""
Eric Lawrey, Marc Hammerton - Australian Institute of Marine Science
This script downloads surface (-2.35 m) Salinity (salt) from the GBR4 AIMS
eReefs THREDDS data service. The data is cropped to the specified bounding
box to reduce the total amount of data. This bounding box was chosen to 
align with the G2G test data corresponding to the southern region of the 
GBR.

This script is designed to cope with cancellation and resumption. If a
data file has already been downloaded then it is skipped on subsequent
runs of this script.

This takes about 1 sec per day to download and 6.2 MB per day.
"""
import sys
import xarray as xr
import pandas as pd
import os
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description="Download surface (-2.35 m) Salinity (salt) from the GBR4 AIMS eReefs "
                                             "THREDDS data service.")

# Define parameters (positional or optional)
parser.add_argument("year", type=str, help="The year for which to download the files.")
parser.add_argument(
    "--start-date",
    type=str,
    default=None,
    help="Optional start date in YYYY-MM-DD format. Defaults to <year>-01-01.",
)
parser.add_argument(
    "--end-date",
    type=str,
    default=None,
    help="Optional end date in YYYY-MM-DD format. Defaults to <year>-12-31.",
)

# Parse the arguments
args = parser.parse_args()

# specify the URL of the OpenDAP endpoint
url = 'https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v4/daily.nc'

# specify the start and end dates
start_date = args.start_date or f"{args.year}-01-01"
end_date = args.end_date or f"{args.year}-12-31"
print(
    f"Downloading daily eReefs Hydro data for {args.year} "
    f"from {start_date} to {end_date} ..."
)

if start_date > end_date:
    print("Start date must be less than or equal to end date.")
    sys.exit(1)

dates = pd.date_range(start=start_date, end=end_date, freq='D') + pd.Timedelta('14h')

# specify the bounding box
north = -10.65
south = -29.30
west = 141.8
east = 155.8

# folder to save the downloaded data to 
destination_folder = os.path.join("src-data", "eReefs-hydro")

# depth layer to download
depth = -1.5

# for reference (GBR4 v4) — verified against remote zc coordinate 2026-03-20
gbr4_depth_to_k_table = {
    -0.5: 16,
    -1.5: 15,
    -3.0: 14,
    -5.55: 13,
    -8.8: 12,
    -12.75: 11,
    -17.75: 10,
    -23.75: 9,
    -31.0: 8,
    -39.5: 7,
    -49.0: 6,
    -60.0: 5,
    -73.0: 4,
    -88.0: 3,
    -103.0: 2,
    -120.0: 1,
    -145.0: 0,
}

k = gbr4_depth_to_k_table.get(depth, None)
if k is None:
    print(f"Depth {depth} not found in the lookup table.")
    sys.exit(1)

if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        
# open the remote dataset
ds = xr.open_dataset(url)

# download the data one day at a time, select the data within the bounding box, and save the result to a new NetCDF file
for i, date in enumerate(dates):
    filename = os.path.join(destination_folder, f'GBR4_H2p0_salt_crop_{date.strftime("%Y-%m-%d")}.nc')
    filename_tmp = f'{filename}.tmp'
    
    # Clean up any temp files left over from previous runs
    if os.path.exists(filename_tmp):
        print('Removing temporary file')
        os.remove(filename_tmp)
        
    # Support cancellations and restarts
    if os.path.exists(filename):
        print(f'Skipping data for {date} ({i+1}/{len(dates)}) file already exists')
    else:
        print(f'Downloading data for {date} ({i+1}/{len(dates)})')
        salt = ds['salt'].sel(k=k, time=date)
        salt = salt.sel(latitude=slice(south, north), longitude=slice(west, east))
        
        # Use a temporary file then rename to make the script more robust to cancellation and restart
        salt.to_netcdf(filename_tmp)
        os.rename(filename_tmp, filename)

print('Download complete.')
