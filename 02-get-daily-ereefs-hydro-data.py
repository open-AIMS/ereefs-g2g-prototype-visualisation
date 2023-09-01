"""
Eric Lawrey - Australian Institute of Marine Science
This script downloads surface (-2.35 m) Salinity (salt) from the GBR1 AIMS
eReefs THREDDs data service. The data is cropped to the specified bounding
box to reduce the total amount of data. This bounding box was chosen to 
align with the G2G test data corresponding to the southern region of the 
GBR.

This script is designed to cope with cancellation and resumption. If a
data file has already been downloaded then it is skipped on subsequent
runs of this script.

This takes about 1 sec per day to download and 6.2 MB per day.
"""
import xarray as xr
import numpy as np
import pandas as pd
import os



# specify the URL of the OpenDAP endpoint
url = 'https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc'

# specify the start and end dates
start_date = '2016-01-01'
end_date = '2022-11-02'
dates = pd.date_range(start=start_date, end=end_date, freq='D') + pd.Timedelta('14h')

# specify the bounding box
north = -20.74
south = -28.31
west = 148.8
east = 154

# folder to save the downloaded data to 
destination_folder = os.path.join("src-data", "eReefs-hydro")

# depth layer to download
depth = -2.35

# for reference (GBR1)
gbr1_depth_to_k_table = {
    -0.5: 15,
    -2.35: 14,
    -5.25: 13,
    -9: 12,
    -13: 11,
    -18: 10,
    -24: 9,
    -31: 8,
    -39.5: 7,
    -49: 6,
    -60.0: 5,
    -73.0: 4,
    -88.0: 3,
    -103.0: 2,
    -120.0: 1,
    -140.0: 0
}

k = gbr1_depth_to_k_table.get(depth, None)
if k is None:
    print(f"Depth {depth_value} not found in the lookup table.")
    sys.exit(1)

if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        
# open the remote dataset
ds = xr.open_dataset(url)

# download the data one day at a time, select the data within the bounding box, and save the result to a new NetCDF file
for i, date in enumerate(dates):
    filename = os.path.join(destination_folder, f'GBR1_H2p0_salt_crop_{date.strftime("%Y-%m-%d")}.nc')
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
