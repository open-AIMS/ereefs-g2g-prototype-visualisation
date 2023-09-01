"""
Eric Lawrey - Australian Institute of Marine Science
This script generates a video for each year showing both the G2G land run off
data and the eReefs GBR1 Hydro salinity data. These videos are intended as a
prototype of this visualisation. The GBR1 Hydro data needed for this visualisations
can be obtained by running `02-get-daily-ereefs-hydro-data.py`. The basemap
data get be obtained by running `01-download-base-map-data.py`. The G2G test
data is not yet public. It was downloaded from the eReefs area on NCI.
"""
import sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm
import geopandas as gpd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import matplotlib.colors as colors
from datetime import datetime
import matplotlib
import os
from shapely.geometry import Point

# Variable to plot from G2G model
v='g2gflow'

export = 'export'

animate = True  # Don't animate when we are prototyping new changes to the plots

if not os.path.exists(export):
    os.makedirs(export)
        
for year in range(2016, 2023):
    print(f'Creating video for {year}')
    year_str = str(year)
    basemap_path = os.path.join('src-data','GBR_AIMS_eReefs-basemap')
    name = f"G2G-GBR1-Hydro-v2.0-Salt_{year_str}"
    video_file = os.path.join(export, f"{name}.mp4")
    video_file_tmp = os.path.join(export, f"{name}.tmp.mp4")
    
    # Support cancellations and restarts
    if os.path.exists(video_file):
        print(f'Skipping video for {year} as {video_file} already exists')
        continue
    
    # Clean up any temp files left over from previous runs
    if os.path.exists(video_file_tmp):
        print('Removing temporary file')
        os.remove(video_file_tmp)
    
    
    
    # Read the shapefile with geopandas
    rivers = gpd.read_file(os.path.join(basemap_path,'GBR_AIMS_eReefs-basemap_GA-topo5m-drainage.shp'))
    land = gpd.read_file(os.path.join(basemap_path,'GBR_AIMS_eReefs-basemap_Land-and-Basins.shp'))
    reefs = gpd.read_file(os.path.join(basemap_path,'GBR_AIMS_eReefs-basemap_Reefs.shp'))

    # =============== Load G2G data ===============
    # Load G2G data
    g2g_root = 'D:/eReefs/eReefs_data_from_BOM/G2G_gridded_test_2015_to_2022/test_gridded_nc'

    g2g_nc_path = f"{g2g_root}/{year_str}/test_grids_{v}_20*.nc"
    print(g2g_nc_path)
    
    g2g_ds = xr.open_mfdataset(g2g_nc_path)

    # Dig out var we want
    g2g_data_raw = g2g_ds["{}".format(v)]

    # Replace -999 vals with nan (I thought we already got rid of those?)
    g2g_data_nans = g2g_data_raw.where(g2g_data_raw!=-999.)
    # Accumulate to daily timescale (actually average...)
    g2g_data = g2g_data_nans.resample(time='1D', skipna=True).mean(skipna=True)
     
    # Set the extent of the G2G to the extent of the raster data
    lat_min, lat_max = g2g_data.lat.min().values, g2g_data.lat.max().values
    lon_min, lon_max = g2g_data.lon.min().values, g2g_data.lon.max().values
    
    # ============= Load Salinity ===============
    gbr1_salt_root = f'src-data\eReefs-hydro\GBR1_H2p0_salt_crop_{year_str}*.nc'
    gbr1_salt_ds = xr.open_mfdataset(gbr1_salt_root, combine='nested', concat_dim='time')

    gbr1_salt = gbr1_salt_ds['salt']
    
    # Set the extent and limits for gbr1_salt visualisation scale
    lat_min_salt, lat_max_salt = gbr1_salt.latitude.min().values, gbr1_salt.latitude.max().values
    lon_min_salt, lon_max_salt = gbr1_salt.longitude.min().values, gbr1_salt.longitude.max().values
    

    # ============= Plot setup =============
    fig, ax = plt.subplots(figsize=(12, 14.5), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    ax.set_aspect('equal')
    
    # ============= Set plot extents =========
    ax.set_extent([min(lon_min,lon_min_salt), max(lon_max,lon_max_salt), 
        min(lat_min,lat_min_salt), max(lat_max, lat_max_salt)], crs=ccrs.PlateCarree())

    extent = [lon_min, lon_max, lat_min, lat_max]

    # ============= Add Land, Rivers and Reefs ===========
    # Plot the rivers as a base map
    land.plot(ax=ax, linewidth=0.5, edgecolor='#898989', facecolor='#F7F7EB', zorder=1, transform=ccrs.PlateCarree())
    rivers.plot(ax=ax, linewidth=0.5, edgecolor='#9090d050', facecolor='none', zorder=1, transform=ccrs.PlateCarree())

    # Place the reefs on top (zorder=3)
    reefs.plot(ax=ax, color='black', zorder=3, alpha=0.4)
    # =========== Add cities ============
    cities = gpd.read_file(os.path.join(basemap_path,'GBR_AIMS_eReefs-basemap_Cities_2023.csv'))
    cities['latitude'] = cities['latitude'].astype(float)
    cities['longitude'] = cities['longitude'].astype(float)


    cities['geometry'] = cities.apply(lambda row: Point(row['longitude'], row['latitude']), axis=1)
    cities_gdf = gpd.GeoDataFrame(cities)
    # Convert ereefs_scale to numeric, replacing errors with np.nan
    cities_gdf['ereefs_scale'] = pd.to_numeric(cities_gdf['ereefs_scale'], errors='coerce')
    cities_gdf_filtered = cities_gdf[cities_gdf['ereefs_scale'] <= 2]
    # Limit the cities to those that are on the plot (otherwise it draws cities
    # outside the plot)
    minx, maxx, miny, maxy = ax.get_extent()
    cities_gdf_filtered = cities_gdf_filtered.cx[minx:maxx, miny:maxy]

    cities_gdf_filtered.plot(ax=ax, color='black', zorder=3, markersize=10)
    for x, y, label in zip(cities_gdf_filtered.geometry.x, cities_gdf_filtered.geometry.y, cities_gdf_filtered['name']):
        ax.text(x, y, label + ' ', verticalalignment='center', horizontalalignment='right', fontsize=14, 
                path_effects=[pe.withStroke(linewidth=3, foreground='white')])
                
    

    # ========== G2G plot ==============
    g2g_d0 = g2g_data[0].values

    # Set the limits on the visualisation scale
    vmin = 1
    vmax = 1e3

    # Define the custom color ramp
    color_ramp = np.array([
        '#f7fbff00',
        '#6baed6ff',
        '#08519cff',
        '#021e44ff'
    ])

    # Create a ListedColormap object using the color ramp
    #custom_cmap = ListedColormap(color_ramp)
    cmap = colors.LinearSegmentedColormap.from_list("", color_ramp)

    # Set transparent color for low values
    low_value_threshold = 1e-1
    colors_ramp = cmap(np.arange(cmap.N))
    colors_ramp[:int(low_value_threshold * cmap.N), -1] = 0
    transparent_cmap = mpl.colors.ListedColormap(colors_ramp)

    # Set the boundary for the colormap and create a BoundaryNorm object
    bounds = np.logspace(np.log10(vmin), np.log10(vmax), cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, transparent_cmap.N)
    
    im = plt.imshow(g2g_d0, cmap=transparent_cmap, norm=norm, extent = extent, zorder=2, transform=ccrs.PlateCarree())
    # ============== Plot Salinity ===============
    extent_salt = [lon_min_salt, lon_max_salt, lat_min_salt, lat_max_salt]
    vmin_salt = gbr1_salt.min().values
    vmax_salt = gbr1_salt.max().values

    # Based on https://github.com/eatlas/GBR_AIMS_eReefs-basemap/blob/main/colour-ramps/styles/RedBlueRainbowSalt_24-36-PSU.pal
    salt_color_ramp = [
        '#5c0035', '#6b0032', '#7a002f', '#89002d', '#98002a', '#a60027', '#b30224',
        '#ba0821', '#c10e1e', '#c8141a', '#cf1a17', '#d62014', '#dd2610', '#e42c0d',
        '#eb3209', '#f23806', '#f93e03', '#ff4500', '#ff4e00', '#ff5700', '#ff6000',
        '#ff6a00', '#ff7300', '#ff7c00', '#ff8500', '#ff9009', '#ff9b14', '#ffa620',
        '#ffb12b', '#ffbd37', '#ffc540', '#ffcb47', '#ffd14d', '#ffd754', '#ffdc5b',
        '#ffe262', '#ffe868', '#ffee6f', '#fff372', '#fff873', '#fffd74', '#f1ffa0',
        '#d4ffdb', '#9cf9dc', '#62e6da', '#39cae2', '#358cf6', '#3146f0', '#3a0db2',
        '#380060'
    ]

    salt_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", salt_color_ramp)
    salt_norm = matplotlib.colors.Normalize(vmin=24, vmax=36)


    # Plot gbr1_salt data. Use origin='lower' to flip the y-direction
    im_salt = plt.imshow(gbr1_salt[0].values, cmap=salt_cmap, vmin=vmin_salt, vmax=vmax_salt, extent=extent_salt, zorder=2, transform=ccrs.PlateCarree(), origin='lower')

    # ============== Legends ===============
    dates = g2g_data.time.values
    date_str = pd.to_datetime(str(dates[0])).strftime('%Y-%m-%d')

    plt.title(date_str, fontsize=24)

    cbar = plt.colorbar(im, orientation='horizontal', fraction=0.05, pad=0.05)
    cbar.ax.set_title('Daily flow (m^3/s)', fontsize=14)
    cbar_salt = plt.colorbar(im_salt, orientation='horizontal', fraction=0.05, pad=0.05)
    cbar_salt.ax.set_title('Salinity (PSU)', fontsize=14)


    # ============== Animate the plots ==============
    if animate:
        fps = 8
        N = min(len(g2g_data),len(gbr1_salt))

        with open("test.log", "w", buffering=1) as f:
            def animate_func(i):
                if i % fps ==0: 
                    message = f'Writing frame {i}/{N} for year {year}\n'
                    f.write(message)
                    sys.stdout.write(message)
                im.set_array(g2g_data[i].values)
                im_salt.set_array(gbr1_salt[i].values)
                date_str = pd.to_datetime(str(dates[i])).strftime('%Y-%m-%d')
                plt.title(date_str, fontsize=22)
                return [im, im_salt]
            
            anim = animation.FuncAnimation(fig, animate_func, frames=N, interval=1000/fps, blit=True)
            anim.save(video_file_tmp, writer=animation.FFMpegWriter(bitrate=5000))
            os.rename(video_file_tmp, video_file)

    print("Done!")
