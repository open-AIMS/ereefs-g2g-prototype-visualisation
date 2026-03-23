"""
Eric Lawrey, Marc Hammerton - Australian Institute of Marine Science
This script generates full-region and zoomed videos showing both G2G land runoff
and eReefs GBR4 Hydro salinity data. It combines the previous separate full and
zoom scripts into one entry point.
"""
import argparse
import os
import sys
from datetime import datetime
from typing import cast

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import geopandas as gpd
import matplotlib
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point

mpl.use("agg")

DATA_BOUNDS = {
    "north": -10.65,
    "south": -29.30,
    "west": 141.8,
    "east": 155.8,
}

ZOOM_REGIONS = {
    "queensland": {
        "bounding_box": {
            "north": -10.65,
            "south": -29.30,
            "west": 141.8,
            "east": 155.8,
        },
        "zoom_level": 1,
    },
    "north": {
        "bounding_box": {
            "north": -10.65,
            "south": -21.5,
            "west": 141.8,
            "east": 149.9,
        },
        "zoom_level": 2,
    },
    "central": {
        "bounding_box": {
            "north": -15.1,
            "south": -25.03,
            "west": 145.15,
            "east": 152.57,
        },
        "zoom_level": 2,
    },
    "south": {
        "bounding_box": {
            "north": -18.95,
            "south": -28.702,
            "west": 148.45,
            "east": 155.8,
        },
        "zoom_level": 2,
    },
}

FLOW_LINE_THICKNESS_BY_ZOOM_LEVEL = {
    1: 0.5,
    2: 0.05,
}


def normalize_coords(g2g: xr.DataArray, salt: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
    """Ensure coordinate axes are ascending for consistent geospatial plotting."""
    if g2g.lat.values[0] > g2g.lat.values[-1]:
        g2g = g2g.sortby("lat")
    if g2g.lon.values[0] > g2g.lon.values[-1]:
        g2g = g2g.sortby("lon")
    if salt.latitude.values[0] > salt.latitude.values[-1]:
        salt = salt.sortby("latitude")
    if salt.longitude.values[0] > salt.longitude.values[-1]:
        salt = salt.sortby("longitude")
    return g2g, salt


def _max_filter_with_nan(data: np.ndarray, radius: int) -> np.ndarray:
    """Apply a square max filter while preserving NaN as no-data."""
    if radius <= 0:
        return data

    rows, cols = data.shape
    padded = np.pad(data, radius, mode="constant", constant_values=np.nan)
    data_for_max = np.where(np.isfinite(padded), padded, -np.inf)

    shifted_views = []
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            row_start = radius + dy
            col_start = radius + dx
            shifted_views.append(
                data_for_max[row_start:row_start + rows, col_start:col_start + cols]
            )

    filtered = np.maximum.reduce(shifted_views)
    filtered[filtered == -np.inf] = np.nan
    return filtered


def _blend_with_nan(lower: np.ndarray, upper: np.ndarray, fraction: float) -> np.ndarray:
    """Blend two rasters while preserving no-data masks."""
    lower_filled = np.where(np.isfinite(lower), lower, 0.0)
    upper_filled = np.where(np.isfinite(upper), upper, 0.0)
    blended = lower_filled + (upper_filled - lower_filled) * fraction
    has_data = np.isfinite(lower) | np.isfinite(upper)
    return np.where(has_data, blended, np.nan)


def thicken_raster_lines(data: np.ndarray, radius: float = 1.0) -> np.ndarray:
    """Visually thicken raster lines; supports fractional positive radius values."""
    if radius <= 0:
        return data

    lower_radius = int(np.floor(radius))
    fraction = radius - lower_radius

    lower = _max_filter_with_nan(data, lower_radius)
    if fraction == 0:
        return lower

    upper = _max_filter_with_nan(data, lower_radius + 1)
    return _blend_with_nan(lower, upper, fraction)


def create_animation(
    region_name: str,
    year: str,
    map_bounds: dict,
    zoom_level: int,
    g2g_data: xr.DataArray,
    gbr4_salt: xr.DataArray,
    rivers: gpd.GeoDataFrame,
    land: gpd.GeoDataFrame,
    reefs: gpd.GeoDataFrame,
    cities_gdf: gpd.GeoDataFrame,
    transparent_cmap: colors.Colormap,
    salt_cmap: colors.Colormap,
    norm: colors.BoundaryNorm,
    ticks: list[float],
    export_dir: str,
    var_name: str,
    animate: bool,
) -> None:
    """Render an animation or preview image for one region."""
    print(f"\n{'=' * 60}")
    print(f"Generating {region_name} animation...")
    print(f"{'=' * 60}")

    north = map_bounds["north"]
    south = map_bounds["south"]
    west = map_bounds["west"]
    east = map_bounds["east"]
    print(f"Bounds: lat [{south}, {north}], lon [{west}, {east}]")

    # Map can extend beyond data bounds.
    map_north, map_south = north, south
    map_west, map_east = west, east

    # Data must remain within downloaded data bounds.
    data_north = min(north, DATA_BOUNDS["north"])
    data_south = max(south, DATA_BOUNDS["south"])
    data_west = max(west, DATA_BOUNDS["west"])
    data_east = min(east, DATA_BOUNDS["east"])

    year_str = str(year)
    name = f"G2G-GBR4-Hydro-v4-Salt_{year_str}-{region_name}"

    video_file = os.path.join(export_dir, f"{name}.mp4")
    video_file_tmp = os.path.join(export_dir, f"{name}.tmp.mp4")
    preview_file = os.path.join(export_dir, f"{name}-preview.png")

    if os.path.exists(video_file_tmp):
        os.remove(video_file_tmp)

    fig, ax = plt.subplots(figsize=(12, 16.5), subplot_kw={"projection": ccrs.PlateCarree()})
    ax = cast(GeoAxes, ax)
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.08, top=0.92)
    ax.set_aspect("equal")
    ax.set_extent([map_west, map_east, map_south, map_north], crs=ccrs.PlateCarree())

    land.plot(
        ax=ax,
        linewidth=0.5,
        edgecolor="#898989",
        facecolor="#F7F7EB",
        zorder=1,
        transform=ccrs.PlateCarree(),
    )
    rivers.plot(
        ax=ax,
        linewidth=0.5,
        edgecolor="#b4b4b450",
        facecolor="none",
        zorder=1,
        transform=ccrs.PlateCarree(),
    )
    reefs.plot(ax=ax, color="black", zorder=3, alpha=0.4)

    cities_gdf_filtered = cities_gdf[cities_gdf["ereefs_scale"] <= 1]
    minx, maxx, miny, maxy = ax.get_extent()
    cities_gdf_filtered = cities_gdf_filtered.cx[minx:maxx, miny:maxy]
    cities_gdf_filtered.plot(ax=ax, color="black", zorder=3, markersize=10)
    for x, y, label in zip(
        cities_gdf_filtered.geometry.x,
        cities_gdf_filtered.geometry.y,
        cities_gdf_filtered["name"],
    ):
        ax.text(
            x,
            y,
            label + " ",
            verticalalignment="center",
            horizontalalignment="right",
            fontsize=14,
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
        )

    g2g_zoom = g2g_data.sel(lat=slice(data_south, data_north), lon=slice(data_west, data_east))
    gbr4_salt_zoom = gbr4_salt.sel(
        latitude=slice(data_south, data_north),
        longitude=slice(data_west, data_east),
    )
    g2g_zoom, gbr4_salt_zoom = normalize_coords(g2g_zoom, gbr4_salt_zoom)

    if (
        g2g_zoom.sizes.get("lat", 0) == 0
        or g2g_zoom.sizes.get("lon", 0) == 0
        or gbr4_salt_zoom.sizes.get("latitude", 0) == 0
        or gbr4_salt_zoom.sizes.get("longitude", 0) == 0
    ):
        print(f"Skipping {region_name}: empty data slice for selected bounds")
        plt.close(fig)
        return

    extent_g2g = (
        float(g2g_zoom.lon.min().values),
        float(g2g_zoom.lon.max().values),
        float(g2g_zoom.lat.min().values),
        float(g2g_zoom.lat.max().values),
    )
    extent_salt = (
        float(gbr4_salt_zoom.longitude.min().values),
        float(gbr4_salt_zoom.longitude.max().values),
        float(gbr4_salt_zoom.latitude.min().values),
        float(gbr4_salt_zoom.latitude.max().values),
    )

    # Apply line thickening by zoom level.
    flow_line_thickness_px = FLOW_LINE_THICKNESS_BY_ZOOM_LEVEL.get(zoom_level, 0.5)
    print(f"Using flow line thickness {flow_line_thickness_px} px for zoom level {zoom_level}")
    river_frame_0 = thicken_raster_lines(g2g_zoom[0].values, radius=flow_line_thickness_px)

    im_river_flow = plt.imshow(
        river_frame_0,
        cmap=transparent_cmap,
        norm=norm,
        extent=extent_g2g,
        zorder=2,
        transform=ccrs.PlateCarree(),
        origin="lower",
    )

    vmin_salt = 24.0
    vmax_salt = 37.0
    im_salt = plt.imshow(
        gbr4_salt_zoom[0].values,
        cmap=salt_cmap,
        vmin=vmin_salt,
        vmax=vmax_salt,
        extent=extent_salt,
        zorder=0,
        transform=ccrs.PlateCarree(),
        origin="lower",
    )

    dates = g2g_zoom.time.values
    date_str = pd.to_datetime(str(dates[0])).strftime("%Y-%m-%d")
    ax_left = ax.get_position().x0
    ax_right = ax.get_position().x1
    region_label = region_name.capitalize()

    fig.text(
        0.5,
        0.978,
        "Daily River Flow and Salinity (-1.5 m)",
        ha="center",
        va="top",
        fontsize=24,
        fontweight="bold",
    )
    fig.text(ax_left, 0.944, region_label, ha="left", va="top", fontsize=20)
    date_text = ax.text(
        0.5,
        0.944,
        date_str,
        transform=fig.transFigure,
        ha="center",
        va="top",
        fontsize=20,
        fontweight="bold",
        clip_on=False,
        path_effects=[pe.withStroke(linewidth=3, foreground="white")],
    )

    fig.text(
        ax_right,
        0.063,
        f"Variable IDs: {var_name}, salt",
        ha="right",
        va="bottom",
        fontsize=11,
        color="gray",
    )

    today_str = datetime.now().strftime("%d-%b-%Y")
    metadata_text = (
        f"Map generation: AIMS {today_str}\n"
        f"Data: BOM G2G (DOI: 10.26274/0K5S-5E13), eReefs CSIRO GBR4 Hydrodynamic Model v4 "
        f"(AIMS daily aggregate product, DOI: 10.26274/Y74K-T032)\n"
        f"Licensing: Creative Commons Attribution 4.0 International "
        f"(https://creativecommons.org/licenses/by/4.0/)"
    )
    fig.text(
        ax_left,
        0.020,
        metadata_text,
        ha="left",
        va="bottom",
        fontsize=9,
        color="black",
        linespacing=1.5,
    )

    cb_ax1 = ax.inset_axes((0.04, 0.03, 0.030, 0.28))
    cb_ax2 = ax.inset_axes((0.18, 0.03, 0.030, 0.28))
    cb1 = fig.colorbar(im_river_flow, cax=cb_ax1, orientation="vertical", ticks=ticks)
    salt_ticks = list(range(24, 38, 2))
    cb2 = fig.colorbar(im_salt, cax=cb_ax2, orientation="vertical", ticks=salt_ticks)
    cb_ax1.set_facecolor((1, 1, 1, 0.6))
    cb_ax2.set_facecolor((1, 1, 1, 0.6))
    cb1.set_label(
        "Daily Mean River Flow (m^3/s)",
        fontsize=12,
        fontweight="bold",
        labelpad=2,
    )
    cb2.set_label(
        "Average Salinity (PSU)",
        fontsize=12,
        fontweight="bold",
        labelpad=8,
    )
    cb1.ax.yaxis.set_tick_params(labelsize=10)
    cb2.ax.yaxis.set_tick_params(labelsize=10)
    cb1.ax.yaxis.set_label_position("right")
    cb1.ax.yaxis.tick_right()
    cb1.minorticks_off()
    cb1.ax.set_yticklabels([f"{tick:,.0f}" for tick in ticks])

    if animate:
        fps = 8
        N = min(len(g2g_zoom), len(gbr4_salt_zoom))
        with open("test.log", "a", buffering=1) as file_obj:
            file_obj.write(f"\n{'=' * 60}\n{region_name.upper()}\n{'=' * 60}\n")

            def animate_func(i: int):
                if i % fps == 0:
                    message = f"Writing frame {i}/{N} for {region_name} ({year})\n"
                    file_obj.write(message)
                    sys.stdout.write(message)
                date_str_i = pd.to_datetime(str(dates[i])).strftime("%Y-%m-%d")
                date_text.set_text(date_str_i)
                river_frame_i = thicken_raster_lines(
                    g2g_zoom[i].values,
                    radius=flow_line_thickness_px,
                )
                im_river_flow.set_array(river_frame_i)
                im_salt.set_array(gbr4_salt_zoom[i].values)
                return [im_river_flow, im_salt, date_text]

            anim = animation.FuncAnimation(
                fig,
                animate_func,
                frames=N,
                interval=1000 / fps,
                blit=True,
            )
            anim.save(video_file_tmp, writer=animation.FFMpegWriter(bitrate=5000))
            os.rename(video_file_tmp, video_file)
            print(f"✓ Saved {video_file}")
    else:
        fig.savefig(preview_file, dpi=100, bbox_inches="tight")
        print(f"✓ Saved {preview_file}")

    plt.close(fig)


def main() -> None:
    """Entrypoint for generating full and/or zoomed animations."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate full and zoomed videos showing G2G land runoff and "
            "eReefs GBR4 Hydro salinity data."
        )
    )
    parser.add_argument("year", type=str, help="The year for which to generate videos.")
    parser.add_argument(
        "--regions",
        type=str,
        default="queensland,north,central,south",
        help=(
            "Comma-separated regions to generate. Supported: "
            "queensland,north,central,south"
        ),
    )
    parser.add_argument(
        "--preview-image",
        action="store_true",
        help="Save one preview PNG per region instead of MP4 animations.",
    )
    args = parser.parse_args()

    year = args.year
    animate = not args.preview_image
    regions_to_process = [region.strip() for region in args.regions.split(",") if region.strip()]

    print(f"Generating outputs for year {year}...")
    print(f"Regions: {', '.join(regions_to_process)}")

    var_name = "g2gflow"
    g2g_root = "./src-data/g2g-data/daily-aggregated"
    export_dir = "export"
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)

    year_str = str(year)
    basemap_path = os.path.join("src-data", "GBR_AIMS_eReefs-basemap")

    print("Loading basemap data...")
    rivers = gpd.read_file(
        os.path.join(basemap_path, "GBR_AIMS_eReefs-basemap_GA-topo5m-drainage.shp")
    )
    land = gpd.read_file(
        os.path.join(basemap_path, "GBR_AIMS_eReefs-basemap_Land-and-Basins.shp")
    )
    reefs = gpd.read_file(
        os.path.join(basemap_path, "GBR_AIMS_eReefs-basemap_Reefs.shp")
    )

    print("Loading G2G data...")
    g2g_nc_path = f"{g2g_root}/{year_str}/sidb2netcdf_{var_name}_daily_20*.nc"
    g2g_ds = xr.open_mfdataset(g2g_nc_path)
    g2g_data = g2g_ds[var_name]

    print("Loading salinity data...")
    gbr4_salt_root = f"src-data/eReefs-hydro/GBR4_H2p0_salt_crop_{year_str}*.nc"
    gbr4_salt_ds = xr.open_mfdataset(gbr4_salt_root, combine="nested", concat_dim="time")
    gbr4_salt = gbr4_salt_ds["salt"]

    g2g_data, gbr4_salt = normalize_coords(g2g_data, gbr4_salt)

    color_ramp = np.array(["#f7fbff00", "#4e95beff", "#03407eff", "#021e44ff"])
    cmap = colors.LinearSegmentedColormap.from_list("", color_ramp)
    colors_ramp = cmap(np.arange(cmap.N))
    colors_ramp[: int(1e-1 * cmap.N), -1] = 0
    transparent_cmap = colors.ListedColormap(colors_ramp)

    salt_color_ramp = [
        "#5c0035", "#7a002f", "#89002d", "#98002a", "#a60027", "#b30224", 
        "#ba0821", "#c8141a", "#cf1a17", "#d62014", "#dd2610", "#eb3209", 
        "#f93e03", "#ff4500", "#ff4e00", "#ff5700", "#ff6000", "#ff7c00", 
        "#ff8500", "#ff9009", "#ffa620", "#ffb12b", "#ffc540", "#ffcb47", 
        "#ffd14d", "#ffd754", "#ffdc5b", "#ffe262", "#ffee6f", "#fff372", 
        "#fff873", "#fffd74", "#f1ffa0", "#d4ffdb", "#9cf9dc", "#62e6da", 
        "#39cae2", "#358cf6", "#3146f0", "#3a0db2", "#3a128e", "#380060",
    ]

    salt_cmap = colors.LinearSegmentedColormap.from_list("", salt_color_ramp)

    vmin = 1.0
    vmax = 3000.0
    bounds = np.logspace(np.log10(vmin), np.log10(vmax), cmap.N)
    norm = colors.BoundaryNorm(bounds, transparent_cmap.N)

    ticks = []
    for decade in np.arange(np.floor(np.log10(vmin)), np.ceil(np.log10(vmax)) + 1):
        scale = 10 ** decade
        tick_value = 3 * scale
        if vmin <= tick_value <= vmax:
            ticks.append(tick_value)

    cities = gpd.read_file(os.path.join(basemap_path, "GBR_AIMS_eReefs-basemap_Cities_2023.csv"))
    cities["latitude"] = cities["latitude"].astype(float)
    cities["longitude"] = cities["longitude"].astype(float)
    cities["geometry"] = cities.apply(
        lambda row: Point(row["longitude"], row["latitude"]),
        axis=1,
    )
    cities_gdf = gpd.GeoDataFrame(cities)
    cities_gdf["ereefs_scale"] = pd.to_numeric(cities_gdf["ereefs_scale"], errors="coerce")

    for region_name in regions_to_process:
        if region_name in ZOOM_REGIONS:
            region_config = ZOOM_REGIONS[region_name]
            bounds_to_use = region_config["bounding_box"]
            zoom_level = region_config["zoom_level"]
        else:
            print(f"Warning: Unknown region '{region_name}'. Skipping.")
            continue

        create_animation(
            region_name=region_name,
            year=year,
            map_bounds=bounds_to_use,
            zoom_level=zoom_level,
            g2g_data=g2g_data,
            gbr4_salt=gbr4_salt,
            rivers=rivers,
            land=land,
            reefs=reefs,
            cities_gdf=cities_gdf,
            transparent_cmap=transparent_cmap,
            salt_cmap=salt_cmap,
            norm=norm,
            ticks=ticks,
            export_dir=export_dir,
            var_name=var_name,
            animate=animate,
        )

    print(f"\n{'=' * 60}")
    print("All requested outputs complete!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
