#%% PLOT REDOUBT

# This script loads in a Redoubt catalog of choice and plots the hypocenters in 3 panels
# (a N-S cross section, an E-W cross section, and a top-down view of the hypocenters)
# The user will have to specify and tweak lat and lon limits, and region settings
# For detailed description of plotting syntax, refer to PyGMT documentation

# Import all dependencies
import numpy as np
import pandas as pd
import pygmt
from matplotlib import pyplot as plt
from toolbox import raster2array, reader
from obspy import UTCDateTime, Catalog
from pandas import read_csv

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
PEC_events_filepath = main_dir + 'data/avo/davidof_december2021_catalog_FI.xml'
# relocated_catalog_filepath = main_dir + 'output/davidof/relocate_catalog/relocated_catalog.xml'
plot_output_dir = main_dir + 'output/davidof/relocate_catalog/'
VOLC_LAT = 51.97
VOLC_LON = 178.33
MAX_DEPTH = 10 # TRY 40  # km
REGION = [VOLC_LON-0.3, VOLC_LON+0.3, VOLC_LAT-0.15, VOLC_LAT+0.1475]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_WE = [VOLC_LON-0.3, VOLC_LON+0.3, -2.5, MAX_DEPTH]
REGION_NS = [-2.5, MAX_DEPTH, VOLC_LAT-0.15, VOLC_LAT+0.1475]
size_by_magnitude = True
plot_temporal = True
plot_FI = True
migration_track = False
swarm_focus = False

#%% Read all catalogs

# (1) Pre-existing catalog
PEC_events = reader(PEC_events_filepath)

# (3) GrowClust relocated catalog
# relocated_catalog = reader(relocated_catalog_filepath)  # min_cc = 0.75
#%% Extract information from catalog, filtered by max depth

# Choose catalog
catalog = PEC_events
if catalog == PEC_events:
    plot_title = '+t\"Original Catalog (N=%d)\"' % len(catalog)
    save_filepath1 = plot_output_dir + 'original_catalog_N' + str(len(PEC_events))
elif catalog == relocated_catalog:
    plot_title = '+t\"Relocated Catalog (N=%d)\"' % len(catalog)
    save_filepath1 = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog))

# Extract hypocenter, time, magnitude and FI information
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
times = np.array([np.datetime64(event.origins[0].time) for event in catalog])
magnitudes = np.array([event.magnitudes[0].mag for event in catalog])
if plot_FI:
    FI_values = []
    for event in catalog:
        if event.comments[-1].text.split('=')[1] == 'None':
            FI_values.append(np.nan)
        else:
            FI_values.append(float(event.comments[-1].text.split('=')[1]))
    FI_values = np.array(FI_values)

# Filter all arrays by depth
valid_index = np.where(depths<MAX_DEPTH)
latitudes = latitudes[valid_index]
longitudes = longitudes[valid_index]
depths = depths[valid_index]
times = times[valid_index]
magnitudes = magnitudes[valid_index]
if plot_FI:
    FI_values = FI_values[valid_index]

# Normalize the magnitudes based on reasonable size values:
if size_by_magnitude:
    # construct reference sizes
    ref_magnitudes = np.array(range(0,501))
    min_mag_marker = 0
    ref_sizes = []
    for magnitude in list(ref_magnitudes):
        normalized_magnitude = (magnitude - np.min(min_mag_marker))/np.ptp(ref_magnitudes[ref_magnitudes > min_mag_marker])
        ref_size = 0.06 + normalized_magnitude*0.24
        ref_sizes.append(ref_size)
    ref_magnitudes = list(ref_magnitudes)
    # now calculate sizes
    sizes = []
    for magnitude in magnitudes:
        if magnitude < min_mag_marker:
            size = ref_sizes[0]
        else:
            size = ref_sizes[ref_magnitudes.index(int(magnitude*100))]
        sizes.append(size)
    sizes = np.array(sizes)
    # define size for mag markers
    size0 = ref_sizes[ref_magnitudes.index(0)]
    size1 = ref_sizes[ref_magnitudes.index(100)]
    size2 = ref_sizes[ref_magnitudes.index(200)]
    size3 = ref_sizes[ref_magnitudes.index(300)]
    size4 = ref_sizes[ref_magnitudes.index(400)]

# # Prepare cross-section data for SW-NE cross section
# W = [VOLC_LAT, -153.83]
# E = [VOLC_LAT, -153.03]
# WE_x = [W[1], E[1]]
# WE_xvec = np.linspace(W[1], E[1], 200)
# WE_y = [W[0], E[0]]
# WE_yvec = np.linspace(W[0], E[0], 200)
#
# # Prepare cross-section data for NW-SE cross section
# N = [59.558, VOLC_LON]
# S = [59.163, VOLC_LON]
# NS_x = [N[1], S[1]]
# NS_xvec = np.linspace(N[1], S[1],600)
# NS_y = [N[0], S[0]]
# NS_yvec = np.linspace(N[0], S[0],600)
# # Load EW and NS elevation profiles (from https://apps.nationalmap.gov/elevation/)
# WE = read_csv(WE_profile_filename)
# NS = read_csv(NS_profile_filename)

# # Calculate down-section distance for all catalog events and elevation profile
# WE_down_distances = calc_down_distances(latitudes,longitudes,WE_yvec,WE_xvec)
# WE_elev_distances = calc_down_distances(WE["Lat"],WE["Lon"],WE_yvec,WE_xvec)
# NS_down_distances = calc_down_distances(latitudes,longitudes,NS_yvec,NS_xvec)
# NS_elev_distances = calc_down_distances(NS["Lat"],NS["Lon"],NS_yvec,NS_xvec)
# WE_down_distances = np.array(WE_down_distances)
# NS_down_distances = np.array(NS_down_distances)

# import geopy.distance
# WE_down_distances = [geopy.distance.GeodesicDistance(W,[lat,lon]).km for lat,lon in zip(list(np.ones(len(latitudes))*W[0]),longitudes)]
# WE_elev_distances = [geopy.distance.GeodesicDistance(W,[lat,lon]).km for lat,lon in zip(WE["Lat"],WE["Lon"])]
# NS_down_distances = [geopy.distance.GeodesicDistance(N,[lat,lon]).km for lat,lon in zip(latitudes,list(np.ones(len(longitudes))*N[1]))]
# NS_elev_distances = [geopy.distance.GeodesicDistance(N,[lat,lon]).km for lat,lon in zip(NS["Lat"],NS["Lon"])]
# WE_down_distances = np.array(WE_down_distances)
# NS_down_distances = np.array(NS_down_distances)

#%% Plot hypocenters colored by FI

if plot_FI:

    # Determine indices of events with calculated FI and events without calculated FI
    valid_FI_index = np.where(~np.isnan(FI_values))

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("17c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="c)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.2f0.05+lLongitude", "yaf+lDepth", "WSne"],
                    panel=[0, 0])

        pygmt.makecpt(cmap="polar", reverse=False, series=[-1.25, 0.5])
        if size_by_magnitude:
            fig.plot(x=longitudes[valid_FI_index], y=depths[valid_FI_index], size=sizes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes[valid_FI_index], y=depths[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c0.07c", pen="black", transparency=20)
        fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa0.5f0.25+l"Frequency Index"')

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)

        pygmt.makecpt(cmap="polar", reverse=False, series=[-1.25, 0.5])
        if size_by_magnitude:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], size=sizes[valid_FI_index],
                     color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=178.12, y=51.845, style="r2.5/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=178.06, y=51.84, color='gray', style="c0.060c", pen="black")
            fig.plot(x=178.09, y=51.84, color='gray', style="c0.108c", pen="black")
            fig.plot(x=178.12, y=51.84, color='gray', style="c0.156c", pen="black")
            fig.plot(x=178.15, y=51.84, color='gray', style="c0.204c", pen="black")
            fig.plot(x=178.18, y=51.84, color='gray', style="c0.252c", pen="black")
            fig.text(x=178.06, y=51.848, text="M0", justify="CB", font="9p,black")
            fig.text(x=178.09, y=51.848, text="M1", justify="CB", font="9p,black")
            fig.text(x=178.12, y=51.848, text="M2", justify="CB", font="9p,black")
            fig.text(x=178.15, y=51.848, text="M3", justify="CB", font="9p,black")
            fig.text(x=178.18, y=51.848, text="M4", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="b)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xaf+lDepth", "ya0.1f0.025+lLatitude", "wsNE"],
                    panel=[0, 0])

        pygmt.makecpt(cmap="polar", reverse=False, series=[-1.25, 0.5])
        if size_by_magnitude:
            fig.plot(x=depths[valid_FI_index], y=latitudes[valid_FI_index], size=sizes[valid_FI_index],
                     color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=depths[valid_FI_index], y=latitudes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

    fig.show(method="external")
    fig.savefig(save_filepath1 + '_FI.jpg')


#######################
style = 'c1' # c1 or c2

if plot_temporal:

    # Determine indices of events with calculated FI and events without calculated FI
    c1_index = np.where((times > np.datetime64('2021-12-08T00:00:00')) & (times < np.datetime64('2021-12-12T00:00:00')))
    c2_index = np.where((times > np.datetime64('2022-01-21T00:00:00')) & (times < np.datetime64('2022-01-27T00:00:00')))
    c3_index = np.where((times > np.datetime64('2021-12-12T00:00:00')) & (times < np.datetime64('2022-01-21T00:00:00')))

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("17c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="c)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.2f0.05+lLongitude", "yaf+lDepth", "WSne"],
                    panel=[0, 0])

        if style == 'c1':
            fig.plot(x=longitudes[c2_index], y=depths[c2_index], size=sizes[c2_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=longitudes[c3_index], y=depths[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2021-12-08T/2021-12-12T/1d", continuous=True)
            fig.plot(x=longitudes[c1_index], y=depths[c1_index], size=sizes[c1_index], color=times[c1_index], cmap=True, style="c", pen="black", transparency=50)
            fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa1Of5+l"Datetime"')
        elif style == 'c2':
            fig.plot(x=longitudes[c1_index], y=depths[c1_index], size=sizes[c1_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=longitudes[c3_index], y=depths[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2022-01-21T/2022-01-27T/1d", continuous=True)
            fig.plot(x=longitudes[c2_index], y=depths[c2_index], size=sizes[c2_index], color=times[c2_index], cmap=True, style="c", pen="black", transparency=50)
            fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa1Of5+l"Datetime"')

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)

        if style == 'c1':
            fig.plot(x=longitudes[c2_index], y=latitudes[c2_index], size=sizes[c2_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=longitudes[c3_index], y=latitudes[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2021-12-08T/2021-12-12T/1d")
            fig.plot(x=longitudes[c1_index], y=latitudes[c1_index], size=sizes[c1_index], color=times[c1_index], cmap=True, style="c", pen="black", transparency=50)
        elif style == 'c2':
            fig.plot(x=longitudes[c1_index], y=latitudes[c1_index], size=sizes[c1_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=longitudes[c3_index], y=latitudes[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2022-01-21T/2022-01-27T/1d")
            fig.plot(x=longitudes[c2_index], y=latitudes[c2_index], size=sizes[c2_index], color=times[c2_index], cmap=True, style="c", pen="black", transparency=50)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=178.12, y=51.845, style="r2.6/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=178.06, y=51.84, color='gray', style="c0.060c", pen="black")
            fig.plot(x=178.09, y=51.84, color='gray', style="c0.108c", pen="black")
            fig.plot(x=178.12, y=51.84, color='gray', style="c0.156c", pen="black")
            fig.plot(x=178.15, y=51.84, color='gray', style="c0.204c", pen="black")
            fig.plot(x=178.18, y=51.84, color='gray', style="c0.252c", pen="black")
            fig.text(x=178.06, y=51.848, text="M0", justify="CB", font="9p,black")
            fig.text(x=178.09, y=51.848, text="M1", justify="CB", font="9p,black")
            fig.text(x=178.12, y=51.848, text="M2", justify="CB", font="9p,black")
            fig.text(x=178.15, y=51.848, text="M3", justify="CB", font="9p,black")
            fig.text(x=178.18, y=51.848, text="M4", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="b)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xaf+lDepth", "ya0.1f0.025+lLatitude", "wsNE"],
                    panel=[0, 0])

        if style == 'c1':
            fig.plot(x=depths[c2_index], y=latitudes[c2_index], size=sizes[c2_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=depths[c3_index], y=latitudes[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2021-12-08T/2021-12-12T/1d")
            fig.plot(x=depths[c1_index], y=latitudes[c1_index], size=sizes[c1_index], color=times[c1_index], cmap=True, style="c", pen="black", transparency=50)
        elif style == 'c2':
            fig.plot(x=depths[c1_index], y=latitudes[c1_index], size=sizes[c1_index], color='gray', style="c", pen="black", transparency=20)
            fig.plot(x=depths[c3_index], y=latitudes[c3_index], size=sizes[c3_index], color='gray', style="c", pen="black", transparency=20)
            pygmt.makecpt(cmap="viridis", series="2022-01-21T/2022-01-27T/1d")
            fig.plot(x=depths[c2_index], y=latitudes[c2_index], size=sizes[c2_index], color=times[c2_index], cmap=True, style="c", pen="black", transparency=50)

    fig.show(method="external")
    fig.savefig(save_filepath1 + '_T.jpg')

if plot_annual_summary:

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("13c", "13c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 13, 0, 13], projection="X13c/13c", frame='+t"Earthquakes between 2021-12-01 and 2021-12-31"')

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("12c", "3.5c"), autolabel="B)"):

        # Create basemap with correct dimensions
        fig.basemap(region=[178.03, 178.63, -3.5, 10], projection="X12c/-3.5c", frame=["xa0.2f0.05+lLongitude", "yaf+lDepth", "WSne"],
                    panel=[0, 0])

        pygmt.makecpt(cmap="viridis", reverse=True, series=[-4, 10])
        fig.plot(x=longitudes, y=depths, size=sizes, color=depths, cmap=True, style="c", pen="black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Bottom plot: Top down
    with fig.subplot(nrows=1, ncols=1, figsize=("12c", "8c"), autolabel="A)"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=[178.03, 178.63, 51.90, 52.1175], projection="X12c/8c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)
        pygmt.makecpt(cmap="viridis", reverse=True, series=[-4, 10])
        # pygmt.makecpt(cmap="polar", reverse=False, series=[-1.2, 0.3])
        fig.plot(x=longitudes, y=latitudes, size=sizes, color=depths, cmap=True, style="c", pen="black")
        fig.colorbar(position="JTR+o0.5c/-12c+w-12c/0.5c+v", frame='xa2f1+l"Depth (km)"')

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=178.14, y=51.925, style="r3.4/1", color="white", pen="0p", transparency=25)
            fig.plot(x=178.08, y=51.92, color='gray', style="c0.060c", pen="black")
            fig.plot(x=178.11, y=51.92, color='gray', style="c0.108c", pen="black")
            fig.plot(x=178.14, y=51.92, color='gray', style="c0.156c", pen="black")
            fig.plot(x=178.17, y=51.92, color='gray', style="c0.204c", pen="black")
            fig.plot(x=178.20, y=51.92, color='gray', style="c0.252c", pen="black")
            fig.text(x=178.08, y=51.928, text="M0", justify="CB", font="10p,black")
            fig.text(x=178.11, y=51.928, text="M1", justify="CB", font="10p,black")
            fig.text(x=178.14, y=51.928, text="M2", justify="CB", font="10p,black")
            fig.text(x=178.17, y=51.928, text="M3", justify="CB", font="10p,black")
            fig.text(x=178.20, y=51.928, text="M4", justify="CB", font="10p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

    fig.show(method="external")
    fig.savefig(main_dir + 'davidof_annual_summary_depth.pdf')

#
# if plot_annual_summary:
#
#     # Initialize figure
#     fig = pygmt.Figure()
#
#     # Craft title
#     with fig.subplot(nrows=1, ncols=1, figsize=("13c", "13c")):
#
#         with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
#             fig.basemap(region=[0, 13, 0, 13], projection="X13c/13c", frame='+t"Earthquakes between 2021-12-01 and 2021-12-31"')
#
#     # Bottom plot: E-W cross section
#     with fig.subplot(nrows=1, ncols=1, figsize=("12c", "3.5c"), autolabel="B)"):
#
#         # Create basemap with correct dimensions
#         fig.basemap(region=[178.03, 178.63, -3.5, 10], projection="X12c/-3.5c", frame=["xa0.2f0.05+lLongitude", "yaf+lDepth", "WSne"],
#                     panel=[0, 0])
#
#         pygmt.makecpt(cmap="viridis", reverse=True, series=[-4, 10])
#         fig.plot(x=longitudes, y=depths, size=sizes, color=depths, cmap=True, style="c", pen="black")
#
#     # Move plot origin to plot top-left plot
#     fig.shift_origin(yshift="h+0.5c")
#
#     # Bottom plot: Top down
#     with fig.subplot(nrows=1, ncols=1, figsize=("12c", "8.4c"), autolabel="A)"):
#
#         # Create basemap with correct dimensions
#         pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org", MAP_FRAME_TYPE="plain",
#                      FORMAT_GEO_MAP="ddd.x")
#         grid = pygmt.datasets.load_earth_relief(resolution="01s", region=[178.03, 178.63, 51.85, 52.11])
#         fig.basemap(region=[178.03, 178.63, 51.85, 52.11], projection="M178.33/51.97/12c",
#                     frame="WsNe", panel=[0, 0], map_scale="jRB+w10M+c+o0.5c/0.5c")
#         fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)
#         pygmt.makecpt(cmap="viridis", reverse=True, series=[-4, 10])
#         # pygmt.makecpt(cmap="polar", reverse=False, series=[-1.2, 0.3])
#         fig.plot(x=longitudes, y=latitudes, size=sizes, color=depths, cmap=True, style="c", pen="black")
#         fig.colorbar(position="JTR+o0.5c/-12.4c+w-12.5c/0.5c+v", frame='xa2f1+l"Depth (km)"')
#
#         # Plot size indicator if sizing by mag:
#         if size_by_magnitude:
#             fig.plot(x=178.14, y=51.925, style="r3.4/1", color="white", pen="0p", transparency=25)
#             fig.plot(x=178.08, y=51.92, color='gray', style="c0.060c", pen="black")
#             fig.plot(x=178.11, y=51.92, color='gray', style="c0.108c", pen="black")
#             fig.plot(x=178.14, y=51.92, color='gray', style="c0.156c", pen="black")
#             fig.plot(x=178.17, y=51.92, color='gray', style="c0.204c", pen="black")
#             fig.plot(x=178.20, y=51.92, color='gray', style="c0.252c", pen="black")
#             fig.text(x=178.08, y=51.928, text="M0", justify="CB", font="10p,black")
#             fig.text(x=178.11, y=51.928, text="M1", justify="CB", font="10p,black")
#             fig.text(x=178.14, y=51.928, text="M2", justify="CB", font="10p,black")
#             fig.text(x=178.17, y=51.928, text="M3", justify="CB", font="10p,black")
#             fig.text(x=178.20, y=51.928, text="M4", justify="CB", font="10p,black")
#
#         # Plot volcano lat/lon
#         fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")
#
#     fig.show(method="external")
#
# if plot_annual_summary:
#
#     # Initialize figure
#     fig = pygmt.Figure()
#
#     # # Craft title
#     # with fig.subplot(nrows=1, ncols=1, figsize=("13c", "13c")):
#     #
#     #     with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
#     #         fig.basemap(region=[0, 13, 0, 13], projection="X13c/13c", frame='+t"Earthquakes between 2021-12-01 and 2021-12-31"')
#
#     # Bottom plot: Top down
#     with fig.subplot(nrows=1, ncols=1, figsize=("12c", "12c")):
#
#         # Create basemap with correct dimensions
#         pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org", MAP_FRAME_TYPE="plain",
#                      FORMAT_GEO_MAP="ddd.x")
#         fig.basemap(region=REGION_CLEAN, projection="M178.33/51.97/12c",
#                     frame="WsNe", panel=[0, 0], map_scale="jRB+w10M+c+o0.5c/0.5c")
#         grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
#         fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)
#         pygmt.makecpt(cmap="viridis", reverse=True, series=[-4, 10])
#         # pygmt.makecpt(cmap="polar", reverse=False, series=[-1.2, 0.3])
#         fig.plot(x=longitudes, y=latitudes, size=sizes, color=depths, cmap=True, style="c", pen="black")
#         fig.colorbar(position="JTR+o0.5c/-12c+w-12c/0.5c+v", frame='xa2f1+l"Depth (km)"')
#
#         # Plot size indicator if sizing by mag:
#         if size_by_magnitude:
#             fig.plot(x=178.12, y=51.845, style="r2.5/0.9", color="white", pen="0p", transparency=25)
#             fig.plot(x=178.06, y=51.84, color='gray', style="c0.060c", pen="black")
#             fig.plot(x=178.09, y=51.84, color='gray', style="c0.108c", pen="black")
#             fig.plot(x=178.12, y=51.84, color='gray', style="c0.156c", pen="black")
#             fig.plot(x=178.15, y=51.84, color='gray', style="c0.204c", pen="black")
#             fig.plot(x=178.18, y=51.84, color='gray', style="c0.252c", pen="black")
#             fig.text(x=178.06, y=51.848, text="M0", justify="CB", font="9p,black")
#             fig.text(x=178.09, y=51.848, text="M1", justify="CB", font="9p,black")
#             fig.text(x=178.12, y=51.848, text="M2", justify="CB", font="9p,black")
#             fig.text(x=178.15, y=51.848, text="M3", justify="CB", font="9p,black")
#             fig.text(x=178.18, y=51.848, text="M4", justify="CB", font="9p,black")
#
#         # Plot volcano lat/lon
#         fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")
#
#     fig.show(method="external")