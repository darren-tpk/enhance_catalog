#%% PLOT AUGUSTINE

# This script loads in a Augustine catalog of choice and plots the hypocenters in 3 panels
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
elev_profile_dir = main_dir + 'data/dem/'
grid_filepath = elev_profile_dir + 'augustine_hillshade_edifice.tif'
WE_profile_filename = elev_profile_dir + 'we_profile_augustine.csv'
NS_profile_filename = elev_profile_dir + 'ns_profile_augustine.csv'
PEC_events_filepath = main_dir + 'output/augustine5/scan_data/PEC_FI.xml'
relocated_catalog_filepath = main_dir + 'output/augustine5/relocate_catalog/hypodd_reloc.xml'
plot_output_dir = main_dir + 'output/augustine5/relocate_catalog/'
VOLC_LAT = 59.363
VOLC_LON = -153.43
MAX_DEPTH = 99
REGION = [-153.4865, -153.3825, 59.336, 59.3885]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_WE = [-153.4865, -153.3825, -2, 1]
REGION_NS = [-2, 1, 59.336, 59.3885]
MAX_PLOT_DEPTH = 1
##
size_by_magnitude = True
plot_temporal = True
plot_FI = True
plot_relocatable = True
plot_2Dhist = True
save_figs = True

#%% Read all catalogs

# Original catlaog
catalog1 = reader(PEC_events_filepath)
plot_title1 = '+t\"Original Catalog (N=%d)\"' % len(catalog1)
# Relocated catalog
catalog2 = reader(relocated_catalog_filepath)
plot_title2 = '+t\"HypoDD Relocated Catalog (N=%d)\"' % len(catalog2)


#%% Extract information from catalog, filtered by max depth

# Extract hypocenter, time, magnitude and FI information for original catalog
latitudes1 = np.array([event.origins[0].latitude for event in catalog1])
longitudes1 = np.array([event.origins[0].longitude for event in catalog1])
depths1 = np.array([event.origins[0].depth for event in catalog1])  # km
utctimes1 = np.array([event.origins[0].time for event in catalog1])
times1 = np.array([np.datetime64(event.origins[0].time) for event in catalog1])
magnitudes1 = np.array([event.magnitudes[0].mag for event in catalog1])
FI_values1 = []
for event in catalog1:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values1.append(np.nan)
    else:
        FI_values1.append(float(event.comments[-1].text.split('=')[1]))
FI_values1 = np.array(FI_values1)

# Filter all arrays by depth
valid_index1 = np.where(depths1<MAX_DEPTH)
latitudes1 = latitudes1[valid_index1]
longitudes1 = longitudes1[valid_index1]
depths1 = depths1[valid_index1]
times1 = times1[valid_index1]
magnitudes1 = magnitudes1[valid_index1]
FI_values1 = FI_values1[valid_index1]

# Extract hypocenter, time, magnitude and FI information for relocated catalog
latitudes2 = np.array([event.origins[0].latitude for event in catalog2])
longitudes2 = np.array([event.origins[0].longitude for event in catalog2])
depths2 = np.array([event.origins[0].depth for event in catalog2])  # km
utctimes2 = np.array([event.origins[0].time for event in catalog2])
times2 = np.array([np.datetime64(event.origins[0].time) for event in catalog2])
magnitudes2 = np.array([event.magnitudes[0].mag for event in catalog2])
FI_values2 = []
for event in catalog2:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values2.append(np.nan)
    else:
        FI_values2.append(float(event.comments[-1].text.split('=')[1]))
FI_values2 = np.array(FI_values2)

# Filter all arrays by depth
valid_index2 = np.where(depths2<MAX_DEPTH)
latitudes2 = latitudes2[valid_index2]
longitudes2 = longitudes2[valid_index2]
depths2 = depths2[valid_index2]
times2 = times2[valid_index2]
magnitudes2 = magnitudes2[valid_index2]
FI_values2 = FI_values2[valid_index2]

# Normalize the magnitudes based on reasonable size values:
if size_by_magnitude:
    # construct reference sizes
    ref_magnitudes = np.array(range(0,261))
    min_mag_marker = 0
    ref_sizes = []
    for magnitude in list(ref_magnitudes):
        normalized_magnitude = (magnitude - np.min(min_mag_marker))/np.ptp(ref_magnitudes[ref_magnitudes > min_mag_marker])
        ref_size = 0.06 + normalized_magnitude*0.24
        ref_sizes.append(ref_size)
    ref_magnitudes = list(ref_magnitudes)
    # now calculate sizes
    sizes1 = []
    sizes2 = []
    for magnitude in magnitudes1:
        if magnitude < min_mag_marker:
            size = ref_sizes[0]
        else:
            size = ref_sizes[ref_magnitudes.index(int(magnitude*100))]
        sizes1.append(size)
    sizes1 = np.array(sizes1)
    for magnitude in magnitudes2:
        if magnitude < min_mag_marker:
            size = ref_sizes[0]
        else:
            size = ref_sizes[ref_magnitudes.index(int(magnitude*100))]
        sizes2.append(size)
    sizes2 = np.array(sizes2)
    # define size for mag markers
    size0 = ref_sizes[ref_magnitudes.index(0)]
    size1 = ref_sizes[ref_magnitudes.index(100)]
    size2 = ref_sizes[ref_magnitudes.index(200)]

# Prepare cross-section data for SW-NE cross section
W = [VOLC_LAT, -153.83]
E = [VOLC_LAT, -153.03]
WE_x = [W[1], E[1]]
WE_xvec = np.linspace(W[1], E[1], 200)
WE_y = [W[0], E[0]]
WE_yvec = np.linspace(W[0], E[0], 200)

# Prepare cross-section data for NW-SE cross section
N = [59.558, VOLC_LON]
S = [59.163, VOLC_LON]
NS_x = [N[1], S[1]]
NS_xvec = np.linspace(N[1], S[1],600)
NS_y = [N[0], S[0]]
NS_yvec = np.linspace(N[0], S[0],600)
# Load EW and NS elevation profiles (from https://apps.nationalmap.gov/elevation/)
WE = read_csv(WE_profile_filename)
NS = read_csv(NS_profile_filename)


if plot_temporal:

    cpt_series = "2005-11-01T/2006-03-01T/1d"

    # Determine indices of events that fall within colorbar time span
    valid_dt_index1 = np.where((times1 > np.datetime64('2005-11-01')) &
                              (times1 < np.datetime64('2006-03-01')))
    invalid_dt_index1 = np.where((times1 < np.datetime64('2005-11-01')) |
                              (times1 > np.datetime64('2006-03-01')))
    valid_dt_index2 = np.where((times2 > np.datetime64('2005-11-01')) &
                              (times2 < np.datetime64('2006-03-01')))
    invalid_dt_index2 = np.where((times2 < np.datetime64('2005-11-01')) |
                              (times2 > np.datetime64('2006-03-01')))

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c","16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title1)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="C"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1[invalid_dt_index1], y=depths1[invalid_dt_index1], size=sizes1[invalid_dt_index1],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1[invalid_dt_index1], y=depths1[invalid_dt_index1], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=longitudes1[valid_dt_index1], y=depths1[valid_dt_index1], size=sizes1[valid_dt_index1],
                     color=times1[valid_dt_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes1[valid_dt_index1], y=depths1[valid_dt_index1], color=times1[valid_dt_index1], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c","10c"), autolabel="A"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1[invalid_dt_index1], y=latitudes1[invalid_dt_index1], size=sizes1[invalid_dt_index1],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1[invalid_dt_index1], y=latitudes1[invalid_dt_index1], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=longitudes1[valid_dt_index1], y=latitudes1[valid_dt_index1], size=sizes1[valid_dt_index1],
                     color=times1[valid_dt_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes1[valid_dt_index1], y=latitudes1[valid_dt_index1], color=times1[valid_dt_index1],
                     cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.418, style="r2.5/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.415, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.415, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.415, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.415, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.42, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.42, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.42, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.42, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c",yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="B"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        if size_by_magnitude:
            fig.plot(x=depths1[invalid_dt_index1], y=latitudes1[invalid_dt_index1], size=sizes1[invalid_dt_index1],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=depths1[invalid_dt_index1], y=latitudes1[invalid_dt_index1], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=depths1[valid_dt_index1], y=latitudes1[valid_dt_index1], size=sizes1[valid_dt_index1],
                     color=times1[valid_dt_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=depths1[valid_dt_index1], y=latitudes1[valid_dt_index1], color=times1[valid_dt_index1], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    # Shift origin
    fig.shift_origin(yshift="-5.5c")

    # Create island with inset indicator
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "5c"), autolabel="D"):
        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=[-153.582, -153.338, 59.304, 59.426], projection="X5c/5c",
                    frame=["xa0.1f0.05+lLongitude", "ya0.05f0.025+lLatitude", 'wSnE'], panel=[0, 0])
        fig.grdimage(grid=elev_profile_dir + 'augustine_hillshade_island.tif', shading=True, cmap="geo",
                     transparency=50)

        # Plot inset indicator (4 lines)
        fig.plot(x=[-153.4865, -153.3825, -153.3825, -153.4865, -153.4865],
                 y=[59.3885, 59.3885, 59.336, 59.336, 59.3885], pen="1.5p,red")

        # Plot inset scale bar
        fig.plot(x=[-153.35, -153.525], y=[59.312, 59.312], pen="2p,black")
        fig.text(x=(-153.35 -153.525) / 2, y=59.316, text="10 km", justify="CB", font="11p,Helvetica,black")

    # Shift origin again to plot
    fig.shift_origin(xshift="7.7c")

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c","16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title2)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="G"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes2[invalid_dt_index2], y=depths2[invalid_dt_index2], size=sizes2[invalid_dt_index2],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes2[invalid_dt_index2], y=depths2[invalid_dt_index2], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=longitudes2[valid_dt_index2], y=depths2[valid_dt_index2], size=sizes2[valid_dt_index2],
                     color=times2[valid_dt_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes2[valid_dt_index2], y=depths2[valid_dt_index2], color=times2[valid_dt_index2], cmap=True,
                     style="c0.07c", pen="black", transparency=20)
        with pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica,black"):
            fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame="xa1Of+lDatetime")

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c","10c"), autolabel="E"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes2[invalid_dt_index2], y=latitudes2[invalid_dt_index2], size=sizes2[invalid_dt_index2],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes2[invalid_dt_index2], y=latitudes2[invalid_dt_index2], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=longitudes2[valid_dt_index2], y=latitudes2[valid_dt_index2], size=sizes2[valid_dt_index2],
                     color=times2[valid_dt_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes2[valid_dt_index2], y=latitudes2[valid_dt_index2], color=times2[valid_dt_index2],
                     cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.418, style="r2.5/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.415, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.415, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.415, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.415, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.42, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.42, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.42, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.42, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c",yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="F"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005+lLatitude", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        if size_by_magnitude:
            fig.plot(x=depths2[invalid_dt_index2], y=latitudes2[invalid_dt_index2], size=sizes2[invalid_dt_index2],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=depths2[invalid_dt_index2], y=latitudes2[invalid_dt_index2], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="viridis", series=cpt_series)
        if size_by_magnitude:
            fig.plot(x=depths2[valid_dt_index2], y=latitudes2[valid_dt_index2], size=sizes2[valid_dt_index2],
                     color=times2[valid_dt_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=depths2[valid_dt_index2], y=latitudes2[valid_dt_index2], color=times2[valid_dt_index2], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")
    if save_figs:
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS5_augustine_time.png')
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS5_augustine_time.pdf')

if plot_FI:

    # Determine indices of events with calculated FI and events without calculated FI
    valid_FI_index1 = np.where(~np.isnan(FI_values1))
    invalid_FI_index1 = np.where(np.isnan(FI_values1))
    valid_FI_index2 = np.where(~np.isnan(FI_values2))
    invalid_FI_index2 = np.where(np.isnan(FI_values2))

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c","16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title1)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="C"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1[invalid_FI_index1], y=depths1[invalid_FI_index1], size=sizes1[invalid_FI_index1], color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1[invalid_FI_index1], y=depths1[invalid_FI_index1], color="darkgray", cmap=False, style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=longitudes1[valid_FI_index1], y=depths1[valid_FI_index1], size=sizes1[valid_FI_index1], color=FI_values1[valid_FI_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes1[valid_FI_index1], y=depths1[valid_FI_index1], color=FI_values1[valid_FI_index1], cmap=True, style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c","10c"), autolabel="A"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1[invalid_FI_index1], y=latitudes1[invalid_FI_index1], size=sizes1[invalid_FI_index1],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1[invalid_FI_index1], y=latitudes1[invalid_FI_index1], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=longitudes1[valid_FI_index1], y=latitudes1[valid_FI_index1], size=sizes1[valid_FI_index1],
                     color=FI_values1[valid_FI_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes1[valid_FI_index1], y=latitudes1[valid_FI_index1], color=FI_values1[valid_FI_index1], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.418, style="r2.5/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.415, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.415, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.415, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.415, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.42, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.42, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.42, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.42, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c",yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="B"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=depths1[invalid_FI_index1], y=latitudes1[invalid_FI_index1], size=sizes1[invalid_FI_index1],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=depths1[invalid_FI_index1], y=latitudes1[invalid_FI_index1], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=depths1[valid_FI_index1], y=latitudes1[valid_FI_index1], size=sizes1[valid_FI_index1],
                     color=FI_values1[valid_FI_index1], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=depths1[valid_FI_index1], y=latitudes1[valid_FI_index1], color=FI_values1[valid_FI_index1], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    # Shift origin
    fig.shift_origin(yshift="-5.5c")

    # Create island with inset indicator
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "5c"), autolabel="D"):
        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=[-153.582, -153.338, 59.304, 59.426], projection="X5c/5c",
                    frame=["xa0.1f0.05+lLongitude", "ya0.05f0.025+lLatitude", 'wSnE'], panel=[0, 0])
        fig.grdimage(grid=elev_profile_dir + 'augustine_hillshade_island.tif', shading=True, cmap="geo",
                     transparency=50)

        # Plot inset indicator (4 lines)
        fig.plot(x=[-153.4865, -153.3825, -153.3825, -153.4865, -153.4865],
                 y=[59.3885, 59.3885, 59.336, 59.336, 59.3885], pen="1.5p,red")

        # Plot inset scale bar
        fig.plot(x=[-153.35, -153.525], y=[59.312, 59.312], pen="2p,black")
        fig.text(x=(-153.35 -153.525) / 2, y=59.316, text="10 km", justify="CB", font="11p,Helvetica,black")

    # Shift origin again to plot
    fig.shift_origin(xshift="7.7c")

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c","16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title2)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="G"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        # if size_by_magnitude:
        #     fig.plot(x=longitudes2[invalid_FI_index2], y=depths2[invalid_FI_index2], size=sizes2[invalid_FI_index2], color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        # else:
        #     fig.plot(x=longitudes2[invalid_FI_index2], y=depths2[invalid_FI_index2], color="darkgray", cmap=False, style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=longitudes2[valid_FI_index2], y=depths2[valid_FI_index2], size=sizes2[valid_FI_index2], color=FI_values2[valid_FI_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes2[valid_FI_index2], y=depths2[valid_FI_index2], color=FI_values2[valid_FI_index2], cmap=True, style="c0.07c", pen="black", transparency=20)
        fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa0.5f0.25+l"Frequency Index"')

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c","10c"), autolabel="E"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        # if size_by_magnitude:
        #     fig.plot(x=longitudes2[invalid_FI_index2], y=latitudes2[invalid_FI_index2], size=sizes2[invalid_FI_index2],
        #              color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        # else:
        #     fig.plot(x=longitudes2[invalid_FI_index2], y=latitudes2[invalid_FI_index2], color="darkgray", cmap=False,
        #              style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=longitudes2[valid_FI_index2], y=latitudes2[valid_FI_index2], size=sizes2[valid_FI_index2],
                     color=FI_values2[valid_FI_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes2[valid_FI_index2], y=latitudes2[valid_FI_index2], color=FI_values2[valid_FI_index2],
                     cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.418, style="r2.5/0.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.415, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.415, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.415, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.415, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.42, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.42, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.42, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.42, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c",yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="F"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005+lLatitude", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        # if size_by_magnitude:
        #     fig.plot(x=depths2[invalid_FI_index2], y=latitudes2[invalid_FI_index2], size=sizes2[invalid_FI_index2],
        #              color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        # else:
        #     fig.plot(x=depths2[invalid_FI_index2], y=latitudes2[invalid_FI_index2], color="darkgray", cmap=False,
        #              style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1, 0.75])
        if size_by_magnitude:
            fig.plot(x=depths2[valid_FI_index2], y=latitudes2[valid_FI_index2], size=sizes2[valid_FI_index2],
                     color=FI_values2[valid_FI_index2], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=depths2[valid_FI_index2], y=latitudes2[valid_FI_index2], color=FI_values2[valid_FI_index2], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")
    if save_figs:
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS5_augustine_FI.png')
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS5_augustine_FI.pdf')

if plot_relocatable:

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title1)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="C"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1, y=depths1, size=sizes1,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1, y=depths1, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="A"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes1, y=latitudes1, size=sizes1,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes1, y=latitudes1, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.427, style="r2.5/1.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.433, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.433, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.433, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.433, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.438, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.438, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.438, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.438, text="M3", justify="CB", font="9p,black")

        # Plot color indicator
        fig.plot(x=-152.935, y=60.425, style="s0.5", color="orange", pen="0p", transparency=15)
        fig.text(x=-152.925, y=60.425, text="Relocated", justify="LM", font="9p,black")
        fig.plot(x=-152.935, y=60.415, style="s0.5", color="darkgray", pen="0p", transparency=15)
        fig.text(x=-152.925, y=60.415, text="Unrelocated", justify="LM", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="B"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=depths1, y=latitudes1, size=sizes1,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=depths1, y=latitudes1, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    # Shift origin
    fig.shift_origin(yshift="-5.5c")

    # Create island with inset indicator
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "5c"), autolabel="D"):
        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=[-153.582, -153.338, 59.304, 59.426], projection="X5c/5c",
                    frame=["xa0.1f0.05+lLongitude", "ya0.05f0.025+lLatitude", 'wSnE'], panel=[0, 0])
        fig.grdimage(grid=elev_profile_dir + 'augustine_hillshade_island.tif', shading=True, cmap="geo",
                     transparency=50)

        # Plot inset indicator (4 lines)
        fig.plot(x=[-153.4865, -153.3825, -153.3825, -153.4865, -153.4865],
                 y=[59.3885, 59.3885, 59.336, 59.336, 59.3885], pen="1.5p,red")

        # Plot inset scale bar
        fig.plot(x=[-153.35, -153.525], y=[59.312, 59.312], pen="2p,black")
        fig.text(x=(-153.35 -153.525) / 2, y=59.316, text="10 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot right plot
    fig.shift_origin(xshift="7.7c")

    # Define function to remove repeats (this is necessary for loc and reloc files)
    def remove_catalog_repeats(catalogA, catalogB):
        catalogA_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogA]
        catalogB_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogB]
        catalogC = catalogA.copy()
        for i, catalogA_evid in reversed(list(enumerate(catalogA_evids))):
            if catalogA_evid in catalogB_evids:
                catalogC.events.pop(i)
        return catalogC

    # Load in loc catalog and obtain unrelocated catalog by removing repeats
    loc_catalog = reader(main_dir + 'output/augustine5/relocate_catalog/hypodd_loc.xml')
    unrelocated_catalog = remove_catalog_repeats(loc_catalog, catalog2)
    all_loc_catalog_comments = [event.comments for event in loc_catalog]
    # Also add events who don't have any dt.cc or dt.ct links
    relocatable_catalog = reader(main_dir + 'output/augustine5/scan_data/relocatable_catalog_FImag.xml')
    for event in relocatable_catalog:
        if event.comments not in all_loc_catalog_comments and event.origins[
            0].latitude is not None and event.magnitudes != []:
            unrelocated_catalog.events.append(event)
    plot_title3 = '+t\"HypoDD Output (%d/%d relocated)\"' % (len(catalog2),len(unrelocated_catalog)+len(catalog2))

    # Extract hypocenter, time, magnitude and FI information from unrelocated events
    unrelocated_latitudes = np.array([event.origins[0].latitude for event in unrelocated_catalog])
    unrelocated_longitudes = np.array([event.origins[0].longitude for event in unrelocated_catalog])
    unrelocated_depths = np.array([event.origins[0].depth for event in unrelocated_catalog])  # km
    unrelocated_utctimes = np.array([event.origins[0].time for event in unrelocated_catalog])
    unrelocated_times = np.array([np.datetime64(event.origins[0].time) for event in unrelocated_catalog])
    unrelocated_magnitudes = np.array([event.magnitudes[0].mag for event in unrelocated_catalog])
    unrelocated_FI_values = []
    for event in unrelocated_catalog:
        if event.comments[-1].text.split('=')[1] == 'None':
            unrelocated_FI_values.append(np.nan)
        else:
            unrelocated_FI_values.append(float(event.comments[-1].text.split('=')[1]))
    unrelocated_FI_values = np.array(unrelocated_FI_values)

    # Filter all arrays by depth
    valid_index = np.where(unrelocated_depths < MAX_DEPTH)
    unrelocated_latitudes = unrelocated_latitudes[valid_index]
    unrelocated_longitudes = unrelocated_longitudes[valid_index]
    unrelocated_depths = unrelocated_depths[valid_index]
    unrelocated_times = unrelocated_times[valid_index]
    unrelocated_magnitudes = unrelocated_magnitudes[valid_index]
    unrelocated_FI_values = unrelocated_FI_values[valid_index]

    # Normalize the magnitudes based on reasonable size values:
    if size_by_magnitude:
        # calculate sizes
        unrelocated_sizes = []
        for magnitude in unrelocated_magnitudes:
            if magnitude < min_mag_marker:
                size = ref_sizes[0]
            else:
                size = ref_sizes[ref_magnitudes.index(int(magnitude * 100))]
            unrelocated_sizes.append(size)
        unrelocated_sizes = np.array(unrelocated_sizes)

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("17c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title3)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="G"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(WE)):
            delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[WE["Lon"][i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=unrelocated_longitudes, y=unrelocated_depths, size=unrelocated_sizes,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
            fig.plot(x=longitudes2, y=depths2, size=sizes2,
                     color="orange", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=unrelocated_longitudes, y=unrelocated_depths, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)
            fig.plot(x=longitudes2, y=depths2, color="orange", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="E"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=unrelocated_longitudes, y=unrelocated_latitudes, size=unrelocated_sizes,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
            fig.plot(x=longitudes2, y=latitudes2, size=sizes2,
                     color="orange", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=unrelocated_longitudes, y=unrelocated_latitudes, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)
            fig.plot(x=longitudes2, y=latitudes2, color="orange", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-152.9, y=60.427, style="r2.5/1.9", color="white", pen="0p", transparency=25)
            fig.plot(x=-152.935, y=60.433, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-152.912, y=60.433, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-152.889, y=60.433, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.866, y=60.433, color='gray', style="c0.24c", pen="black")
            fig.text(x=-152.935, y=60.438, text="M0", justify="CB", font="9p,black")
            fig.text(x=-152.912, y=60.438, text="M1", justify="CB", font="9p,black")
            fig.text(x=-152.889, y=60.438, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.866, y=60.438, text="M3", justify="CB", font="9p,black")

        # Plot color indicator
        fig.plot(x=-152.935, y=60.425, style="s0.5", color="orange", pen="0p", transparency=15)
        fig.text(x=-152.925, y=60.425, text="Relocated", justify="LM", font="9p,black")
        fig.plot(x=-152.935, y=60.415, style="s0.5", color="darkgray", pen="0p", transparency=15)
        fig.text(x=-152.925, y=60.415, text="Unrelocated", justify="LM", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="F"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005+lLatitude", "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NS)):
            delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[3] - REGION_NS[2]) * 10  # inter-measurement width
            height_km = (MAX_PLOT_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_PLOT_DEPTH+2) * 4.85  # ratio * figure height
            midpoint = MAX_PLOT_DEPTH - 0.5 * height_km
            data = [[midpoint, NS["Lat"][i], height_cm, delta]]
            fig.plot(data=data, style="r", color="gray90", pen="3p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=unrelocated_depths, y=unrelocated_latitudes, size=unrelocated_sizes,
                     color="darkgray", cmap=False, style="c", pen="0.7p,black", transparency=20)
            fig.plot(x=depths2, y=latitudes2, size=sizes2,
                     color="orange", cmap=False, style="c", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=unrelocated_depths, y=unrelocated_latitudes, color="darkgray", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)
            fig.plot(x=depths2, y=latitudes2, color="orange", cmap=False,
                     style="c0.12c", pen="0.7p,black", transparency=20)

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")

    if save_figs:
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/fig5_augustine_relocatable.png')
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/fig5_augustine_relocatable.pdf')

if plot_2Dhist:

    def pygmt_hist2d(X, Y, FI_values, REGION, orientation, include_zeros=False):
        import matplotlib.pyplot as plt
        if orientation == 'h':
            out = plt.hist2d(X, Y, bins=[100, 50], range=[[REGION[0], REGION[1]], [REGION[2], REGION[3]]])
        elif orientation == 'v':
            out = plt.hist2d(X, Y, bins=[50, 100], range=[[REGION[0], REGION[1]], [REGION[2], REGION[3]]])
        elif orientation == 't':
            out = plt.hist2d(X, Y, bins=[100, 100], range=[[REGION[0], REGION[1]], [REGION[2], REGION[3]]])
        evgrid = out[0]
        xedges = out[1]
        yedges = out[2]
        delta_x = xedges[1] - xedges[0]
        delta_y = yedges[1] - yedges[0]
        figrid = np.zeros(np.shape(evgrid))
        xgrid = np.zeros(np.shape(evgrid))
        ygrid = np.zeros(np.shape(evgrid))
        for i, xedge in enumerate(xedges[:-1]):
            for j, yedge in enumerate(yedges[:-1]):
                logic = (X > xedge) & (X < (xedge + delta_x)) & (Y > yedge) & (Y < (yedge + delta_y))
                mean_fi = np.nanmedian(FI_values[logic])
                xgrid[i, j] = xedge + delta_x / 2
                ygrid[i, j] = yedge + delta_y / 2
                figrid[i, j] = mean_fi
        if not include_zeros:
            evvec = evgrid.flatten()[~np.isnan(figrid.flatten())]
            xvec = xgrid.flatten()[~np.isnan(figrid.flatten())]
            yvec = ygrid.flatten()[~np.isnan(figrid.flatten())]
            fivec = figrid.flatten()[~np.isnan(figrid.flatten())]
            return evvec, xvec, yvec, fivec
        else:
            evvec = evgrid.flatten()
            xvec = xgrid.flatten()
            yvec = ygrid.flatten()
            return evvec, xvec, yvec

    # # For PEC catalog, dont need to use unrelocated info
    hist_longitudes1 = longitudes1
    hist_latitudes1 = latitudes1
    hist_depths1 = depths1
    hist_FI_values1 = FI_values1
    # # Merge unrelocated and relocated info for HypoDD output
    hist_longitudes2 = np.concatenate((longitudes2,unrelocated_longitudes))
    hist_latitudes2 = np.concatenate((latitudes2,unrelocated_latitudes))
    hist_depths2 = np.concatenate((depths2,unrelocated_depths))
    hist_FI_values2 = np.concatenate((FI_values2,unrelocated_FI_values))

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title1)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="C"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_longitudes1, hist_depths1, hist_FI_values1, REGION_WE, 'h', include_zeros=True)
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="A"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_longitudes1, hist_latitudes1, hist_FI_values1, REGION, 't', include_zeros=True)
        c = pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="B"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005", "wsNE"],
                    panel=[0, 0])

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_depths1, hist_latitudes1, hist_FI_values1, REGION_NS, 'v', include_zeros=True)
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    # Move plot origin to plot inset
    fig.shift_origin(yshift="-5.5c")

    # Create island with inset indicator
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "5c"), autolabel="D"):
        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=[-153.582, -153.338, 59.304, 59.426], projection="X5c/5c",
                    frame=["xa0.1f0.05+lLongitude", "ya0.05f0.025+lLatitude", 'wSnE'], panel=[0, 0])
        fig.grdimage(grid=elev_profile_dir + 'augustine_hillshade_island.tif', shading=True, cmap="geo",
                     transparency=50)

        # Plot inset indicator (4 lines)
        fig.plot(x=[-153.4865, -153.3825, -153.3825, -153.4865, -153.4865],
                 y=[59.3885, 59.3885, 59.336, 59.336, 59.3885], pen="1.5p,red")

        # Plot inset scale bar
        fig.plot(x=[-153.35, -153.525], y=[59.312, 59.312], pen="2p,black")
        fig.text(x=(-153.35 -153.525) / 2, y=59.316, text="10 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot right plot
    fig.shift_origin(xshift="7.7c")

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("16.5c", "16.5c")):

        with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
            fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title3)

    # Bottom plot: E-W cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="G"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.02f0.01+lLongitude", "ya1f0.25+lDepth", "WSne"],
                    panel=[0, 0])

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_longitudes2, hist_depths2, hist_FI_values2, REGION_WE, 'h', include_zeros=True)
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')
        fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa10f5+l"Number of Events"')

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="E"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["xa0.02f0.01+lLongitude", "ya0.01f0.005+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=50)

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_longitudes2, hist_latitudes2, hist_FI_values2, REGION, 't', include_zeros=True)
        c = pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.5c", pen="black")

        # Plot scale bar
        fig.plot(x=[-153.3875,-153.4225], y=[59.34,59.34], pen="2p,black")
        fig.text(x=(-153.3875-153.4225)/2, y=59.3412, text="2 km", justify="CB", font="11p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="F"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xa1f0.25+lDepth", "ya0.01f0.005+lLatitude", "wsNE"],
                    panel=[0, 0])

        # Plot event count histogram
        evvec, xvec, yvec = pygmt_hist2d(hist_depths2, hist_latitudes2, hist_FI_values2, REGION_NS, 'v', include_zeros=True)
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 50])
        fig.plot(x=xvec[evvec>0], y=yvec[evvec>0], color=evvec[evvec>0], cmap=True, style="s0.14c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec==0], y=yvec[evvec==0], style="s0.14c", pen='0.1p,darkgrey')

        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")
    if save_figs:
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS9_augustine_2dhist.png')
        fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/figS9_augustine_2dhist.pdf')

