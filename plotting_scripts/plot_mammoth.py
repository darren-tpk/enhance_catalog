#%% PLOT MAMMOTH

# This script loads in a MM catalog of choice and plots the hypocenters in 3 panels
# (a NW-SE cross section, a NE-SW cross section, and a top-down view of the hypocenters)
# The user will have to specify and tweak lat and lon limits, and region settings
# For detailed description of plotting syntax, refer to PyGMT documentation

# Import all dependencies
from phase_processing.read_hypoddpha import read_hypoddpha
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
SWNE_profile_filename = elev_profile_dir + 'sw_ne_profile_MM.csv'
NWSE_profile_filename = elev_profile_dir + 'nw_se_profile_MM.csv'
PEC_events_dir = main_dir + 'data/ncedc/'
cores_filename = main_dir + 'output/mammoth3/convert_redpy/core_catalog_picked.xml'
unmatched_PEC_events_filepath = main_dir + 'output/mammoth3/convert_redpy/unmatched_PEC_events_redpy.xml'
relocated_catalog_filepath = main_dir + 'output/mammoth3/relocate_catalog/hypodd_reloc.xml'
plot_output_dir = main_dir + 'output/mammoth2/relocate_catalog/'
VOLC_LAT = 37.631
VOLC_LON = -119.032
MAX_DEPTH = 30  # km
REGION = [-119.20, -119.00, 37.55, 37.70]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_CLEAN[0] = REGION_CLEAN[0] - 0.025
REGION_CLEAN[1] = REGION_CLEAN[1] + 0.025
REGION_CLEAN[2] = REGION_CLEAN[2] - 0.0125
REGION_CLEAN[3] = REGION_CLEAN[3] + 0.0125
REGION_SWNE = [0, 20, -5, MAX_DEPTH]
REGION_NWSE = [0, 20, -5, MAX_DEPTH]
START_TIME = UTCDateTime(2012,10,1,0,0,0)  # reference time for "days" colorbar
size_by_magnitude = False
plot_temporal = True
plot_FI = False
migration_track = True
swarm_focus = False

#%% Define functions

# Calculate down-cross-section distances
def calc_down_distances(lats,lons,section_yvec,section_xvec):
    import geopy.distance
    down_distances = []
    for lat, lon in zip (lats,lons):
        all_distances = [geopy.distance.GeodesicDistance([lat,lon],[section_yvec[i],section_xvec[i]]).km for i in range(len(section_xvec))]
        shortest_distance_index = np.argmin(all_distances)
        down_distance = geopy.distance.GeodesicDistance([section_yvec[0],section_xvec[0]],[section_yvec[shortest_distance_index],section_xvec[shortest_distance_index]]).km
        down_distances.append(down_distance)
    return down_distances

#%% Read all catalogs

# (1) Pre-existing catalog
PEC_hypoi = PEC_events_dir + 'mammoth_20121001_20130131_hypoi.txt'
PEC_hypoddpha = PEC_events_dir + 'mammoth_20121001_20130131_hypoddpha.txt'
PEC_events = read_hypoddpha(PEC_hypoi, PEC_hypoddpha, channel_convention=True)

# (2) All located templates
cores = reader(cores_filename)
unmatched_PEC_events = reader(unmatched_PEC_events_filepath)
templates = cores + unmatched_PEC_events
located_templates = Catalog()
for template in templates:
    if template.origins[0].latitude is not None:
        located_templates += template

# (3) GrowClust relocated catalog
growclust_relocated_catalog = reader(relocated_catalog_filepath)  # min_cc = 0.75
#%% Extract information from catalog, filtered by max depth

# (4) HypoDD relocated catalog 1
hypodd_loc_filtered = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth2/relocate_catalog/hypodd_alicia_settings_adapted/hypodd_loc_filtered.xml')
hypodd_reloc = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth2/relocate_catalog/hypodd_alicia_settings_adapted/hypodd_reloc.xml')
hypodd_relocated_catalog_1 = hypodd_loc_filtered + hypodd_reloc
hypodd_relocated_catalog = hypodd_reloc

# (5) HypoDD relocated catalog 2
hypodd_loc_filtered = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth2/relocate_catalog/hypodd_jeremy_settings_adapted/hypodd_loc_filtered.xml')
hypodd_reloc = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth2/relocate_catalog/hypodd_jeremy_settings_adapted/hypodd_reloc.xml')
hypodd_relocated_catalog_2 = hypodd_loc_filtered + hypodd_reloc

# (6) HypoDD relocated catalog 3
hypodd_loc_filtered = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth3/relocate_catalog/hypodd_loc_filtered.xml')
hypodd_reloc = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth3/relocate_catalog/hypodd_reloc.xml')
hypodd_relocated_catalog_3 = hypodd_loc_filtered + hypodd_reloc

# (7) Alicia's catalog
hypodd_reloc = reader('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/mammoth2/alicia_reloc.xml')
hypodd_relocated_catalog_AHE = hypodd_reloc

# Choose catalog
catalog = hypodd_relocated_catalog_3
if catalog == PEC_events:
    plot_title = '+t\"Original Catalog (N=%d)\"' % len(catalog)
    save_filepath = plot_output_dir + 'original_catalog_N' + str(len(PEC_events)) + '.jpg'
elif catalog == located_templates:
    plot_title = '+t\"Located Templates (N=%d)\"' % len(catalog)
    save_filepath = plot_output_dir + 'located_templates_N' + str(len(located_templates)) + '.jpg'
elif catalog == growclust_relocated_catalog:
    plot_title = '+t\"GrowClust Relocated Catalog (N=%d)\"' % len(catalog)
    # save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'
elif catalog == hypodd_relocated_catalog_1:
    plot_title = '+t\"HypoDD Relocated Catalog 1 (N=%d)\"' % len(catalog)
    # save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'
elif catalog == hypodd_relocated_catalog_2:
    plot_title = '+t\"HypoDD Relocated Catalog 2 (N=%d)\"' % len(catalog)
    # save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'
elif catalog == hypodd_relocated_catalog_3:
    plot_title = '+t\"HypoDD Relocated Catalog (N=%d)\"' % len(catalog)
    # save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'
elif catalog == hypodd_relocated_catalog_AHE:
    plot_title = '+t\"AHE Relocated Catalog (N=%d)\"' % len(catalog)
    # save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'

# Extract hypocenter, time, magnitude and FI information
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog])  # km
utctimes = np.array([event.origins[0].time for event in catalog])
times = np.array([np.datetime64(event.origins[0].time) for event in catalog])
days = ((utctimes - START_TIME) / 86400).astype(int)
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
days = days[valid_index]
magnitudes = magnitudes[valid_index]
if plot_FI:
    FI_values = FI_values[valid_index]

# Normalize the magnitudes based on reasonable size values:
if size_by_magnitude:
    # # Sanity check via histogram
    # _ = plt.hist(magnitudes, bins='auto')  # arguments are passed to np.histogram
    # plt.title("Histogram with 'auto' bins")
    # plt.show()
    sizes = []
    min_mag_marker = 0
    for magnitude in magnitudes:
        if magnitude < min_mag_marker:
            size = 0.06
        else:
            normalized_magnitude = (magnitude - np.min(min_mag_marker))/np.ptp(magnitudes[magnitudes > min_mag_marker])
            size = 0.06 + normalized_magnitude*0.14
        sizes.append(size)
    sizes = np.array(sizes)
    # This is only valid from the current spread of magnitudes

# Prepare cross-section data for SW-NE cross section
SW = [37.55, -119.15]
NE = [37.6675, -118.98]
SWNE_x = [SW[1], NE[1]]
SWNE_xvec = np.linspace(SW[1], NE[1],600)
SWNE_y = [SW[0], NE[0]]
SWNE_yvec = np.linspace(SW[0], NE[0],600)

# Prepare cross-section data for NW-SE cross section
NW = [37.71, -119.15]
SE = [37.595, -118.98]
NWSE_x = [NW[1], SE[1]]
NWSE_xvec = np.linspace(NW[1], SE[1],600)
NWSE_y = [NW[0], SE[0]]
NWSE_yvec = np.linspace(NW[0], SE[0],600)

# Load EW and NS elevation profiles (from https://apps.nationalmap.gov/elevation/)
SWNE = read_csv(SWNE_profile_filename)
NWSE = read_csv(NWSE_profile_filename)

# Calculate down-section distance for all catalog events and elevation profile
SWNE_down_distances = calc_down_distances(latitudes,longitudes,SWNE_yvec,SWNE_xvec)
SWNE_elev_distances = calc_down_distances(SWNE["Lat"],SWNE["Lon"],SWNE_yvec,SWNE_xvec)
NWSE_down_distances = calc_down_distances(latitudes,longitudes,NWSE_yvec,NWSE_xvec)
NWSE_elev_distances = calc_down_distances(NWSE["Lat"],NWSE["Lon"],NWSE_yvec,NWSE_xvec)
SWNE_down_distances = np.array(SWNE_down_distances)
NWSE_down_distances = np.array(NWSE_down_distances)

#%% Plot hypocenters using PyGMT
if plot_temporal:

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("17.5c","9c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17.5, 0, 9.5], projection="X17.5c/9.5c", frame=plot_title)

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("8.5c","8.5c")):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=REGION_CLEAN, projection="X8.5c/8.5c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)

        # Plot earthquakes
        pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
        if size_by_magnitude:
            fig.plot(x=longitudes, y=latitudes, color=times, cmap=True, size=sizes, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes, y=latitudes, color=times, cmap=True, style="c0.08c", pen="black", transparency=20)
        with pygmt.config(FONT_ANNOT_PRIMARY="7p,Helvetica,black"):
            fig.colorbar(position="JBL+o-8.5c/0.5c+w8.5c/0.4c+h",frame="xa1Of")

        # # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-119.195, y=37.553, style="r1.8/0.8", color="white", pen="0p", transparency=25)
            fig.plot(x=-119.21, y=37.55, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-119.195, y=37.55, color='gray', style="c0.13c", pen="black")
            fig.plot(x=-119.18, y=37.55, color='gray', style="c0.20c", pen="black")
            fig.text(x=-119.21, y=37.555, text="M0", justify="CB", font="9p,black")
            fig.text(x=-119.195, y=37.555, text="M1", justify="CB", font="9p,black")
            fig.text(x=-119.18, y=37.555, text="M2", justify="CB", font="9p,black")

        # Plot cross-section indicators
        fig.plot(x=SWNE_x, y=SWNE_y, pen="0.9p,black,-", transparency=40)
        fig.text(x=SW[1], y=SW[0], text="A", justify="RT", font="11p,black")
        fig.text(x=NE[1], y=NE[0], text="B", justify="CB", font="11p,black")
        fig.plot(x=NWSE_x, y=NWSE_y, pen="0.9p,black,-", transparency=40)
        fig.text(x=NW[1], y=NW[0], text="C", justify="RT", font="11p,black")
        fig.text(x=SE[1], y=SE[0], text="D", justify="CT", font="11p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

        # Plot scale bar
        fig.plot(x=[-118.985,-119.0975], y=[37.547,37.547], pen="2p,black")
        fig.text(x=(-118.985-119.0975)/2, y=37.552, text="10 km", justify="CB", font="8.5p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="9c",yshift="-1.5c")

    with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_SWNE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', "ya5f1", "wsNe"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(SWNE)):
            delta = (SWNE_elev_distances[i] - SWNE_elev_distances[i - 1]) / (REGION_SWNE[1]-REGION_SWNE[0]) * 4  # inter-measurement width
            height_km = (MAX_DEPTH + 0.001 * SWNE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            midpoint = MAX_DEPTH - 0.5 * height_km
            data = [[SWNE_elev_distances[i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
        if size_by_magnitude:
            fig.plot(x=SWNE_down_distances, y=depths, size=sizes, style="c", color=times, cmap=True, pen="black", transparency=20)
        else:
            fig.plot(x=SWNE_down_distances, y=depths, style="c0.08c", color=times, cmap=True, pen="black", transparency=20)

        # Plot shallow-deep separator
        fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
        fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
        fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

        # Plot elevation profile
        fig.plot(x=SWNE_elev_distances, y=-0.001*SWNE["Elev(m)"], pen="1.5p,black")
        fig.text(x=0.5, y=-4.5, text="A", justify="LT", font="13p,black")
        fig.text(x=19.5, y=-4.5, text="B", justify="RT", font="13p,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="4.5c")

    with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NWSE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', 'ya5f1+l"Depth (km)"', "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NWSE)):
            delta = (NWSE_elev_distances[i] - NWSE_elev_distances[i - 1]) / (REGION_NWSE[1]-REGION_NWSE[0]) * 4  # inter-measurement width
            height_km = (MAX_DEPTH + 0.001 * NWSE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            midpoint = MAX_DEPTH - 0.5 * height_km
            data = [[NWSE_elev_distances[i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
        if size_by_magnitude:
            fig.plot(x=NWSE_down_distances, y=depths, size=sizes, style="c", color=times, cmap=True, pen="black", transparency=20)
        else:
            fig.plot(x=NWSE_down_distances, y=depths, style="c0.08c", color=times, cmap=True, pen="black", transparency=20)

        # Plot shallow-deep separator
        fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
        fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
        fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

        # Plot elevation profile
        fig.plot(x=NWSE_elev_distances, y=-0.001 * NWSE["Elev(m)"], pen="1.5p,black")
        fig.text(x=0.5, y=-4.5, text="C", justify="LT", font="13p,black")
        fig.text(x=19.5, y=-4.5, text="D", justify="RT", font="13p,black")

    fig.show(method="external")
    fig.savefig('/Users/darrentpk/Desktop/figures/paper/mammoth_PEC_nomag.pdf')
    fig.savefig('/Users/darrentpk/Desktop/figures/paper/mammoth_PEC_nomag.png')

#%% Plot hypocenters colored by FI

if plot_FI:

    # Determine indices of events with calculated FI and events without calculated FI
    valid_FI_index = np.where(~np.isnan(FI_values))
    invalid_FI_index = np.where(np.isnan(FI_values))

    # # Bin the events with FI into a LF-Hybrid-HF color scheme
    # FI_class = []
    # for FI_value in FI_values[valid_FI_index]:
    #     if FI_value < -0.45:
    #         FI_class.append('LF')
    #     elif FI_value > -0.10:
    #         FI_class.append('HF')
    #     else:
    #         FI_class.append('Hybrid')
    # FI_class = pd.Categorical(FI_class)

    # Initialize figure
    fig = pygmt.Figure()

    # Craft title
    with fig.subplot(nrows=1, ncols=1, figsize=("17.5c","9c")):

        with pygmt.config(MAP_FRAME_PEN='white',MAP_TICK_PEN_PRIMARY='white',FONT_ANNOT='white'):
            fig.basemap(region=[0, 17.5, 0, 9.5], projection="X17.5c/9.5c", frame=plot_title)

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("8.5c","8.5c")):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=REGION_CLEAN, projection="X8.5c/8.5c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=40)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], size=sizes[invalid_FI_index], color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], color="darkgray", cmap=False, style="x0.12c", pen="0.7p,black", transparency=20)
        #pygmt.makecpt(cmap="seis", color_model="+cLF,Hybrid,HF", series=(0,2,1))
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1.00, 1.00])
        if size_by_magnitude:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], size=sizes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c0.07c", pen="black", transparency=20)
        fig.colorbar(position="JBL+o-8.5c/0.5c+w8.5c/0.4c+h", frame='xa0.5f0.25+l"Frequency Index"')

        # # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-119.195, y=37.553, style="r1.8/0.8", color="white", pen="0p", transparency=25)
            fig.plot(x=-119.21, y=37.55, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-119.195, y=37.55, color='gray', style="c0.13c", pen="black")
            fig.plot(x=-119.18, y=37.55, color='gray', style="c0.20c", pen="black")
            fig.text(x=-119.21, y=37.555, text="M0", justify="CB", font="9p,black")
            fig.text(x=-119.195, y=37.555, text="M1", justify="CB", font="9p,black")
            fig.text(x=-119.18, y=37.555, text="M2", justify="CB", font="9p,black")

        # Plot cross-section indicators
        fig.plot(x=SWNE_x, y=SWNE_y, pen="0.9p,black,-", transparency=40)
        fig.text(x=SW[1], y=SW[0], text="A", justify="RT", font="11p,black")
        fig.text(x=NE[1], y=NE[0], text="B", justify="CB", font="11p,black")
        fig.plot(x=NWSE_x, y=NWSE_y, pen="0.9p,black,-", transparency=40)
        fig.text(x=NW[1], y=NW[0], text="C", justify="RT", font="11p,black")
        fig.text(x=SE[1], y=SE[0], text="D", justify="CT", font="11p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

        # Plot scale bar
        fig.plot(x=[-118.985,-119.0975], y=[37.547,37.547], pen="2p,black")
        fig.text(x=(-118.985-119.0975)/2, y=37.552, text="10 km", justify="CB", font="8.5p,Helvetica,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="9c",yshift="-1.5c")

    with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_SWNE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', "ya5f1", "wsNe"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(SWNE)):
            delta = (SWNE_elev_distances[i] - SWNE_elev_distances[i - 1]) / (REGION_SWNE[1]-REGION_SWNE[0]) * 4  # inter-measurement width
            height_km = (MAX_DEPTH + 0.001 * SWNE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            midpoint = MAX_DEPTH - 0.5 * height_km
            data = [[SWNE_elev_distances[i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=SWNE_down_distances[invalid_FI_index], y=depths[invalid_FI_index], size=sizes[invalid_FI_index], color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=SWNE_down_distances[invalid_FI_index], y=depths[invalid_FI_index], color="darkgray", cmap=False, style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1.00,1.00])
        if size_by_magnitude:
            fig.plot(x=SWNE_down_distances[valid_FI_index], y=depths[valid_FI_index], size=sizes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=SWNE_down_distances[valid_FI_index], y=depths[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c0.07c", pen="black", transparency=20)

        # Plot shallow-deep separator
        fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
        fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
        fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

        # Plot elevation profile
        fig.plot(x=SWNE_elev_distances, y=-0.001*SWNE["Elev(m)"], pen="1.5p,black")
        fig.text(x=0.5, y=-4.5, text="A", justify="LT", font="13p,black")
        fig.text(x=19.5, y=-4.5, text="B", justify="RT", font="13p,black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="4.5c")

    with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NWSE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', 'ya5f1+l"Depth (km)"', "wsNE"],
                    panel=[0, 0])

        # Plot landfill
        for i in range(1, len(NWSE)):
            delta = (NWSE_elev_distances[i] - NWSE_elev_distances[i - 1]) / (REGION_NWSE[1]-REGION_NWSE[0]) * 4  # inter-measurement width
            height_km = (MAX_DEPTH + 0.001 * NWSE["Elev(m)"][i])  # depth - negative elevation
            height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            midpoint = MAX_DEPTH - 0.5 * height_km
            data = [[NWSE_elev_distances[i], midpoint, delta, height_cm]]
            fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=NWSE_down_distances[invalid_FI_index], y=depths[invalid_FI_index], size=sizes[invalid_FI_index], color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=NWSE_down_distances[invalid_FI_index], y=depths[invalid_FI_index], color="darkgray", cmap=False, style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=True, series=[-1.00,1.00])
        if size_by_magnitude:
            fig.plot(x=NWSE_down_distances[valid_FI_index], y=depths[valid_FI_index], size=sizes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=NWSE_down_distances[valid_FI_index], y=depths[valid_FI_index], color=FI_values[valid_FI_index], cmap=True, style="c0.07c", pen="black", transparency=20)

        # Plot shallow-deep separator
        fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
        fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
        fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

        # Plot elevation profile
        fig.plot(x=NWSE_elev_distances, y=-0.001 * NWSE["Elev(m)"], pen="1.5p,black")
        fig.text(x=0.5, y=-4.5, text="C", justify="LT", font="13p,black")
        fig.text(x=19.5, y=-4.5, text="D", justify="RT", font="13p,black")

    fig.show(method="external")
    fig.savefig('/Users/darrentpk/Desktop/figures/paper/figS2x_mammoth_relocations_FI.pdf')

#%% Plot hypocenters using distance along A-B as coloration
if migration_track:

    # Define time markings and durations on Alicia's plot
    key_times = [UTCDateTime(2012, 10, 29), UTCDateTime(2012, 10, 31),
                 UTCDateTime(2012, 11, 3), UTCDateTime(2012, 11, 18),
                 UTCDateTime(2012, 12, 7), UTCDateTime(2012, 12, 15),
                 UTCDateTime(2012, 12, 26), UTCDateTime(2013, 1, 4)]
    key_durations = [2.0, 1.5, 1.5, 1.0, 1.5, 2.0, 2.0, 1.0]  # days from respective key_time in key_times
    key_cumulative_detections = key_durations.copy()
    key_cumulative_detections.insert(0,0.0)
    key_cumulative_durations = np.cumsum(key_cumulative_detections)  # also calculate cumulative for plotting

    # Construct lists storing plotting points in cumulative days
    migration_x = []
    migration_y = []
    migration_z = []

    for i, key_time in enumerate(key_times):
        key_end = key_time + key_durations[i]*86400
        key_index = np.where((utctimes > key_time) & (utctimes < key_end))

        migration_x.extend(key_cumulative_durations[i] + ((utctimes[key_index] - key_time) / 86400))
        migration_y.extend(depths[key_index])
        migration_z.extend(SWNE_down_distances[key_index])

    # Create scatter plot
    fig, ax1 = plt.subplots(figsize=(15,4))
    main = ax1.scatter(migration_x,migration_y,s=8,c=migration_z,cmap='viridis')
    main.set_clim([9, 14.5])
    cb = fig.colorbar(main)
    cb.set_label('Distance along A-B (km)',fontsize=13)
    for key_cumulative_duration in key_cumulative_durations:
        ax1.axvline(key_cumulative_duration,color='grey',alpha=0.5)
    ax1.set_xlim([0,key_cumulative_durations[-1]])
    ax1.set_xticks(np.arange(0,key_cumulative_durations[-1].round()+1,1))
    ax2 = ax1.twiny()
    ax2.set_xticks(key_cumulative_durations)
    ax2_labels = key_times.copy()
    ax2_labels.append(key_times[-1]+key_durations[-1]*86400)
    ax2_labels = [ax2_label.strftime('%Y-%m-%d') for ax2_label in ax2_labels]
    ax2.set_xticklabels(ax2_labels, rotation=30, ha='left')
    ax1.set_ylim([10,30])
    ax1.set_yticks(np.arange(10,30+5,5))
    ax1.set_ylabel('Depth (km)',fontsize=15)
    ax1.set_xlabel('Concatenated swarm time (days)',fontsize=15)
    ax1.set_title('Hypocenter migration in depth with time',fontsize=15)
    ax1.invert_yaxis()
    fig.savefig(plot_output_dir + 'migration_track_cc75.jpg',bbox_inches='tight')

#%% Plot hypocenters using PyGMT
if swarm_focus:

    # Define start and end of hypocenter times, time step, swarm duration to be plotted
    ec_start = UTCDateTime(2012,10,1)
    ec_end = UTCDateTime(2013,2,1)
    time_step = 60*60 # 1 hour
    swarm_hours = 40  # hours
    ec_time_blocks = int((ec_end - ec_start) / time_step)

    # Loop over time blocks
    for i in range(ec_time_blocks+1):

        # Define swarm reference
        swarm_reference = ec_start + (i*time_step)

        # Define some variables
        swarm_destination = "/Users/darrentpk/Desktop/frames/"
        swarm_framename = swarm_destination + swarm_reference.strftime("%Y_%m_%d_%H_%M_%S") +'.png'
        print('Now at ' + swarm_reference.strftime("%Y_%m_%d_%H_%M_%S") + '...')

        # Determine swarm start and swarm end
        swarm_start = swarm_reference - (swarm_hours * 3600)
        swarm_end = swarm_reference

        # Find out which hypocenters to color
        swarm_index = np.where((utctimes>swarm_start) & (utctimes<swarm_end))
        regular_index = np.where((utctimes<swarm_start))

        # Now compute hours from swarm end
        hours = -(utctimes[swarm_index] - swarm_reference) / 3600

        # Initialize figure
        fig = pygmt.Figure()

        # Craft title
        with fig.subplot(nrows=1, ncols=1, figsize=("17.5c", "9c")):

            with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
                fig.basemap(region=[0, 17.5, 0, 9.5], projection="X17.5c/9.5c", frame=plot_title)

        # Left plot: Top-down view
        with fig.subplot(nrows=1, ncols=1, figsize=("8.5c","8.5c")):

            # Create basemap with correct dimensions
            grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
            fig.basemap(region=REGION_CLEAN, projection="X8.5c/8.5c",
                        frame=["x+lLongitude", "y+lLatitude", "WsNe"], panel=[0, 0])
            fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=25)

            # Plot earthquakes
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=longitudes[regular_index], y=latitudes[regular_index], color="lightgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            pygmt.makecpt(cmap="inferno", reverse=False, series=[0,swarm_hours])
            if len(longitudes[swarm_index]) > 0:
                fig.plot(x=longitudes[swarm_index], y=latitudes[swarm_index], color=hours.tolist(), cmap=True, style="c0.07c", pen="black", transparency=20)
            fig.colorbar(position="JBL+o-8.5c/0.5c+w-8.5c/0.4c+h",frame='xa10f5+l"Age (hours)"')
            fig.text(x=-118.98,y=37.54, text=swarm_reference.strftime("%Y/%m/%d %H:%M:%S"), justify="RB", font="15p,yellow")

            # Plot cross-section indicators
            fig.plot(x=SWNE_x, y=SWNE_y, pen="0.9p,black,-", transparency=40)
            fig.text(x=SW[1], y=SW[0], text="A", justify="RT", font="11p,black")
            fig.text(x=NE[1], y=NE[0], text="B", justify="CB", font="11p,black")
            fig.plot(x=NWSE_x, y=NWSE_y, pen="0.9p,black,-", transparency=40)
            fig.text(x=NW[1], y=NW[0], text="C", justify="RT", font="11p,black")
            fig.text(x=SE[1], y=SE[0], text="D", justify="CT", font="11p,black")

            # Plot volcano lat/lon
            fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

        # Move plot origin to plot cross-section plot
        fig.shift_origin(xshift="9c",yshift="-1.5c")

        with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

            # Create basemap with correct dimensions
            fig.basemap(region=REGION_SWNE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', "ya5f1", "wsNe"],
                        panel=[0, 0])

            # # Plot landfill
            # for i in range(1, len(SWNE)):
            #     delta = (SWNE_elev_distances[i] - SWNE_elev_distances[i - 1]) / (REGION_SWNE[1]-REGION_SWNE[0]) * 4  # inter-measurement width
            #     height_km = (MAX_DEPTH + 0.001 * SWNE["Elev(m)"][i])  # depth - negative elevation
            #     height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            #     midpoint = MAX_DEPTH - 0.5 * height_km
            #     data = [[SWNE_elev_distances[i], midpoint, delta, height_cm]]
            #     fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

            # Plot earthquakes
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=SWNE_down_distances[regular_index], y=depths[regular_index], color="lightgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            if len(longitudes[swarm_index]) > 0:
                pygmt.makecpt(cmap="inferno", reverse=False, series=[0,swarm_hours])
                fig.plot(x=SWNE_down_distances[swarm_index], y=depths[swarm_index], color=hours.tolist(), cmap=True, style="c0.07c", pen="black", transparency=20)

            # Plot shallow-deep separator
            fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
            fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
            fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

            # Plot elevation profile
            fig.plot(x=SWNE_elev_distances, y=-0.001*SWNE["Elev(m)"], pen="1.5p,black")
            fig.text(x=0.5, y=-4.5, text="A", justify="LT", font="13p,black")
            fig.text(x=19.5, y=-4.5, text="B", justify="RT", font="13p,black")

        # Move plot origin to plot cross-section plot
        fig.shift_origin(xshift="4.5c")

        with fig.subplot(nrows=1, ncols=1, figsize=("4c", "10c")):

            # Create basemap with correct dimensions
            fig.basemap(region=REGION_NWSE, projection="X4c/-10c", frame=['xa5f1+l"Distance (km)"', 'ya5f1+l"Depth (km)"', "wsNE"],
                        panel=[0, 0])

            # # Plot landfill
            # for i in range(1, len(NWSE)):
            #     delta = (NWSE_elev_distances[i] - NWSE_elev_distances[i - 1]) / (REGION_NWSE[1]-REGION_NWSE[0]) * 4  # inter-measurement width
            #     height_km = (MAX_DEPTH + 0.001 * NWSE["Elev(m)"][i])  # depth - negative elevation
            #     height_cm = height_km / (MAX_DEPTH+5) * 9.9  # ratio * figure height
            #     midpoint = MAX_DEPTH - 0.5 * height_km
            #     data = [[NWSE_elev_distances[i], midpoint, delta, height_cm]]
            #     fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

            # Plot earthquakes
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=NWSE_down_distances[regular_index], y=depths[regular_index], color="lightgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            if len(longitudes[swarm_index]) > 0:
                pygmt.makecpt(cmap="inferno", reverse=False, series=[0,swarm_hours])
                fig.plot(x=NWSE_down_distances[swarm_index], y=depths[swarm_index], color=hours.tolist(), cmap=True, style="c0.07c", pen="black", transparency=20)

            # Plot shallow-deep separator
            fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
            fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
            fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

            # Plot elevation profile
            fig.plot(x=NWSE_elev_distances, y=-0.001 * NWSE["Elev(m)"], pen="1.5p,black")
            fig.text(x=0.5, y=-4.5, text="C", justify="LT", font="13p,black")
            fig.text(x=19.5, y=-4.5, text="D", justify="RT", font="13p,black")

        fig.savefig(swarm_framename)
