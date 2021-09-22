#%% PLOT HYPOCENTERS

# This script loads in a catalog of choice and plots the hypocenters in 3 panels
# (a N-S cross section, an E-W cross section, and a top-down view of the hypocenters)
# The user will have to specify and tweak lat and lon limits, and region settings
# For detailed description of plotting syntax, refer to PyGMT documentation

# Import all dependencies
from phase_processing.read_hypoddpha import read_hypoddpha
import numpy as np
import pygmt
import os, glob
from PIL import Image
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
cores_filename = main_dir + 'output/mammoth/convert_redpy/core_catalog_picked.xml'
unmatched_PEC_events_filepath = main_dir + 'output/mammoth/convert_redpy/unmatched_PEC_events.xml'
relocated_catalog_filepath = main_dir + 'output/mammoth/relocate_catalog/relocated_catalog.xml'
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
relocated_catalog = reader(relocated_catalog_filepath)

#%% Extract information from catalog, filtered by max depth

# Choose catalog
catalog = relocated_catalog

# Extract hypocenter and time information
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
times = np.array([np.datetime64(event.origins[0].time) for event in catalog])
days = ((utctimes - START_TIME) / 86400).astype(int)

# Filter all arrays by depth
valid_index = np.where(depths<MAX_DEPTH)
latitudes = latitudes[valid_index]
longitudes = longitudes[valid_index]
depths = depths[valid_index]
times = times[valid_index]
days = days[valid_index]

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
ec_start = UTCDateTime(2012,10,1)
ec_end = UTCDateTime(2013,2,1)
time_step = 4*60*60 # 4 hours
ec_time_blocks = int((ec_end - ec_start) / time_step)

for i in range(ec_time_blocks+1):

    # tidy up some things
    swarm_reference = ec_start + (i*time_step)
    # plot only earthquakes that have happened

    # Choose duration of interest
    if swarm_focus:

        # Define some variables
        swarm_hours = 40  # hours
        swarm_destination = "/Users/darrentpk/Desktop/frames/"
        swarm_framename = swarm_destination + swarm_reference.strftime("%Y_%m_%d_%H_%M_%S") +'.png'

        # Determine swarm start and swarm end
        swarm_start = swarm_reference
        swarm_end = swarm_reference + (swarm_hours * 3600)

        # Find out which hypocenters to color
        swarm_index = np.where((utctimes>swarm_start) & (utctimes<swarm_end))
        regular_index = np.where((utctimes<swarm_start))

        # Now compute hours from swarm start
        hours = (utctimes[swarm_index] - swarm_reference) / (3600)

    # Initialize figure
    fig = pygmt.Figure()

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("8.5c","8.5c")):

        # Create basemap with correct dimensions
        grid = pygmt.datasets.load_earth_relief(resolution="01s", region=REGION_CLEAN)
        fig.basemap(region=REGION_CLEAN, projection="X8.5c/8.5c",
                    frame=["x+lLongitude", "y+lLatitude", "WsNe"], panel=[0, 0])
        fig.grdimage(grid=grid, shading=False, cmap="geo")

        # Plot earthquakes
        if not swarm_focus:
            pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
            fig.plot(x=longitudes, y=latitudes, color=times, cmap=True, style="c0.07c", pen="black", transparency=20)
            with pygmt.config(FONT_ANNOT_PRIMARY="7p,Helvetica,black"):
                fig.colorbar(position="JBL+o-8.5c/0.5c+w8.5c/0.4c+h",frame="xa1Of")
        else:
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=longitudes[regular_index], y=latitudes[regular_index], color="darkgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            pygmt.makecpt(cmap="inferno", reverse=True, series=[0,swarm_hours])
            if len(longitudes[swarm_index]) > 0:
                fig.plot(x=longitudes[swarm_index], y=latitudes[swarm_index], color=hours.tolist(), cmap=True, style="c0.07c", pen="black", transparency=20)
            fig.colorbar(position="JBL+o-8.5c/0.5c+w8.5c/0.4c+h",frame='xa10f5+l"Hours from time stamp"')
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
        if not swarm_focus:
            pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
            fig.plot(x=SWNE_down_distances, y=depths, style="c0.07c", color=times, cmap=True, pen="black", transparency=20)
        else:
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=SWNE_down_distances[regular_index], y=depths[regular_index], color="darkgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            if len(longitudes[swarm_index]) > 0:
                pygmt.makecpt(cmap="inferno", reverse=True, series=[0,swarm_hours])
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
        if not swarm_focus:
            pygmt.makecpt(cmap="viridis", series="2012-10-01T/2013-02-01T/1d")
            fig.plot(x=NWSE_down_distances, y=depths, style="c0.07c", color=times, cmap=True, pen="black", transparency=20)
        else:
            if len(longitudes[regular_index]) > 0:
                fig.plot(x=NWSE_down_distances[regular_index], y=depths[regular_index], color="darkgray", cmap=False, style="c0.07c", pen="black", transparency=50)
            if len(longitudes[swarm_index]) > 0:
                pygmt.makecpt(cmap="inferno", reverse=True, series=[0,swarm_hours])
                fig.plot(x=NWSE_down_distances[swarm_index], y=depths[swarm_index], color=hours.tolist(), cmap=True, style="c0.07c", pen="black", transparency=20)

        # Plot shallow-deep separator
        fig.plot(x=[0, 20], y=[8, 8], pen="1p,black,-")
        fig.text(x=0.5, y=7.5, text="Shallow", justify="LM", font="6.5p,black")
        fig.text(x=0.5, y=8.5, text="Deep", justify="LM", font="6.5p,black")

        # Plot elevation profile
        fig.plot(x=NWSE_elev_distances, y=-0.001 * NWSE["Elev(m)"], pen="1.5p,black")
        fig.text(x=0.5, y=-4.5, text="C", justify="LT", font="13p,black")
        fig.text(x=19.5, y=-4.5, text="D", justify="RT", font="13p,black")

    fig.show(method="external")

    fig.savefig(swarm_framename)