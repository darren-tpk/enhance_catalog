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
from obspy import UTCDateTime
from pandas import read_csv

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
elev_profile_dir = main_dir + 'data/dem/'
EW_profile_filename = elev_profile_dir + 'ew_profile.csv'
NS_profile_filename = elev_profile_dir + 'ns_profile.csv'
PEC_events_dir = main_dir + 'data/avo/'
located_cores_filename = main_dir + 'output/greatsitkin2/convert_redpy/core_catalog_picked.xml'
unmatched_PEC_events_filename = main_dir + 'output/greatsitkin2/convert_redpy/unmatched_PEC_events.xml'
relocated_catalog_filename = main_dir + 'output/greatsitkin2/relocate_catalog/relocated_catalog.xml'
VOLC_LAT = 60.4852
VOLC_LON = -152.7438
MAX_DEPTH = 15  # km
REGION = [VOLC_LON - 0.4, VOLC_LON + 0.4,
          VOLC_LAT - 0.2, VOLC_LAT + 0.2]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_CLEAN[3] = 60.68  # hard coded
MAIN_WIDTH = 10
REGION_EW = [VOLC_LON - 0.4, VOLC_LON + 0.4, -5, MAX_DEPTH]
REGION_NS = [-5, MAX_DEPTH, VOLC_LAT - 0.2, 60.68]  # hard coded
START_TIME = UTCDateTime(2009,1,1,0,0,0)  # reference time for "days" colorbar

#%% Define functions

# Function to craft a PyGMT projection quote
def projection(region, width, cm=False):
    if cm:
        unit = 'c'
    else:
        unit = 'i'
    return 'S{}/90/{}{}'.format(np.mean(region[:2]), width, unit)

#%% Read all catalogs

# (1) Pre-existing catalog
PEC_hypoi = PEC_events_dir + 'redoubt_20090101_20090501_hypoi.txt'
PEC_hypoddpha = PEC_events_dir + 'redoubt_20090101_20090501_hypoddpha.txt'
PEC_events = read_hypoddpha(PEC_hypoi, PEC_hypoddpha, channel_convention=True)

# (2) All located templates
located_cores = reader(located_cores_filename)
unmatched_PEC_events = reader(unmatched_PEC_events_filename)
located_templates = located_cores + unmatched_PEC_events

# (3) GrowClust relocated catalog
relocated_catalog = reader(relocated_catalog_filename)

#%% Extract information from catalog, filtered by max depth

# Choose catalog
catalog = located_templates

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

#%% Plot hypocenters using PyGMT

# Load EW and NS elevation profiles (from https://apps.nationalmap.gov/elevation/)
EW = read_csv(EW_profile_filename)
NS = read_csv(NS_profile_filename)

# Initialize figure
fig = pygmt.Figure()

# Bottom plot: E-W cross section
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "4c"), autolabel="c)"):

    # Create basemap with correct dimensions
    fig.basemap(region=REGION_EW, projection="X10c/-4c", frame=["xa0.2f0.05+lLongitude","yaf+lDepth", "WSne"], panel=[0, 0])

    # Plot landfill
    for i in range(1, len(EW)):
        delta = abs(EW["Lon"][i] - EW["Lon"][i - 1]) / (REGION_EW[1]-REGION_EW[0]) * 10  # inter-measurement width
        height_km = (15 + 0.001 * EW["Elev(m)"][i])  # depth - negative elevation
        height_cm = height_km / 20 * 3.9  # ratio * figure height
        midpoint = 15 - 0.5 * height_km
        data = [[EW["Lon"][i], midpoint, delta, height_cm]]
        fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

    # Plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2009-01-01T/2009-05-01T/1d")
    fig.plot(x=longitudes, y=depths, color=times, cmap=True, style="c0.07c", pen="black", transparency=20)
    fig.colorbar(position="JTR+o1c/-4c+w4c/0.4c+v",frame="xa1Of+lDatetime")

    # Plot elevation profile
    fig.plot(x=EW["Lon"], y=-0.001*EW["Elev(m)"], pen="1.5p,black")

# Move plot origin to plot top-left plot
fig.shift_origin(yshift="h+1c")

# Top-left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):

    # Create basemap with correct dimensions
    fig.grdimage("@earth_relief_15s",region=REGION_CLEAN,
                 projection=projection(REGION_CLEAN, MAIN_WIDTH, cm=True),
                 shading=True, t=30, cmap="geo")
    fig.basemap(region=REGION_CLEAN, projection="X10c/10c", frame=["xa0.2f0.05+lLongitude","ya0.1f0.025+lLatitude", "WsNe"], panel=[0, 0])

    # Plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2009-01-01T/2009-05-01T/1d")
    fig.plot(x=longitudes, y=latitudes, color=times, cmap=True, style="c0.07c", pen="black", transparency=20)

# Move plot origin to plot top-right plot
fig.shift_origin(xshift="11c",yshift="0.04c")

# Top-right plot: N-S cross section
with fig.subplot(nrows=1, ncols=1, figsize=("4c", "9.95c"), autolabel="b)"):

    # Create basemap with correct dimensions
    fig.basemap(region=REGION_NS, projection="X4c/9.95c", frame=["xaf+lDepth","ya0.1f0.025+lLatitude", "wsNE"], panel=[0, 0])

    # Plot landfill
    for i in range(1, len(NS)):
        delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[1]-REGION_NS[0]) * 10 # inter-measurement width
        height_km = (15 + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
        height_cm = height_km / 20 * 3.9  # ratio * figure height
        midpoint = 15 - 0.5 * height_km
        data = [[midpoint, NS["Lat"][i], height_cm, delta]]
        fig.plot(data=data, style="r", color="gray90", pen="1.5p,gray90")

    # Plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2009-01-01T/2009-05-01T/1d")
    fig.plot(x=depths, y=latitudes, style="c0.07c", color=times, cmap=True, pen="black", transparency=20)

    # Plot elevation profile
    fig.plot(x=-0.001*NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

# Show figure
fig.show(method='external')
