#%% PLOT REDOUBT TEMPORAL

# This script loads in a Redoubt catalog of choice and plots the temporal progression
# of FI and depth of events where possible.

# Import all dependencies
from obspy import UTCDateTime
import numpy as np
from toolbox import reader
from pandas import read_csv
import matplotlib.pyplot as plt

# Define general paths
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
elev_profile_dir = main_dir + 'data/dem/'
WE_profile_filename = elev_profile_dir + 'we_profile_redoubt.csv'
NS_profile_filename = elev_profile_dir + 'ns_profile_redoubt.csv'
add_labels = True

# Load color code labels
color_code_table = read_csv(main_dir + 'data/supporting/rdt_color_codes.csv')
color_codes = []
color_code_times = []
for i in range(len(color_code_table)):
    color_code_time = UTCDateTime(color_code_table.year[i],color_code_table.month[i],color_code_table.day[i],color_code_table.hh[i],color_code_table.mm[i]) \
                      - (color_code_table.utc_shift[i] * 3600)
    color_codes.append(color_code_table.color[i])
    color_code_times.append(color_code_time)


# Load explosion times
explosion_table = read_csv(main_dir + 'data/supporting/rdt_explosion_list.csv')
explosion_label = []
explosion_times = []
for i in range(len(explosion_table)):
    explosion_label.append(str(explosion_table.event_number[i]))
    explosion_times.append(UTCDateTime(explosion_table.year[i],explosion_table.month[i],explosion_table.day[i],explosion_table.hh[i],explosion_table.mm[i]))

# load catalogs
evlist_filepath = main_dir + 'output/redoubt2/scan_data/party_catalog_full_reFI.xml'
evlist = reader(evlist_filepath)
# relocated_catalog_filepath = main_dir + 'output/redoubt2/relocate_catalog/relocated_catalog_reFImag.xml'
# relocated_catalog = reader(relocated_catalog_filepath)
original_filepath = main_dir + 'output/redoubt2/scan_data/PEC_events_reFI.xml'
original = reader(original_filepath)

# event list
utctimes = np.array([UTCDateTime(event.resource_id.id.split('_')[-1]) for event in evlist])
FI_values = []
for event in evlist:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)

original_utctimes = np.array([event.origins[0].time for event in original])
original_FI_values = []
for event in original:
    if event.comments[-1].text.split('=')[1] == 'None':
        original_FI_values.append(np.nan)
    else:
        original_FI_values.append(float(event.comments[-1].text.split('=')[1]))
original_FI_values = np.array(original_FI_values)

# remove overlapping
remove_index = []
for i, original_utctime in enumerate(original_utctimes):
    if np.min(np.abs(utctimes-original_utctime))<1:
        remove_index.append(np.argmin(np.abs(utctimes-original_utctime)))
    print(i,'/',len(original_utctimes))
sorted_remove_index = list(np.sort(remove_index))
utctimes = list(utctimes)
FI_values = list(FI_values)
for i in reversed(sorted_remove_index):
    utctimes.pop(i)
    FI_values.pop(i)
utctimes = np.array(utctimes)
FI_values = np.array(FI_values)

# combine
utctimes = np.concatenate((utctimes,original_utctimes))
FI_values = np.concatenate((FI_values,original_FI_values))

# full duration
base_time = UTCDateTime(2008,5,1)
end_time = UTCDateTime(2009,9,1)
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
zoom1_base_time = UTCDateTime(2009,3,14)
zoom1_end_time = UTCDateTime(2009,4,5)
zoom2_base_time = UTCDateTime(2009,5,2)
zoom2_end_time = UTCDateTime(2009,5,9)
tick_times = [UTCDateTime(2008, 5, 1), UTCDateTime(2008, 6, 1),
              UTCDateTime(2008, 7, 1), UTCDateTime(2008, 8, 1),
              UTCDateTime(2008, 9, 1), UTCDateTime(2008, 10, 1),
              UTCDateTime(2008, 11, 1), UTCDateTime(2008, 12, 1),
              UTCDateTime(2009, 1, 1), UTCDateTime(2009, 2, 1),
              UTCDateTime(2009, 3, 1), UTCDateTime(2009, 4, 1),
              UTCDateTime(2009, 5, 1), UTCDateTime(2009, 6, 1),
              UTCDateTime(2009, 7, 1), UTCDateTime(2009, 8, 1),
              UTCDateTime(2009, 9, 1)]
fig, ax = plt.subplots(figsize=(22, 4))
ax.hist2d(days, FI_values,bins=[int(num_days),280],range=[[0, num_days], [-2,0.8]],cmap='hot',cmin=1)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %y') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI with time (Temporally Complete Event List, N=%d)' % len(FI_values),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
    ax.axhspan(ymin=-2, ymax=-1.9, xmin=(zoom1_base_time-base_time)/(86400*num_days), xmax=(zoom1_end_time-base_time)/(86400*num_days),color='blue', alpha=0.5)
    ax.axhspan(ymin=-2, ymax=-1.9, xmin=(zoom2_base_time-base_time)/(86400*num_days), xmax=(zoom2_end_time-base_time)/(86400*num_days), color='blue', alpha=0.5)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/combined_2dhist.jpg')

base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(16, 4))
ax.hist2d(days, FI_values,bins=[int(num_days*12),280],range=[[0, num_days], [-2,0.8]],cmap='hot',cmin=1)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Temporally Complete Event List, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/combined_2dhist_ZOOM1.jpg')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(5, 4))
ax.hist2d(days, FI_values,bins=[int(num_days*24),280],range=[[0, num_days], [-2,0.8]],cmap='hot',cmin=1)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
#ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('(N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/combined_2dhist_ZOOM2.jpg')

#### now try pygmt

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
grid_filepath = elev_profile_dir + 'redoubt_hillshade.tif'
WE_profile_filename = elev_profile_dir + 'we_profile_redoubt.csv'
NS_profile_filename = elev_profile_dir + 'ns_profile_redoubt.csv'
PEC_events_filepath = main_dir + 'output/redoubt2/scan_data/PEC_events_reFI.xml'
cores_filename = main_dir + 'output/redoubt2/convert_redpy/core_catalog_picked.xml'
unmatched_PEC_events_filepath = main_dir + 'output/redoubt2/convert_redpy/unmatched_PEC_events.xml'
relocated_catalog_filepath = main_dir + 'output/redoubt2/relocate_catalog/relocated_catalog_reFImag.xml'
plot_output_dir = main_dir + 'output/redoubt2/relocate_catalog/'
VOLC_LAT = 60.4852
VOLC_LON = -152.7438
MAX_DEPTH = 15 # TRY 40  # km
REGION = [-153.1438, -152.3438, 60.2852, 60.6800]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_WE = [-153.1438, -152.3438, -5, MAX_DEPTH]
REGION_NS = [-5, MAX_DEPTH, 60.2852, 60.6800]
size_by_magnitude = True
plot_temporal = True
plot_FI = True
migration_track = False
swarm_focus = False

#%% Read all catalogs

# (1) Pre-existing catalog
PEC_events = reader(PEC_events_filepath)

# (2) All located templates
cores = reader(cores_filename)
unmatched_PEC_events = reader(unmatched_PEC_events_filepath)
templates = cores + unmatched_PEC_events
located_templates = Catalog()
for template in templates:
    if template.origins[0].latitude is not None:
        located_templates += template

# (3) GrowClust relocated catalog
relocated_catalog = reader(relocated_catalog_filepath)  # min_cc = 0.75
#%% Extract information from catalog, filtered by max depth

# Choose catalog
catalog = PEC_events
if catalog == PEC_events:
    plot_title = '+t\"Original Catalog (N=%d)\"' % len(catalog)
    save_filepath = plot_output_dir + 'original_catalog_N' + str(len(PEC_events)) + '.jpg'
elif catalog == located_templates:
    plot_title = '+t\"Located Templates (N=%d)\"' % len(catalog)
    save_filepath = plot_output_dir + 'located_templates_N' + str(len(located_templates)) + '.jpg'
elif catalog == relocated_catalog:
    plot_title = '+t\"Relocated Catalog (N=%d)\"' % len(catalog)
    save_filepath = plot_output_dir + 'relocated_catalog_N' + str(len(relocated_catalog)) + '.jpg'

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

# Prepare cross-section data for SW-NE cross section
W = [60.4852, -153.1000]
E = [60.4852, -152.4000]
WE_x = [W[1], E[1]]
WE_xvec = np.linspace(W[1], E[1], 600)
WE_y = [W[0], E[0]]
WE_yvec = np.linspace(W[0], E[0], 600)

# Prepare cross-section data for NW-SE cross section
N = [60.66000, -152.7438]
S = [60.3000, -152.7438]
NS_x = [N[1], S[1]]
NS_xvec = np.linspace(N[1], S[1],600)
NS_y = [N[0], S[0]]
NS_yvec = np.linspace(N[0], S[0],600)
# Load EW and NS elevation profiles (from https://apps.nationalmap.gov/elevation/)
WE = read_csv(WE_profile_filename)
NS = read_csv(NS_profile_filename)

# Determine indices of events with calculated FI and events without calculated FI
valid_FI_index = np.where(~np.isnan(FI_values))
invalid_FI_index = np.where(np.isnan(FI_values))

# Normalize the magnitudes based on reasonable size values:
if size_by_magnitude:
    # construct reference sizes
    ref_magnitudes = np.array(range(0,401))
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

def pygmt_hist2d(X,Y,REGION,orientation):
    import matplotlib.pyplot as plt
    if orientation == 'h':
        out = plt.hist2d(X,Y,bins=[160,80],range=[[REGION[0],REGION[1]],[REGION[2],REGION[3]]])
    elif orientation == 'v':
        out = plt.hist2d(X,Y, bins=[80,160], range=[[REGION[0], REGION[1]], [REGION[2], REGION[3]]])
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
            logic = (X>xedge) & (X<(xedge+delta_x)) & (Y>yedge) & (Y<(yedge+delta_y))
            mean_fi = np.nanmedian(FI_values[logic])
            xgrid[i,j] = xedge + delta_x/2
            ygrid[i,j] = yedge + delta_y/2
            figrid[i,j] = mean_fi
    evvec = evgrid.flatten()[~np.isnan(figrid.flatten())]
    xvec = xgrid.flatten()[~np.isnan(figrid.flatten())]
    yvec = ygrid.flatten()[~np.isnan(figrid.flatten())]
    fivec = figrid.flatten()[~np.isnan(figrid.flatten())]
    return evvec,xvec,yvec,fivec

if plot_FI:

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
        evvec, xvec, yvec, fivec = pygmt_hist2d(longitudes, depths, REGION_WE,'h')
        pygmt.makecpt(cmap="polar", reverse=False, series=[-2, 0.5])
        num = 30
        fig.plot(x=xvec, y=yvec, color=fivec, cmap=True,style="s0.09c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec >= num], y=yvec[evvec >= num],style="d0.06c",pen='0.1p,black')
        fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa0.5f0.25+l"Frequency Index"')

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=40)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], size=sizes[invalid_FI_index],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=False, series=[-2, 0.5])
        if size_by_magnitude:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], size=sizes[valid_FI_index],
                     color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-153.025, y=60.325, style="r2.4/0.8", color="white", pen="0p", transparency=25)
            fig.plot(x=-153.1, y=60.32, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-153.05, y=60.32, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-153.00, y=60.32, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.95, y=60.32, color='gray', style="c0.24c", pen="black")
            fig.text(x=-153.1, y=60.328, text="M0", justify="CB", font="9p,black")
            fig.text(x=-153.05, y=60.328, text="M1", justify="CB", font="9p,black")
            fig.text(x=-153.00, y=60.328, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.95, y=60.328, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="b)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xaf+lDepth", "ya0.1f0.025+lLatitude", "wsNE"],
                    panel=[0, 0])

        evvec, xvec, yvec, fivec = pygmt_hist2d(depths, latitudes, REGION_NS,'v')
        pygmt.makecpt(cmap="polar", reverse=False, series=[-2, 0.5])
        num = 30
        fig.plot(x=xvec, y=yvec, color=fivec, cmap=True,style="s0.09c", pen='0.1p,darkgrey')
        fig.plot(x=xvec[evvec >= num], y=yvec[evvec >= num],style="d0.06c",pen='0.1p,black')
        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")
    fig.savefig('/Users/darrentpk/Desktop/mag_plots/original_catalog_reFI_hist2d.jpg')

#####

if plot_EV:

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
        evvec, xvec, yvec, fivec = pygmt_hist2d(longitudes, depths, REGION_WE,'h')
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 30])
        fig.plot(x=xvec, y=yvec, color=evvec, cmap=True,style="s0.09c", pen='0.1p,darkgrey')
        #fig.plot(x=xvec[evvec >= num], y=yvec[evvec >= num],style="d0.06c",pen='0.1p,black')
        fig.colorbar(position="JTR+o0.5c/-5c+w5c/0.4c+v", frame='xa5f1+l"Number of Events"')

        # Plot elevation profile
        fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Move plot origin to plot top-left plot
    fig.shift_origin(yshift="h+0.5c")

    # Left plot: Top-down view
    with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):

        # Create basemap with correct dimensions
        pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
        fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                    frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
        fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=40)

        # Plot earthquakes
        if size_by_magnitude:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], size=sizes[invalid_FI_index],
                     color="darkgray", cmap=False, style="x", pen="0.7p,black", transparency=20)
        else:
            fig.plot(x=longitudes[invalid_FI_index], y=latitudes[invalid_FI_index], color="darkgray", cmap=False,
                     style="x0.12c", pen="0.7p,black", transparency=20)
        pygmt.makecpt(cmap="polar", reverse=False, series=[-2, 0.5])
        if size_by_magnitude:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], size=sizes[valid_FI_index],
                     color=FI_values[valid_FI_index], cmap=True, style="c", pen="black", transparency=20)
        else:
            fig.plot(x=longitudes[valid_FI_index], y=latitudes[valid_FI_index], color=FI_values[valid_FI_index], cmap=True,
                     style="c0.07c", pen="black", transparency=20)

        # # Plot size indicator if sizing by mag:
        if size_by_magnitude:
            fig.plot(x=-153.025, y=60.325, style="r2.4/0.8", color="white", pen="0p", transparency=25)
            fig.plot(x=-153.1, y=60.32, color='gray', style="c0.06c", pen="black")
            fig.plot(x=-153.05, y=60.32, color='gray', style="c0.12c", pen="black")
            fig.plot(x=-153.00, y=60.32, color='gray', style="c0.18c", pen="black")
            fig.plot(x=-152.95, y=60.32, color='gray', style="c0.24c", pen="black")
            fig.text(x=-153.1, y=60.328, text="M0", justify="CB", font="9p,black")
            fig.text(x=-153.05, y=60.328, text="M1", justify="CB", font="9p,black")
            fig.text(x=-153.00, y=60.328, text="M2", justify="CB", font="9p,black")
            fig.text(x=-152.95, y=60.328, text="M3", justify="CB", font="9p,black")

        # Plot volcano lat/lon
        fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

    # Move plot origin to plot cross-section plot
    fig.shift_origin(xshift="10.5c", yshift="0c")

    # Top-right plot: N-S cross section
    with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="b)"):

        # Create basemap with correct dimensions
        fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xaf+lDepth", "ya0.1f0.025+lLatitude", "wsNE"],
                    panel=[0, 0])

        evvec, xvec, yvec, fivec = pygmt_hist2d(depths, latitudes, REGION_NS,'v')
        pygmt.makecpt(cmap="magma", reverse=True, series=[1, 30])
        fig.plot(x=xvec, y=yvec, color=evvec, cmap=True,style="s0.09c", pen='0.1p,darkgrey')
        #fig.plot(x=xvec[evvec >= num], y=yvec[evvec >= num],style="d0.06c",pen='0.1p,black')
        # Plot elevation profile
        fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    fig.show(method="external")
    fig.savefig('/Users/darrentpk/Desktop/mag_plots/original_catalog_numev_hist2d.jpg')
