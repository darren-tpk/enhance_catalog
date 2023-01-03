with open('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/redoubt2/relocate_catalog/node_bootstrap_results.txt') as f:
    lines = f.readlines()

evid = []
rlat = []
rlon = []
rdep = []
lats = []
lons = []
deps = []
maderr_h = []
maderr_z = []
maderr_t = []
stdev_h = []
stdev_z = []
stdev_t = []

for i, line in enumerate(lines):
    # Skip first line
    if i == 0:
        continue
    fields = line.split()
    # Check if it is a header
    if len(fields) == 19:
        evid.append(int(fields[1]))
        rlat.append(float(fields[2]))
        rlon.append(float(fields[3]))
        rdep.append(float(fields[4]))
        lats.append([])
        lons.append([])
        deps.append([])
        maderr_h.append(float(fields[13]))
        maderr_z.append(float(fields[14]))
        maderr_t.append(float(fields[15]))
        stdev_h.append(float(fields[16]))
        stdev_z.append(float(fields[17]))
        stdev_t.append(float(fields[18]))
    # Else it is an interation
    elif len(fields) == 5:
        lats[-1].append(float(fields[0]))
        lons[-1].append(float(fields[1]))
        deps[-1].append(float(fields[2]))

rlat_x = [rlat[i] for i in range(len(rlat)) if maderr_h[i] == -1]
rlon_x = [rlon[i] for i in range(len(rlon)) if maderr_h[i] == -1]
rdep_x = [rdep[i] for i in range(len(rdep)) if maderr_h[i] == -1]
rlat_t = [rlat[i] for i in range(len(rlat)) if maderr_h[i] != -1]
rlon_t = [rlon[i] for i in range(len(rlon)) if maderr_h[i] != -1]
rdep_t = [rdep[i] for i in range(len(rdep)) if maderr_h[i] != -1]
maderr_h_t = [maderr_h[i] for i in range(len(rdep)) if maderr_h[i] != -1]
maderr_z_t = [maderr_z[i] for i in range(len(rdep)) if maderr_h[i] != -1]
maderr_t_t = [maderr_t[i] for i in range(len(rdep)) if maderr_h[i] != -1]


from matplotlib import pyplot as plt
# plt.hist(stdev_h,100,color='navy')
# plt.xlim([0,0.9])
# plt.title('horizontal standard deviation (m), N/A: ' + str(stdev_h.count(-1)))
# plt.show()
# plt.hist(stdev_z,100,color='maroon')
# plt.xlim([0,1.2])
# plt.title('vertical standard deviation (m), N/A: ' + str(stdev_z.count(-1)))
# plt.show()
# plt.hist(stdev_t,100,color='darkgreen')
# plt.xlim([-0.05,0.15])
# plt.title('temporal standard deviation (s), N/A: ' + str(stdev_t.count(-1)))
# plt.show()

plt.hist(maderr_h_t,50,color='navy')
plt.xlim([0,0.9])
plt.title('horizontal median abs deviation (m), N/A: ' + str(maderr_h.count(-1)))
plt.show()
plt.hist(maderr_z_t,50,color='maroon')
plt.xlim([0,1.2])
plt.title('vertical median abs deviation (m), N/A: ' + str(maderr_z.count(-1)))
plt.show()
plt.hist(maderr_t_t,250,color='darkgreen')
plt.xlim([0,0.1])
plt.title('temporal median abs deviation (s), N/A: ' + str(maderr_t.count(-1)))
plt.show()

################### plot redoubt
import numpy as np
import pandas as pd
import pygmt
from matplotlib import pyplot as plt
from toolbox import raster2array, reader
from obspy import UTCDateTime, Catalog
from pandas import read_csv
import geopy.distance
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
elev_profile_dir = main_dir + 'data/dem/'
grid_filepath = elev_profile_dir + 'redoubt_hillshade.tif'
WE_profile_filename = elev_profile_dir + 'we_profile_redoubt.csv'
NS_profile_filename = elev_profile_dir + 'ns_profile_redoubt.csv'
VOLC_LAT = 60.490      # dome: 60.490, -152.7626
VOLC_LON = -152.7626    # summit: 60.4852, -152.7438
MAX_DEPTH = 15 # TRY 40  # km
REGION = [-153.1438, -152.3438, 60.2852, 60.6800]
REGION_CLEAN = REGION  # edit region to fix spill at top
REGION_WE = [-153.1438, -152.3438, -5, MAX_DEPTH]
REGION_NS = [-5, MAX_DEPTH, 60.2852, 60.6800]
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

# Calcuulate the multiplier for errors (cm per km)
h_multiplier = 10 / geopy.distance.GeodesicDistance((VOLC_LAT, REGION_WE[0]), (VOLC_LAT, REGION_WE[1])).km
z_multiplier = 5 / 20

# Initialize figure
fig = pygmt.Figure()

# Craft title
with fig.subplot(nrows=1, ncols=1, figsize=("17c", "16.5c")):
    plot_title = '+t\"MAD of Hypocenter Locations (N=%d)\"' % len(rlat)
    with pygmt.config(MAP_FRAME_PEN='white', MAP_TICK_PEN_PRIMARY='white', FONT_ANNOT='white'):
        fig.basemap(region=[0, 17, 0, 16.5], projection="X17.5c/16.5c", frame=plot_title)

# Bottom plot: E-W cross section
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "5c"), autolabel="c)"):
    # Create basemap with correct dimensions
    fig.basemap(region=REGION_WE, projection="X10c/-5c", frame=["xa0.2f0.05+lLongitude", "yaf+lDepth", "WSne"],
                panel=[0, 0])

    # # Plot landfill
    # for i in range(1, len(WE)):
    #     delta = abs(WE["Lon"][i] - WE["Lon"][i - 1]) / (REGION_WE[1] - REGION_WE[0]) * 10  # inter-measurement width
    #     height_km = (MAX_DEPTH + 0.001 * WE["Elev(m)"][i])  # depth - negative elevation
    #     height_cm = height_km / (MAX_DEPTH+5) * 4.9  # ratio * figure height
    #     midpoint = MAX_DEPTH - 0.5 * height_km
    #     data = [[WE["Lon"][i], midpoint, delta, height_cm]]
    #     fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")

    # Plot earthquakes that have no bootstrap error(?)
    fig.plot(x=rlon_x, y=rdep_x, style="c0.07c", pen="gray15", transparency=30)

    # Now plot error crosshairs for earthquakes that have bootstrap error
    for k in range(len(rlon_t)):
        lon_spec = rlon_t[k]
        dep_spec = rdep_t[k]
        errh_spec = "%.3f" % (maderr_h_t[k] * h_multiplier)
        errz_spec = "%.3f" % (maderr_z_t[k] * z_multiplier)
        style_r1 = "r0.001/" + errz_spec
        style_r2 = "r" + errh_spec + "/0.001"
        fig.plot(x=lon_spec, y=dep_spec, style=style_r1, pen="0.01p,red4")
        fig.plot(x=lon_spec, y=dep_spec, style=style_r2, pen="0.01p,red4")

    # Plot elevation profile
    fig.plot(x=WE["Lon"], y=-0.001 * WE["Elev(m)"], pen="1.5p,black")

    # Plot volcano lon/elev
    fig.plot(x=VOLC_LON, y=-2.8, color="red", style="t0.3c", pen="black")

# Move plot origin to plot top-left plot
fig.shift_origin(yshift="h+0.5c")

# Left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):
    # Create basemap with correct dimensions
    pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
    fig.basemap(region=REGION_CLEAN, projection="X10c/10c",
                frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
    fig.grdimage(grid=grid_filepath, shading=True, cmap="geo", transparency=40)

    # Plot earthquakes that have no bootstrap error(?)
    fig.plot(x=rlon_x, y=rlat_x, style="c0.07c", pen="gray15", transparency=30)

    # Now plot error crosshairs for earthquakes that have bootstrap error
    for k in range(len(rlon_t)):
        lon_spec = rlon_t[k]
        lat_spec = rlat_t[k]
        errh_spec = "%.3f" % (maderr_h_t[k] * h_multiplier)
        style_r1 = "r0.001/" + errh_spec
        style_r2 = "r" + errh_spec + "/0.001"
        fig.plot(x=lon_spec, y=lat_spec, style=style_r1, pen="0.01p,red4")
        fig.plot(x=lon_spec, y=lat_spec, style=style_r2, pen="0.01p,red4")

    # Plot volcano lat/lon
    fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

# Move plot origin to plot cross-section plot
fig.shift_origin(xshift="10.5c", yshift="0c")

# Top-right plot: N-S cross section
with fig.subplot(nrows=1, ncols=1, figsize=("5c", "10c"), autolabel="b)"):
    # Create basemap with correct dimensions
    fig.basemap(region=REGION_NS, projection="X5c/10c", frame=["xaf+lDepth", "ya0.1f0.025+lLatitude", "wsNE"],
                panel=[0, 0])

    # # Plot landfill
    # for i in range(1, len(NS)):
    #     delta = abs(NS["Lat"][i] - NS["Lat"][i - 1]) / (REGION_NS[1] - REGION_NS[0]) * 10  # inter-measurement width
    #     height_km = (MAX_DEPTH + 0.001 * NS["Elev(m)"][i])  # depth - negative elevation
    #     height_cm = height_km / (MAX_DEPTH+5) * 4.9  # ratio * figure height
    #     midpoint = MAX_DEPTH - 0.5 * height_km
    #     data = [[midpoint, NS["Lat"][i], height_cm, delta]]
    #     fig.plot(data=data, style="r", color="gray90", pen="1.5p,gray90")

    # Plot earthquakes
    fig.plot(x=rdep_x, y=rlat_x, style="c0.07c", pen="gray15", transparency=30)

    # Now plot error crosshairs for earthquakes that have bootstrap error
    for k in range(len(rlon_t)):
        dep_spec = rdep_t[k]
        lat_spec = rlat_t[k]
        errz_spec = "%.3f" % (maderr_z_t[k] * z_multiplier)
        errh_spec = "%.3f" % (maderr_h_t[k] * h_multiplier)
        style_r1 = "r0.001/" + errh_spec
        style_r2 = "r" + errz_spec + "/0.001"
        fig.plot(x=dep_spec, y=lat_spec, style=style_r1, pen="0.01p,red4")
        fig.plot(x=dep_spec, y=lat_spec, style=style_r2, pen="0.01p,red4")

    # Plot elevation profile
    fig.plot(x=-0.001 * NS["Elev(m)"], y=NS["Lat"], pen="1.5p,black")

    # Plot volcano elev/lat
    fig.plot(x=-2.8, y=VOLC_LAT, color="red", style="t0.3c", pen="black")

fig.show(method="external")
fig.savefig('/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/redoubt2/relocate_catalog/MAD_layer_crosshairs.png',dpi=720)