###################################

# Plot the stations but hard code it

import numpy as np
import pandas as pd
import pygmt
from matplotlib import pyplot as plt
from toolbox import raster2array, reader
from obspy import UTCDateTime, Catalog
from pandas import read_csv

# Directories
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
elev_profile_dir = main_dir + 'data/dem/'
# Augustine
PAV_LAT = 55.4173
PAV_LON = -161.8937
PAV_REGION = [-162.123823, -161.613437, 55.281719, 55.552062]
# Stations
PV6A = (55.507,-161.9714,"infrasound","PV6A")
PN7A = (55.4329,-161.9973,"infrasound","PN7A")
PS4A = (55.346,-161.8567,"infrasound","PS4A")
PVV = (55.3732,-161.7919,"infrasound","PVV")
HAG = (55.317,-161.9045,"infrasound","HAG")
PS1A = (55.4201,-161.7437,"infrasound","PS1A")

fig = pygmt.Figure()

# Left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("8c","8c")):

    # Create basemap with correct dimensions
    grid = pygmt.datasets.load_earth_relief(resolution="01s", region=PAV_REGION)
    fig.basemap(region=PAV_REGION, projection="X8c/8c",
                frame=["x+lLongitude", "y+lLatitude", 'WsNe'], panel=[0, 0])
    fig.grdimage(grid=grid, shading=True, cmap="geo", transparency=10)
    # Plot volcano lat/lon
    fig.plot(x=PAV_LON, y=PAV_LAT, color="yellow", style="t0.6c", pen="2p,black")
    fig.plot(x=-161.8544, y=55.4569, color="white", style="t0.6c", pen="2p,black")
    fig.text(x=PAV_LON, y=PAV_LAT, text='Pavlof', justify="RM", font="13p,Helvetica-Bold,black", offset="0.6/0.6")
    fig.text(x=-161.8544, y=55.4569, text='Pavlof Sister', justify="RM", font="13p,Helvetica-Bold,black", offset="1.4/0.6")

    # Plot stations
    fig.plot(x=PV6A[1], y=PV6A[0], style="t0.45c", color="black", pen="1p,black")
    fig.plot(x=PN7A[1], y=PN7A[0], style="s0.45c", color="black", pen="1p,black")
    fig.plot(x=PS4A[1], y=PS4A[0], style="s0.45c", color="black", pen="1p,black")
    fig.plot(x= PVV[1], y= PVV[0], style="t0.45c", color="black", pen="1p,black")
    fig.plot(x= HAG[1], y= HAG[0], style="t0.45c", color="black", pen="1p,black")
    fig.plot(x=PS1A[1], y=PS1A[0], style="s0.45c", color="black", pen="1p,black")
    fig.text(x=PV6A[1], y=PV6A[0], text=PV6A[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")
    fig.text(x=PN7A[1], y=PN7A[0], text=PN7A[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")
    fig.text(x=PS4A[1], y=PS4A[0], text=PS4A[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")
    fig.text(x= PVV[1], y= PVV[0], text= PVV[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")
    fig.text(x= HAG[1], y= HAG[0], text= HAG[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")
    fig.text(x=PS1A[1], y=PS1A[0], text=PS1A[3], justify="RM", font="13p,Helvetica-Bold,black", offset="0.45/0.45")

    fig.plot(x=-161.7, y=55.297, style="r3.1/1.0", color="white", pen="0p", transparency=25)
    fig.plot(x=-161.77, y=55.303, style="t0.45c", color="black", pen="1p,black")
    fig.plot(x=-161.77, y=55.290, style="s0.45c", color="black", pen="1p,black")
    fig.text(x=-161.75, y=55.305,  text="Short Period", justify="LM", font="10p,Helvetica,black")
    fig.text(x=-161.75, y=55.290, text="Broadband", justify="LM", font="10p,Helvetica,black")

fig.savefig('/Users/darrentpk/Desktop/pavlof_dem.png',transparent=True)
fig.show(method="external")