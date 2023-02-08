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
AUG_LAT = 59.363
AUG_LON = -153.43
AUG_REGION = [-153.582, -153.338, 59.304, 59.426]
# Redoubt
RED_LAT = 60.490  # dome
RED_LON = -152.7626  # dome
RED_REGION = [-153.1438, -152.3438, 60.2852, 60.6800]
# Stations
AU10 = (59.3496,-153.3854,"campaign","AU10")
AU11 = (59.3596,-153.4803,"campaign","AU11")
AU12 = (59.3835,-153.4519,"campaign","AU12")
AU13 = (59.3464,-153.4341,"campaign","AU13")
AU14 = (59.3711,-153.3969,"campaign","AU14")
AU15 = (59.3507,-153.4856,"campaign","AU15")
AU20 = (59.3703,-153.3541,"campaign","AU20")
AUE = (59.3718,-153.3751,"permanent","AUE")
AUH = (59.3639,-153.4432,"permanent","AUH")
AUI = (59.3352,-153.4277,"permanent","AUI")
AUL = (59.3823,-153.4357,"permanent","AUL")
AUNW = (59.3782,-153.4768,"permanent","AUNW")
AUP = (59.3634,-153.4201,"permanent","AUP")
AUR = (59.3628,-153.4313,"permanent","AUR")
AUS = (59.3600,-153.4306,"permanent","AUS")
AUSE = (59.3414,-153.3975,"permanent","AUSE")
AUW = (59.3701,-153.4708,"permanent","AUW")
RDN = (60.5224,-152.7402,"permanent","RDN")
REF = (60.4886,-152.7039,"permanent","REF")
RSO = (60.4616,-152.7561,"permanent","RSO")
DFR = (60.5913,-152.6883,"permanent","DFR")
NCT = (60.5615,-152.9316,"permanent","NCT")
RDDF = (60.5912,-152.6883,"permanent","RDDF")
RDDR = (60.5843,-152.5887,"permanent","RDDR")
RDE = (60.5869,-152.5925,"permanent","RDE")
RDJH = (60.5905,-152.8058,"permanent","RDJH")
RDSO = (60.4536,-152.7430,"permanent","RDSO")
RDT = (60.5726,-152.4075,"permanent","RDT")
RDWB = (60.4874,-152.8423,"permanent","RDWB")
RED = (60.4196,-152.7742,"permanent","RED")
RD01 = (60.4885,-152.7033,"campaign","RD01")
RD02 = (60.5208,-152.7373,"campaign","RD02")
RD03 = (60.4705,-152.8203,"campaign","RD03")
RDW = (60.483,-152.810,"campaign","RDW")

fig = pygmt.Figure()

# Left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("7c","7c")):

    # Create basemap with correct dimensions
    pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
    fig.basemap(region=RED_REGION, projection="X7c/7c",
                frame=["xa0.2f0.1+lRedoubt", "ya0.1f0.05", 'WsNe'],
                panel=[0, 0])
    fig.grdimage(grid=elev_profile_dir + 'redoubt_hillshade.tif', shading=True, cmap="geo", transparency=50)

    # Plot volcano lat/lon
    fig.plot(x=RED_LON, y=RED_LAT, color="red", style="t0.3c", pen="black")

    # Plot stations

    # Campaign first
    fig.plot(x=RD01[1], y=RD01[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=RD02[1], y=RD02[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=RD03[1], y=RD03[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=RDW[1], y=RDW[0], style="s0.3c", pen="1p,darkred")
    fig.text(x=RD01[1], y=RD01[0], text=RD01[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="0.875/0.125")
    fig.text(x=RD02[1], y=RD02[0], text=RD02[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="0.875/0")
    fig.text(x=RD03[1], y=RD03[0], text=RD03[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    fig.text(x=RDW[1], y=RDW[0], text=RDW[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="0.3/0.3")

    # Permanent after
    fig.plot(x=RDN[1] , y=RDN[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=REF[1] , y=REF[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=RSO[1] , y=RSO[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=DFR[1] , y=DFR[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=NCT[1] , y=NCT[0] , style="i0.3c", pen="1p,black")
    # fig.plot(x=RDDF[1], y=RDDF[0], style="i0.3c", pen="1p,black")
    # fig.plot(x=RDDR[1], y=RDDR[0], style="i0.3c", pen="1p,black")
    fig.plot(x=RDE[1] , y=RDE[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=RDJH[1], y=RDJH[0], style="s0.3c", pen="1p,black")
    # fig.plot(x=RDSO[1], y=RDSO[0], style="i0.3c", pen="1p,black") # Not present in 2010
    fig.plot(x=RDT[1] , y=RDT[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=RDWB[1], y=RDWB[0], style="s0.3c", pen="1p,black")
    fig.plot(x=RED[1] , y=RED[0] , style="i0.3c", pen="1p,black")
    fig.text(x=RDN[1] , y=RDN[0] , text=RDN[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=REF[1] , y=REF[0] , text=REF[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="0.75/-0.125")
    fig.text(x=RSO[1] , y=RSO[0] , text=RSO[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.1/-0.1")
    fig.text(x=DFR[1] , y=DFR[0] , text=DFR[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=NCT[1] , y=NCT[0] , text=NCT[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    # fig.text(x=RDDF[1], y=RDDF[0], text=RDDF[3], justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    # fig.text(x=RDDR[1], y=RDDR[0], text=RDDR[3], justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/-0.1")
    fig.text(x=RDE[1] , y=RDE[0] , text=RDE[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.15/-0.1")
    fig.text(x=RDJH[1], y=RDJH[0], text=RDJH[3], justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    # fig.text(x=RDSO[1], y=RDSO[0], text=RDSO[3], justify="RM", font="8p,Helvetica-Bold,black", offset="-0.1/-0.2")
    fig.text(x=RDT[1] , y=RDT[0] , text=RDT[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=RDWB[1], y=RDWB[0], text=RDWB[3], justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=RED[1] , y=RED[0] , text=RED[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")

    # Plot legend
    fig.plot(x=-152.875, y=60.325, style="r4.45/1.1", color="white", pen="0p", transparency=25)
    fig.plot(x=-153.09, y=60.345, style="i0.3c", pen="1p,black")
    fig.plot(x=-153.09, y=60.325, style="s0.3c", pen="1p,black")
    fig.plot(x=-153.09, y=60.305, style="s0.3c", pen="1p,darkred")
    fig.text(x=-153.05, y=60.345, text="Permanent short-period stations", justify="LM", font="7p,Helvetica,black")
    fig.text(x=-153.05, y=60.325, text="Permanent broadband stations", justify="LM", font="7p,Helvetica,black")
    fig.text(x=-153.05, y=60.305, text="Temporary broadband stations", justify="LM", font="7p,Helvetica,black")
    fig.plot(x=[-152.571,-152.39], y=[60.325,60.325], pen="2p,black")
    fig.text(x=(-152.571-152.39)/2, y=60.335, text="10 km", justify="CB", font="8.5p,Helvetica,black")

# Move plot origin to plot cross-section plot
fig.shift_origin(xshift="8.5c",yshift="0c")

# Left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("7c","7c")):

    # Create basemap with correct dimensions
    pygmt.config(GMT_DATA_SERVER="https://oceania.generic-mapping-tools.org")
    fig.basemap(region=AUG_REGION, projection="X7c/7c",
                frame=["xa0.05f0.025+lAugustine", "ya0.04f0.02", 'WsNe'], panel=[0, 0])
    fig.grdimage(grid=elev_profile_dir + 'augustine_hillshade_island.tif', shading=True, cmap="geo", transparency=50)

    # Plot volcano lat/lon
    fig.plot(x=AUG_LON, y=AUG_LAT, color="red", style="t0.3c", pen="black")

    # Plot stations

    # Plot campaign first
    fig.plot(x=AU10[1], y=AU10[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=AU11[1], y=AU11[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=AU12[1], y=AU12[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=AU13[1], y=AU13[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=AU14[1], y=AU14[0], style="s0.3c", pen="1p,darkred")
    fig.plot(x=AU15[1], y=AU15[0], style="s0.3c", pen="1p,darkred")
    # fig.plot(x=AU20[1], y=AU20[0], style="s0.3c", pen="1p,darkred")
    fig.text(x=AU10[1], y=AU10[0], text=AU10[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    fig.text(x=AU11[1], y=AU11[0], text=AU11[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    fig.text(x=AU12[1], y=AU12[0], text=AU12[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    fig.text(x=AU13[1], y=AU13[0], text=AU13[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    fig.text(x=AU14[1], y=AU14[0], text=AU14[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.175/0.05")
    fig.text(x=AU15[1], y=AU15[0], text=AU15[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="-0.2/0")
    # fig.text(x=AU20[1], y=AU20[0], text=AU20[3], justify="RM", font="8p,Helvetica-Bold,darkred", offset="0.2/-0.25")

    # Plot permanent next

    fig.plot(x=AUE[1] , y=AUE[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUH[1] , y=AUH[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUI[1] , y=AUI[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUL[1] , y=AUL[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUNW[1], y=AUNW[0], style="i0.3c", pen="1p,black")
    fig.plot(x=AUP[1] , y=AUP[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUR[1] , y=AUR[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUS[1] , y=AUS[0] , style="i0.3c", pen="1p,black")
    fig.plot(x=AUSE[1], y=AUSE[0], style="i0.3c", pen="1p,black")
    fig.plot(x=AUW[1] , y=AUW[0] , style="i0.3c", pen="1p,black")
    fig.text(x=AUE[1] , y=AUE[0] , text=AUE[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="0.3/0.25")
    fig.text(x=AUH[1] , y=AUH[0] , text=AUH[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=AUI[1] , y=AUI[0] , text=AUI[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=AUL[1] , y=AUL[0] , text=AUL[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="0.3/0.25")
    fig.text(x=AUNW[1], y=AUNW[0], text=AUNW[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=AUP[1] , y=AUP[0] , text=AUP[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="0.73/0")
    fig.text(x=AUR[1] , y=AUR[0] , text=AUR[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="0.22/0.27")
    fig.text(x=AUS[1] , y=AUS[0] , text=AUS[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.1/-0.2")
    fig.text(x=AUSE[1], y=AUSE[0], text=AUSE[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")
    fig.text(x=AUW[1] , y=AUW[0] , text=AUW[3] , justify="RM", font="8p,Helvetica-Bold,black", offset="-0.2/0")

    # Plot scale bar
    fig.plot(x=[-153.35,-153.4378], y=[59.316,59.316], pen="2p,black")
    fig.text(x=(-153.35-153.4378)/2, y=59.32, text="5 km", justify="CB", font="8.5p,Helvetica,black")

fig.show(method="external")
fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/fig3_stations_v2.png')
fig.savefig('/Users/darrentpk/Desktop/figures/paper/newest/fig3_stations_v2.pdf')