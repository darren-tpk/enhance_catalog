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

####################################################################

# Load the complete event list and plot FI vs time
catalog_filepath = main_dir + 'output/redoubt2/scan_data/party_catalog_full_reFI.xml'
catalog = reader(catalog_filepath)
utctimes = np.array([UTCDateTime(event.resource_id.id.split('_')[-1]) for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)

# Full duration
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
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
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
fig.savefig('/Users/darrentpk/Desktop/event_list_reFI_time.jpg')

# Zoom to shorter period (Mar 14 to April 5)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
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
fig.savefig('/Users/darrentpk/Desktop/event_list_reFI_time_ZOOM1.jpg')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
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
fig.savefig('/Users/darrentpk/Desktop/event_list_reFI_time_ZOOM2.jpg')

####################################################################

# Load the relocated catalog and plot FI and Depth vs time
catalog_filepath = main_dir + 'output/redoubt2/relocate_catalog/relocated_catalog_reFI.xml'
catalog = reader(catalog_filepath)
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)
MAX_DEPTH = 40
valid_index = np.where(depths<MAX_DEPTH)
latitudes = latitudes[valid_index]
longitudes = longitudes[valid_index]
depths = depths[valid_index]
FI_values = FI_values[valid_index]
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
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
#ax.grid(axis='x')
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %y') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI and depth with time (Relocated Catalog, N=%d)' % len(FI_values),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='brown',alpha=0.5)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
    ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom1_base_time-base_time)/(86400*num_days), xmax=(zoom1_end_time-base_time)/(86400*num_days),color='blue', alpha=0.5)
    ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom2_base_time-base_time)/(86400*num_days), xmax=(zoom2_end_time-base_time)/(86400*num_days), color='blue', alpha=0.5)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/relocated_catalog_reFIdepth_time.jpg')

# Zoom to shorter period (Mar 14 to April 5)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI and depth with time (Relocated Catalog, N=%d)' % len(FI_values),fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Relocated Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/relocated_catalog_reFIdepth_time_ZOOM1.jpg')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30)
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI and depth with time (Relocated Catalog, N=%d)' % len(FI_values),fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Relocated Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
ax.text(1,-0.3,'Event pair correlation must still be cross-correlating this huge cluster(!)', fontsize=15)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/relocated_catalog_reFIdepth_time_ZOOM2.jpg')


####################################################################

# Load the original catalog and plot FI vs time
catalog_filepath = main_dir + 'output/redoubt2/scan_data/PEC_events_reFI.xml'
catalog = reader(catalog_filepath)
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)
MAX_DEPTH = 40
valid_index = np.where(depths<MAX_DEPTH)
latitudes = latitudes[valid_index]
longitudes = longitudes[valid_index]
depths = depths[valid_index]
FI_values = FI_values[valid_index]

# Full duration
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
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %y') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI with time (Original Catalog, N=%d)' % len(FI_values),fontsize=15)
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
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFI_time.jpg')

# Zoom to shorter period (Mar 14 to April 5)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Original Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFI_time_ZOOM1.jpg')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,FI_values,s=10,c='red',alpha=0.5)
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Frequency Index',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Original Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFI_time_ZOOM2.jpg')

####################################################################

# Load the original catalog and plot FI and Depth vs time

# Full duration
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
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
#ax.grid(axis='x')
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %y') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('Variation of earthquake FI and depth with time (Original Catalog, N=%d)' % len(FI_values),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='brown',alpha=0.5)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
    ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom1_base_time-base_time)/(86400*num_days), xmax=(zoom1_end_time-base_time)/(86400*num_days),color='blue', alpha=0.5)
    ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom2_base_time-base_time)/(86400*num_days), xmax=(zoom2_end_time-base_time)/(86400*num_days), color='blue', alpha=0.5)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFIdepth_time.jpg')

# Zoom to shorter period (Mar 14 to April 5)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Original Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFIdepth_time_ZOOM1.jpg')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(15, 4))
main = ax.scatter(days,depths,c=FI_values,edgecolors='k',cmap='seismic',alpha=0.5)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=13)
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
ax.set_ylabel('Depth',fontsize=15)
ax.set_xlabel('UTC Time',fontsize=15)
included_utctimes = utctimes[(utctimes > base_time) & (utctimes < end_time)]
ax.set_title('Variation of earthquake FI with time (Original Catalog, N=%d)' % len(included_utctimes),fontsize=15)
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.2)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/original_catalog_reFIdepth_time_ZOOM2.jpg')

#####################

# combined plot

catalog_filepath = main_dir + 'output/redoubt2/scan_data/party_catalog_full_reFI.xml'
catalog = reader(catalog_filepath)
relocated_catalog_filepath = main_dir + 'output/redoubt2/relocate_catalog/relocated_catalog_reFImag.xml'
relocated_catalog = reader(relocated_catalog_filepath)
original_filepath = main_dir + 'output/redoubt2/scan_data/PEC_events_reFI.xml'
original = reader(original_filepath)
size_by_mag = True

# event list
utctimes = np.array([UTCDateTime(event.resource_id.id.split('_')[-1]) for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)

# relocated catalog
relocated_utctimes = np.array([event.origins[0].time for event in relocated_catalog])
relocated_FI_values = []
for event in relocated_catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        relocated_FI_values.append(np.nan)
    else:
        relocated_FI_values.append(float(event.comments[-1].text.split('=')[1]))
relocated_FI_values = np.array(relocated_FI_values)

# remove overlapping
remove_index = []
for i, relocated_utctime in enumerate(relocated_utctimes):
    remove_index.append(np.argmin(np.abs(utctimes-relocated_utctime)))
    print(i)
sorted_remove_index = list(np.sort(remove_index))
utctimes = list(utctimes)
FI_values = list(FI_values)
for i in reversed(sorted_remove_index):
    utctimes.pop(i)
    FI_values.pop(i)
utctimes = np.array(utctimes)
FI_values = np.array(FI_values)

original_utctimes = np.array([event.origins[0].time for event in original])
original_FI_values = []
for event in original:
    if event.comments[-1].text.split('=')[1] == 'None':
        original_FI_values.append(np.nan)
    else:
        original_FI_values.append(float(event.comments[-1].text.split('=')[1]))
original_FI_values = np.array(original_FI_values)

def calibrate_sizes(magnitudes):
    # construct reference sizes
    ref_magnitudes = np.array(range(0, 401))
    min_mag_marker = 1
    ref_sizes = []
    for magnitude in list(ref_magnitudes):
        normalized_magnitude = (magnitude - np.min(min_mag_marker)) / np.ptp(
            ref_magnitudes[ref_magnitudes > min_mag_marker])
        ref_size = 50 + normalized_magnitude * (5*160)
        ref_sizes.append(ref_size)
    ref_magnitudes = list(ref_magnitudes)
    # now calculate sizes
    sizes = []
    for magnitude in magnitudes:
        if magnitude < min_mag_marker:
            size = ref_sizes[0]
        else:
            size = ref_sizes[ref_magnitudes.index(int(magnitude * 100))]
        sizes.append(size)
    sizes = np.array(sizes)
    return sizes

if size_by_mag:
    relocated_magnitudes = np.array([event.magnitudes[0].mag for event in relocated_catalog])
    relocated_sizes = calibrate_sizes(relocated_magnitudes)
    original_magnitudes = np.array([event.magnitudes[0].mag for event in original])
    original_sizes = calibrate_sizes(original_magnitudes)

# Full duration
base_time = UTCDateTime(2008,5,1)
end_time = UTCDateTime(2009,9,1)
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
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
fig, ax = plt.subplots(figsize=(44, 6.5))
if size_by_mag:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days,relocated_FI_values,s=relocated_sizes,c='indianred',alpha=0.4)
    main = ax.scatter(original_days,original_FI_values,s=original_sizes,c='darkred',alpha=0.7)
else:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days, relocated_FI_values, s=50, c='indianred', alpha=0.5)
    main = ax.scatter(original_days, original_FI_values, s=50, c='darkred')
#### PLOT LEGEND
ax.scatter(8, -1.7, s=500, marker='s', c='#505050', alpha=0.5)
ax.scatter(8, -1.5, s=500, marker='s', c='indianred', alpha=0.5)
ax.scatter(8, -1.3, s=500, marker='s', c='darkred')
ax.scatter(8, -1.7, s=500, marker='s', c='None', edgecolor='k')
ax.scatter(8, -1.5, s=500, marker='s', c='None', edgecolor='k')
ax.scatter(8, -1.3, s=500, marker='s', c='None', edgecolor='k')
ax.text(10.7, -1.7, ': Enhanced Catalog w/o Magnitudes', fontsize=24, va='center')
ax.text(10.7, -1.5, ': Enhanced Catalog with Magnitudes', fontsize=24, va='center')
ax.text(10.7, -1.3, ': Original AVO Catalog', fontsize=24, va='center')
ax.scatter(8, -1.0, s=50, marker='o', c='None', edgecolor='k')
ax.scatter(17, -1.0, s=250, marker='o', c='None', edgecolor='k')
ax.scatter(26, -1.0, s=450, marker='o', c='None', edgecolor='k')
ax.scatter(35, -1.0, s=650, marker='o', c='None', edgecolor='k')
ax.text(8, -0.8, 'M0', fontsize=24, va='center', ha='center')
ax.text(17, -0.8, 'M1', fontsize=24, va='center', ha='center')
ax.text(26, -0.8, 'M2', fontsize=24, va='center', ha='center')
ax.text(35, -0.8, 'M3', fontsize=24, va='center', ha='center')
####
ax.set_xlim([0,num_days])
ax.set_ylim([-2,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b \'%y') for tick_time in tick_times],fontsize=26,rotation=30)
ax.tick_params(axis='y',labelsize=26)
ax.tick_params(which='major',length=10,width=2,direction='in')
ax.set_ylabel('Frequency Index',fontsize=26)
# ax.set_xlabel('UTC Time',fontsize=28)
ax.set_title('Variation of earthquake FI with time (Temporally Complete Event List, N=%d)' % len(catalog),fontsize=30,fontweight='bold')
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.4)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
    ax.axhspan(ymin=-2, ymax=-1.95, xmin=(zoom1_base_time-base_time)/(86400*num_days), xmax=(zoom1_end_time-base_time)/(86400*num_days),color='blue', alpha=0.5)
    ax.axhspan(ymin=-2, ymax=-1.95, xmin=(zoom2_base_time-base_time)/(86400*num_days), xmax=(zoom2_end_time-base_time)/(86400*num_days), color='blue', alpha=0.5)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/combined_poster_FINAL.jpg',bbox_inches='tight')

# Zoom to shorter period (Mar 14 to April 5)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(30.3, 6.5))
if size_by_mag:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days,relocated_FI_values,s=relocated_sizes,c='indianred',alpha=0.4)
    main = ax.scatter(original_days,original_FI_values,s=original_sizes,c='darkred',alpha=0.7)
else:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days, relocated_FI_values, s=50, c='indianred', alpha=0.5)
    main = ax.scatter(original_days, original_FI_values, s=50, c='darkred')
ax.set_xlim([0,num_days])
ax.set_ylim([-1.999,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times],fontsize=26,rotation=30,ha='right')
ax.tick_params(axis='y',labelsize=26)
ax.tick_params(which='major',length=10,width=2,direction='in')
ax.set_ylabel('Frequency Index',fontsize=26)
# ax.set_xlabel('UTC Time',fontsize=15)
num_events = len(utctimes[(utctimes > base_time) & (utctimes < end_time)]) + len(relocated_utctimes[(relocated_utctimes > base_time) & (relocated_utctimes < end_time)])
ax.set_title('(N=%d)' % num_events,fontsize=30,fontweight='bold')
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.4)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/ZOOM1_poster_FINAL.jpg',bbox_inches='tight')

# Zoom to second short period (May 2 to May 9)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
fig, ax = plt.subplots(figsize=(9.65, 6.5))
if size_by_mag:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days,relocated_FI_values,s=relocated_sizes,c='indianred',alpha=0.4)
    main = ax.scatter(original_days,original_FI_values,s=original_sizes,c='darkred',alpha=0.7)
else:
    main = ax.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
    main = ax.scatter(relocated_days, relocated_FI_values, s=50, c='indianred', alpha=0.5)
    main = ax.scatter(original_days, original_FI_values, s=50, c='darkred')
ax.set_xlim([0,num_days])
ax.set_ylim([-1.999,0.8])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times],fontsize=26,rotation=30,ha='right')
ax.tick_params(axis='y',labelsize=26)
ax.tick_params(which='major',length=10,width=2,direction='in')
# ax.set_ylabel('Frequency Index',fontsize=26)
# ax.set_xlabel('UTC Time',fontsize=15)
num_events = len(utctimes[(utctimes > base_time) & (utctimes < end_time)]) + len(relocated_utctimes[(relocated_utctimes > base_time) & (relocated_utctimes < end_time)])
ax.set_title('(N=%d)' % num_events,fontsize=30,fontweight='bold')
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='k',alpha=0.4)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=0.6,ymax=0.79,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/ZOOM2_poster_FINAL.jpg',bbox_inches='tight')

# Plot truncated

# Load the relocated catalog and plot FI and Depth vs time
catalog_filepath = main_dir + 'output/redoubt2/relocate_catalog/relocated_catalog_reFImag.xml'
catalog = reader(catalog_filepath)
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)
MAX_DEPTH = 40
valid_index = np.where(depths<MAX_DEPTH)
latitudes = latitudes[valid_index]
longitudes = longitudes[valid_index]
depths = depths[valid_index]
FI_values = FI_values[valid_index]
base_time = UTCDateTime(2009,1,1)
end_time = UTCDateTime(2009,7,1)
days = (utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [UTCDateTime(2009, 1, 1), UTCDateTime(2009, 2, 1),
              UTCDateTime(2009, 3, 1), UTCDateTime(2009, 4, 1),
              UTCDateTime(2009, 5, 1), UTCDateTime(2009, 6, 1),
              UTCDateTime(2009, 7, 1)]
relocated_sizes = calibrate_sizes(relocated_magnitudes)
fig, ax = plt.subplots(figsize=(22, 10))
main = ax.scatter(days,depths,c=FI_values,s=relocated_sizes,edgecolors='k',cmap='seismic',alpha=0.6)
main.set_clim([-2, 0.5])
cb = fig.colorbar(main)
cb.set_label('Frequency Index',fontsize=28)
cb.ax.tick_params(labelsize=28)
##### CREATE MAG LEGEND
ax.scatter(8, 12.6, s=50, marker='o', c='None', edgecolor='k')
ax.scatter(17, 12.6, s=250, marker='o', c='None', edgecolor='k')
ax.scatter(26, 12.6, s=450, marker='o', c='None', edgecolor='k')
ax.scatter(35, 12.6, s=650, marker='o', c='None', edgecolor='k')
ax.text(8, 11.6, 'M0', fontsize=28, va='center', ha='center')
ax.text(17, 11.6, 'M1', fontsize=28, va='center', ha='center')
ax.text(26, 11.6, 'M2', fontsize=28, va='center', ha='center')
ax.text(35, 11.6, 'M3', fontsize=28, va='center', ha='center')
######
#ax.grid(axis='x')
ax.set_xlim([0,num_days])
ax.set_ylim([-5,14])
ax.invert_yaxis()
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%b %y') for tick_time in tick_times],fontsize=28,rotation=30,ha='right')
ax.tick_params(axis='y',labelsize=28)
ax.tick_params(which='major',length=10,width=2,direction='in')
ax.set_ylabel('Depth',fontsize=28)
# ax.set_xlabel('UTC Time',fontsize=15)
ax.set_title('FI and Depth vs Time (Relocated Catalog, N=%d)' % len(FI_values),fontsize=30,fontweight='bold')
if add_labels:
    explosion_days = [(explosion_time-base_time)/86400 for explosion_time in explosion_times]
    for explosion_day in explosion_days:
        ax.axvline(x=explosion_day,color='black',alpha=0.4,zorder=0)
    color_code_days = [(color_code_time-base_time)/86400 for color_code_time in color_code_times]
    for k in range(len(color_code_days)-1):
        ax.axhspan(ymin=-4.95,ymax=-4.2,xmin=color_code_days[k]/num_days,xmax=color_code_days[k+1]/num_days,color=color_codes[k],zorder=10)
    # ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom1_base_time-base_time)/(86400*num_days), xmax=(zoom1_end_time-base_time)/(86400*num_days),color='blue', alpha=0.5)
    # ax.axhspan(ymin=13, ymax=13.9, xmin=(zoom2_base_time-base_time)/(86400*num_days), xmax=(zoom2_end_time-base_time)/(86400*num_days), color='blue', alpha=0.5)
fig.show()
fig.savefig('/Users/darrentpk/Desktop/mag_plots/FI_depth_mag_time_FINAL.jpg',bbox_inches='tight')
