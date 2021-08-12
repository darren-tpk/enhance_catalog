#%% RELOCATE CATALOG

# This script loads in a catalog object obtained from a party object, and gives each detection a location from their parent template's event.
# After that, the script prepares a stream dictionary that corresponds to the catalog, and correlates the catalog's detections with their
# parent templates once again to create a .cc cross correlation file. The .cc file, along with an event list, station list and velocity model,
# are then parsed into GrowClust to relocate the catalog.
# All growclust outputs are then saved in the user-defined output_dir.

# Import all dependencies
import os
import shutil
import pandas as pd
import subprocess
from eqcorrscan.utils.catalog_to_dd import write_correlations
from obspy import Catalog, UTCDateTime
from obspy.core.event import Origin, Event
from toolbox import reader, prepare_stream_dict, writer

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
data_dir = None
output_dir = main_dir + 'output/greatsitkin2/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
catalog_filename = 'party_catalog.xml'
relocate_catalog_output_dir = output_dir + 'relocate_catalog/'
relocated_catalog_filename = 'relocated_catalog.xml'
raw_station_list_dir = main_dir + 'data/avo/'
raw_station_list_filename = 'station_list.csv'
raw_vzmodel_dir = main_dir + 'data/avo/'
raw_vzmodel_filename = 'adak_vzmodel.txt'
ratio_provided = False
growclust_exe = main_dir + 'growclust/SRC/growclust'

length_actual = 8       # same as template
length_excess = 15      # in excess for stream_dict
pre_pick_actual = 1     # same as template
pre_pick_excess = 3     # in excess for stream_dict
shift_len = 1           # width of search for max_cc
lowcut = 1              # same as template
highcut = 10            # same as template
max_sep = 8             # max separation tolerated (8km)
min_link = 3            # minimum number of matching pick stachans (3)
min_cc = 0.4            # minimum cc to be considered pick

#%% Prepare GrowClust .inp file (configurations)

# define GrowClust input, travel time and output directories
growclust_in = relocate_catalog_output_dir + 'IN/'
growclust_tt = relocate_catalog_output_dir + 'TT/'
growclust_out = relocate_catalog_output_dir + 'OUT/'

# Event list format & name
# (0 = evlist, 1 = phase, 2 = GrowClust, 3 = HypoInverse)
evlist_format = '1'
evlist_filename = growclust_in + 'evlist.txt'

# Station list format & name
# (0 = SEED channel, 1 = station name)
stlist_format = '1'
stlist_filename = growclust_in + 'stlist.txt'

# Cross correlation format & name
# (0 = binary, 1 = text), (21 = tt2-tt1, 12 = tt1-tt2)
xcordata_format = '1  12'
xcordata_filename = growclust_in + 'xcordata.txt'

# Velocity model input and output name
vzmodel_filename = growclust_in + 'vzmodel.txt'
vzmodel_fileout = growclust_tt + 'vzfine.txt'

# Output travel time table P & S phase names
output_travel_time_table_P = growclust_tt + 'tt.pg'
output_travel_time_table_S = growclust_tt + 'tt.sg'

# Relocated catalog filename
outcat_filename = growclust_out + 'out.growclust_cat'

# Relocated cluster filename
outclust_filename = growclust_out + 'out.growclust_clust'

# Program log filename
outlog_filename = growclust_out + 'out.growclust_log'

# Bootstrap distribution filename
outboot_filename = growclust_out + 'out.growclust_boot'

# Other GrowClust inputs (separate by space)
# For detailed description look at GrowClust manual

# vpvs_factor, rayparam_min (-1 = default)
vpvs_rayparam = '1.732 -1'

# tt_dep0, tt_dep1, tt_ddep
tt_deps = '0. 50. 1.'

# tt_del0, tt_del1, tt_ddel
tt_dels = '0. 200. 2.'

# rmin, delmax, rmsmax
tt_thres = '0.6 80 0.2'

# rpsavgmin, rmincut, ngoodmin, iponly
tt_thres2 = '0 0 0 0'

# nboot, nbranch_min
nboot_nbranch = '10 1'


#%% Define functions

# Nil

#%% Load in tribe and detected party catalog to derive a catalog of detections WITH locations
# Note that some detections in the detected party catalog are obtained from location-less templates

# Load in tribe and catalog
tribe = reader(create_tribe_output_dir + tribe_filename)
catalog = reader(scan_data_output_dir + catalog_filename)

# Use the tribe to get a list of template names
template_names = [template.name for template in tribe]

# Initialize new catalog that stores detections that can be relocated
new_catalog = Catalog()

# Loop over the detected catalog, copying over events if their templates have locations
for event in catalog:

    # Get source event
    source_template_name = event.comments[0].text.split(' ')[1]
    source_template_index = template_names.index(source_template_name)
    source_template = tribe[source_template_index]
    source_event = source_template.event

    # Check if the source event has a location. Only continue if it does
    if source_event.origins[0].latitude is not None:

        # We extract the source event's location (lat, lon, dep)
        lat = source_event.origins[0].latitude
        lon = source_event.origins[0].longitude
        dep = source_event.origins[0].depth

        # Reformat the event's detection time to UTC
        time_text = event.resource_id.id.split('_')[-1]
        UTCdt = UTCDateTime(time_text)

        # Create an ObsPy Origin object and store it in the event
        event.origins = [Origin(time= UTCdt, latitude = lat, longitude = lon, depth = dep)]

        # Append to new catalog
        new_catalog += event

# Replace catalog in main workflow
catalog = new_catalog

#%% Give each event pick phase information based on channel convention

# Loop over events in catalog
for event in catalog:

    # Loop over picks in event
    for pick in event.picks:

        # Check for non-vertical channel, and give an 'S' hint
        if pick.waveform_id.channel_code[-1] == 'N' or pick.waveform_id.channel_code[-1] == 'E':
            pick.phase_hint = 'S'

        # Also check for a vertical channel, and give an 'P' hint
        elif pick.waveform_id.channel_code[-1] == 'Z':
            pick.phase_hint = 'P'

        # Throw out a ValueError if the channel code is not recognized
        else:
            raise ValueError

#%% Generate a stream dictionary using prepare_stream_dict

# For detailed commenting, refer to toolbox.py
# Note that we want pre_pick and length to be in excess, since write_correlations trims the data for us
stream_dict = prepare_stream_dict(catalog,pre_pick=pre_pick_excess,length=length_excess,local=False,data_dir=None)

#%% Execute cross correlations and write out a .cc file using write_correlations

# For detailed documentation, refer to EQcorrscan docs
# Note this stores a file called "dt.cc" in your current working directory
event_id_mapper = write_correlations(catalog=catalog, stream_dict=stream_dict, extract_len=length_actual,
                                     pre_pick=pre_pick_actual, shift_len=shift_len, lowcut=lowcut,
                                     highcut=highcut, max_sep=max_sep, min_link=min_link, min_cc=min_cc,
                                     interpolate=False, max_workers=None, parallel_process=True)

#%% Prepare all files for GrowClust

# Craft the .inp file
inp_file = open(relocate_catalog_output_dir + 'config.inp', 'w')
inp_file.write(evlist_format + '\n' +
               evlist_filename + '\n' +
               stlist_format + '\n' +
               stlist_filename + '\n' +
               xcordata_format + '\n' +
               xcordata_filename + '\n' +
               vzmodel_filename + '\n' +
               vzmodel_fileout + '\n' +
               output_travel_time_table_P + '\n' +
               output_travel_time_table_S + '\n' +
               vpvs_rayparam + '\n' +
               tt_deps + '\n' +
               tt_dels + '\n' +
               tt_thres + '\n' +
               tt_thres2 + '\n' +
               nboot_nbranch + '\n' +
               outcat_filename + '\n' +
               outclust_filename + '\n' +
               outlog_filename + '\n' +
               outboot_filename + '\n')
inp_file.close()

# Create all necessary sub-directories for GrowClust
print('Creating subdirectories for growclust run...')
try:
    os.mkdir(growclust_in)
    os.mkdir(growclust_tt)
    os.mkdir(growclust_out)
except FileExistsError:
    print('Subdirectories already exist')

# Copy dt.cc file to its appropriate directory and rename
original_dt_dir = os.getcwd() + '/dt.cc'
target_dt_dir = xcordata_filename
shutil.copyfile(original_dt_dir, target_dt_dir)

# Prepare station list file
raw_station_list = pd.read_csv(raw_station_list_dir + raw_station_list_filename)

# Loop through catalog event picks to get a unique list of stations used
stations_used = []
for event in catalog:
    for pick in event.picks:
        stations_used.append(pick.waveform_id.station_code)
stations_used = list(set(stations_used))

# Write GrowClust's station list input
stlist_file = open(stlist_filename, 'w')
for station_used in stations_used:
    index = list(raw_station_list.station).index(station_used)
    station_lat = raw_station_list.latitude[index]
    station_lon = raw_station_list.longitude[index]
    entry = '%s %9.4f %9.4f' % (station_used, station_lat, station_lon)
    stlist_file.write(entry + '\n')
stlist_file.close()

# Format velocity model to layer cake and write
# (Current velocity model is from Toth and Kisslinger, 1984 for Adak, and Power et al., 2012 for Redoubt)
vzmodel = pd.read_csv(raw_vzmodel_dir + raw_vzmodel_filename, header=None)
vzmodel_file = open(vzmodel_filename, 'w')
if ratio_provided:
    for i in range(len(vzmodel.values)):
        if i != len(vzmodel.values)-1:
            vz_str1 = vzmodel.values[i][0].split(' ')
            entry1 = '%5.1f %4.2f %4.2f' % (float(vz_str1[0]), float(vz_str1[1]), float(vz_str1[1])/float(vz_str1[2]))
            vz_str2 = vzmodel.values[i + 1][0].split(' ')
            entry2 = '%5.1f %4.2f %4.2f' % (float(vz_str2[0]), float(vz_str1[1]), float(vz_str1[1])/float(vz_str1[2]))
            vzmodel_file.write(entry1 + '\n' + entry2 + '\n')
        else:
            vz_str1 = vzmodel.values[i][0].split(' ')
            entry1 = '%5.1f %4.2f %4.2f' % (float(vz_str1[0]), float(vz_str1[1]), float(vz_str1[1])/float(vz_str1[2]))
            vzmodel_file.write(entry1)
else:
    for i in range(len(vzmodel.values)):
        if i != len(vzmodel.values)-1:
            vz_str1 = vzmodel.values[i][0].split(' ')
            entry1 = '%5.1f %4.2f 0.0' % (float(vz_str1[0]), float(vz_str1[1]))
            vz_str2 = vzmodel.values[i + 1][0].split(' ')
            entry2 = '%5.1f %4.2f 0.0' % (float(vz_str2[0]), float(vz_str1[1]))
            vzmodel_file.write(entry1 + '\n' + entry2 + '\n')
        else:
            vz_str1 = vzmodel.values[i][0].split(' ')
            entry1 = '%5.1f %4.2f 0.0' % (float(vz_str1[0]), float(vz_str1[1]))
            vzmodel_file.write(entry1)
vzmodel_file.close()

# Create vzfine file in TT directory
vzfine_file = open(vzmodel_fileout, 'w')
vzfine_file.close()

# Format event list from post-processed detected catalog
evlist_file = open(evlist_filename, 'w')
for event in catalog:
    time = event.origins[0].time
    yr = time.year
    mon = time.month
    day = time.day
    hr = time.hour
    min = time.minute
    sec = time.second + time.microsecond/(10**6)
    lat = event.origins[0].latitude
    lon = event.origins[0].longitude
    dep = event.origins[0].depth / 1000 # km
    mag = 0.0  # magnitude not calculated in workflow as of May 28
    eh = 0.0
    ez = 0.0
    rms = 0.0
    evid = event_id_mapper[event.resource_id.id]
    entry = '%4d %2d %2d %2d %2d %6.3f %9.5f %10.5f %7.3f %7.2f %6.3f %6.3f %6.3f %3d' % \
            (yr, mon, day, hr, min, sec, lat, lon, dep, mag, eh, ez, rms, evid)
    evlist_file.write(entry + '\n')
evlist_file.close()

#%% Run GrowClust

# Craft command and call subprocess in shell
growclust_command = growclust_exe + ' ' + relocate_catalog_output_dir + 'config.inp'
subprocess.call(growclust_command, shell=True)

#%% Get GrowClust events in catalog form and save

# Read in output catalog using pandas
relocated_event_table = pd.read_csv(outcat_filename, header=None, sep='\s+',
                                    names=['yr','mon','day','hr','mm','sec','evid','lat','lon','dep','mag','qID','cID',
                                           'nbranch','qnpair','qndiffP','qndiffS','rmsP','rmsS','eh','ez','et','latC',
                                           'lonC','depC'])

# Initialize ObsPy catalog object and populate with events on a row-by-row basis
relocated_catalog = Catalog()
for i in range(len(relocated_event_table)):
    time = UTCDateTime(relocated_event_table.yr[i],relocated_event_table.mon[i],relocated_event_table.day[i],relocated_event_table.hr[i],relocated_event_table.mm[i],relocated_event_table.sec[i])
    event = Event(origins=[Origin(time=time,latitude=relocated_event_table.lat[i],longitude=relocated_event_table.lon[i],depth=relocated_event_table.dep[i]*1000)])
    relocated_catalog += event
writer(relocate_catalog_output_dir + relocated_catalog_filename, relocated_catalog)


