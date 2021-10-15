#%% RELOCATE CATALOG

# This script loads in a catalog object obtained from a party object, and gives each detection a location from their parent template's event.
# After that, the script prepares a stream dictionary that corresponds to the catalog, and correlates the catalog's detections with their
# parent templates once again to create a .cc cross correlation file. The .cc file, along with an event list, station list and velocity model,
# are then parsed into GrowClust to relocate the catalog.
# All growclust outputs are then saved in the user-defined output_dir.

# Import all dependencies
import os
import shutil
import pickle
import numpy as np
import pandas as pd
import subprocess
from eqcorrscan.utils.catalog_to_dd import write_correlations, _generate_event_id_mapper
from obspy import Catalog, UTCDateTime
from obspy.core.event import Origin, Event, Comment
from toolbox import reader, writer, prepare_stream_dict, adjust_weights

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/'  # '/home/ptan/enhance_catalog/'
data_dir = '/home/ptan/enhance_catalog/data/mammoth/'
output_dir = main_dir + 'output/mammoth2/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
relocatable_catalog_filename = 'relocatable_catalog_mag.xml'
relocate_catalog_output_dir = '/Users/darrentpk/Desktop/OUTPUT/'  # '/home/ptan/enhance_catalog/output/mammoth2/relocate_catalog/'
relocated_catalog_filename = 'relocated_catalog_mag.xml'
raw_station_list_dir = main_dir + 'data/stations/'
raw_station_list_filename = 'mammoth_station_list.csv'
raw_vzmodel_dir = main_dir + 'data/vz/'
raw_vzmodel_filename = 'mammoth_vzmodel.txt'
ratio_provided = False
by_cluster = True
growclust_exe = main_dir + 'growclust/SRC/growclust'

length_actual = 8       # same as template
length_excess = 15      # in excess for stream_dict
pre_pick_actual = 1     # same as template
pre_pick_excess = 3     # in excess for stream_dict
shift_len = 1.5         # width of search for max_cc
lowcut = 1              # same as template
highcut = 10            # same as template
max_sep = 8             # max separation tolerated (8km)
min_link = 3            # minimum number of matching pick stachans (3)
min_cc = 0.75            # minimum cc to be considered pick

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

#%% Generate a stream dictionaries using prepare_stream_dict and execute cross correlations

# Load in catalog
catalog = reader(scan_data_output_dir + relocatable_catalog_filename)

# # If we are writing correlations in bulk of the whole catalog, we process the entire catalog at one go
# if not by_cluster:
#
#     # Generate stream dictionary (refer to toolbox.py)
#     # Note that we want pre_pick and length to be in excess, since write_correlations trims the data for us
#     stream_dict = prepare_stream_dict(catalog,pre_pick=pre_pick_excess,length=length_excess,local=True,data_dir=data_dir)
#
#     # Execute cross correlations and write out a .cc file using write_correlations (refer to EQcorrscan docs)
#     # Note this stores a file called "dt.cc" in your current working directory
#     event_id_mapper = write_correlations(catalog=catalog, stream_dict=stream_dict, extract_len=length_actual,
#                                          pre_pick=pre_pick_actual, shift_len=shift_len, lowcut=lowcut,
#                                          highcut=highcut, max_sep=max_sep, min_link=min_link, min_cc=min_cc,
#                                          interpolate=False, max_workers=None, parallel_process=False)
#
#     # Define source and target cc filepaths
#     original_dt_dir = os.getcwd() + '/dt.cc'
#     target_dt_dir = os.getcwd() + '/master_dt.cc'
#     adjust_weights(original_dt_dir, target_dt_dir, append=False)
#
# # If we are correlating by cluster
# else:
#
#     # Start by generating event id mapper for full catalog
#     event_id_mapper = _generate_event_id_mapper(catalog, event_id_mapper=None)
#
#     # Obtain a unique list of source templates
#     templates = []
#     for event in catalog:
#         templates.append(event.comments[0].text.split()[1])
#     unique_templates = np.unique(templates)
#
#     # Initialize a list to store the catalog index of every template's self-detection
#     template_indices = []
#
#     # Define source and target cc filepaths
#     original_dt_dir = os.getcwd() + '/dt.cc'
#     target_dt_dir = os.getcwd() + '/master_dt.cc'
#
#     # Loop through each unique template
#     for i, unique_template in enumerate(unique_templates):
#
#         # Find index of catalog events that correspond to this template
#         template_detection_indices = [i for i, template in enumerate(templates) if unique_template in template]
#
#         # Craft sub-catalog
#         sub_catalog = Catalog()
#         detect_vals = []
#         for template_detection_index in template_detection_indices:
#             template_detection = catalog[template_detection_index]
#             detect_val = float(template_detection.comments[2].text.split('=')[1])
#             detect_vals.append(detect_val)
#             sub_catalog += template_detection
#
#         # Store the index of the template's self-detection, this will be used for our final inter-cluster step
#         template_indices.append(template_detection_indices[np.argmax(detect_vals)])
#
#         # Now craft stream dictionary
#         stream_dict = prepare_stream_dict(sub_catalog, pre_pick=pre_pick_excess, length=length_excess, local=True, data_dir=data_dir)
#
#         # Execute cross correlations
#         _ = write_correlations(catalog=sub_catalog, stream_dict=stream_dict, event_id_mapper=event_id_mapper,
#                                extract_len=length_actual, pre_pick=pre_pick_actual, shift_len=shift_len,
#                                lowcut=lowcut, highcut=highcut, max_sep=max_sep, min_link=min_link,
#                                min_cc=min_cc, interpolate=False, max_workers=None, parallel_process=False)
#
#         # Write/append dt.cc to target cc file
#         if i == 0:
#             adjust_weights(original_dt_dir, target_dt_dir, append=False)
#         else:
#             adjust_weights(original_dt_dir, target_dt_dir, append=True)
#
#     # Now execute inter-template cross correlation by generating a sub-catalog containing template self-detections
#     sub_catalog = Catalog()
#     for template_index in template_indices:
#         sub_catalog += catalog[template_index]
#
#     # Now craft stream dictionary
#     stream_dict = prepare_stream_dict(sub_catalog, pre_pick=pre_pick_excess, length=length_excess, local=True,
#                                       data_dir=data_dir)
#
#     # Execute cross correlations
#     _ = write_correlations(catalog=sub_catalog, stream_dict=stream_dict, event_id_mapper=event_id_mapper,
#                            extract_len=length_actual, pre_pick=pre_pick_actual, shift_len=shift_len,
#                            lowcut=lowcut, highcut=highcut, max_sep=max_sep, min_link=min_link,
#                            min_cc=min_cc, interpolate=False, max_workers=None, parallel_process=False)
#
#     # Append dt.cc to master dt.cc
#     adjust_weights(original_dt_dir, target_dt_dir, append=True)
#
# # Write event_id_mapper
# with open(relocate_catalog_output_dir + 'event_id_mapper.pkl', 'wb') as evid_pickle:  # Pickling
#     pickle.dump(event_id_mapper, evid_pickle)
#
# print('Cross correlations done!')

# Read event_id mapper
with open(relocate_catalog_output_dir + 'event_id_mapper.pkl', 'rb') as evid_pickle:
    event_id_mapper = pickle.load(evid_pickle)

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
for dir in [growclust_in, growclust_tt, growclust_out]:
    try:
        os.mkdir(dir)
    except FileExistsError:
        print('Subdirectory %s already exists.' % dir)

# Also create all necessary output destinations for GrowClust
print('Creating necessary output destinations for growclust run...')

for output_file in ['out.growclust_boot', 'out.growclust_cat', 'out.growclust_clust', 'out.growclust_log']:
    output_fullfile = growclust_out + output_file
    try:
        open(output_fullfile, 'a').close()
    except FileExistsError:
        print('Output file %s already exists.' % output_file)

# Copy dt.cc file to its appropriate directory and rename
original_dt_dir = os.getcwd() + '/master_dt.cc'
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
    FI_comment = catalog[i].comments[-1]
    magnitude = catalog[i].magnitudes[-1]
    event = Event(origins=[Origin(time=time,latitude=relocated_event_table.lat[i],longitude=relocated_event_table.lon[i],depth=relocated_event_table.dep[i]*1000)],comments=[FI_comment],magnitudes=[magnitude])
    relocated_catalog += event
writer(relocate_catalog_output_dir + relocated_catalog_filename, relocated_catalog)


