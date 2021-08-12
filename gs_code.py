#%% CONVERT REDPY

# This script aims to convert the output from the REDPy tool into ObsPy catalog objects.
# In the process, the REDPy output is compared to a pre-existing catalog to generate a list of associated cores.
# Lastly, the code can also revisit the seismic data to use the STA/LTA coincidence trigger as P pick times.
# All .xml catalog and .txt files created are saved in the user-defined output_dir.
# Note that [PEC = pre-existing catalog]

# Import all dependencies
import os
import pickle
import time
import pandas as pd
import numpy as np
from obspy import Catalog, UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Comment, Pick, WaveformStreamID, Arrival, ResourceIdentifier
from obspy.signal.trigger import coincidence_trigger, classic_sta_lta, trigger_onset
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from toolbox import read_trace, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = None #'/home/ptan/enhance_catalog/data/mammoth/'  # redoubt data directory on local
output_dir = main_dir + 'output/greatsitkin2/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
redpy_results_dir = main_dir + 'redpy_results/greatsitkin_0731_0810/'
PEC_dir = main_dir + 'data/avo/'
hypoi_file = 'greatsitkin_20210731_20210810_hypoi.txt'
hypoddpha_file = 'greatsitkin_20210731_20210810_hypoddpha.txt'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks
redpy_network_list = ['AV','AV','AV','AV','AV']
redpy_station_list = ['GSCK','GSIG','GSMY','GSSP','GSTR']
redpy_channel_list = [['BHZ'],['BHZ'],['BHZ'],['BHZ'],['BHZ']]
redpy_location_list = ['--','--','--','--','--']
tolerance = 4e4  # tolerance for boxcar removal from data (as a factor to median)
local = False
client_name = 'IRIS' # IRIS/NCEDC

#%% Define functions

# Function to pull a core catalog out of a full catalog
def pull_cores(full_catalog):
    print('\nExtracting core events...')
    time_start = time.time()

    # Initialize core catalog object
    all_cores = Catalog()

    # Loop through full catalog to check for cores
    for detection in full_catalog:
        if detection.origins[0].comments[0].text.split(' ')[2] == 'CORE':
            all_cores.append(detection)

    # Find repeats and save their indices
    remove_index = []
    for ii in range(len(all_cores)):
        if all_cores[ii].origins[0].time == all_cores[ii - 1].origins[0].time:
            remove_index.append(ii)

    # Reappend detections to core catalog, avoiding repeat indices
    core_catalog = Catalog()
    for ii, core in enumerate(all_cores):
        if ii not in remove_index:
            core_catalog.append(core)

    # Conclude process
    time_stop = time.time()
    print('Core extraction complete, processing time: %.2f s' % (time_stop - time_start))
    return core_catalog

#%% Prepare REDPy output for catalog creation

# Read redpy text files using pandas
redpy_detections = pd.read_csv(redpy_results_dir + 'catalog.txt', sep=' ', names=['Cluster', 'DateTime'])
redpy_cores = pd.read_csv(redpy_results_dir + 'cores.txt', sep=' ', names=['Cluster', 'DateTime'])
redpy_orphans = pd.read_csv(redpy_results_dir + 'orphancatalog.txt', names=['DateTime'])

# Prepare core and orphan datetimes as an array
core_event_times = np.array([UTCDateTime(t) for t in redpy_cores.DateTime])
orphan_event_times = np.array([UTCDateTime(t) for t in redpy_orphans.DateTime])

# Prepare AVO catalog with known picks
hypoi_path = PEC_dir + hypoi_file
hypoddpha_path = PEC_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path, channel_convention=True)
PEC_events = read_hypoddpha(hypoi_path, hypoddpha_path, channel_convention=True)
PEC_event_times = np.array([PEC_event.origins[0].time for PEC_event in PEC_events])

#%% Create ObsPy catalog for REDPy results

# Prepare output directory
print('Creating subdirectories for workflow outputs...')
output_subdirs = [output_dir, output_dir + 'convert_redpy/', output_dir + 'create_tribe/',
                  output_dir + 'scan_data/', output_dir + 'relocate_catalog/']
for output_subdir in output_subdirs:
    try:
        os.mkdir(output_subdir)
    except FileExistsError:
        print('%s already exists. Resuming...' % output_subdir)
print('All output subdirectories created.')

# Initialize catalog object and lists before commencing loop
redpy_catalog = Catalog()  # empty catalog to populate\
associated_cluster_list = []  # store unique cluster numbers that have associated catalog events
unmatched_indices = list(range(len(PEC_events)))  # list of avo event indices. matched indices are removed

# Loop through redpy detection list
print('\nCreating ObsPy catalog for REDPy results...')
time_start = time.time()
for i in range(len(redpy_detections)):

    # Extract values
    cluster = redpy_detections.Cluster[i]
    detection_time = UTCDateTime(redpy_detections.DateTime[i])

    # Check if event is a core
    if min(abs(core_event_times - detection_time)) == 0:
        event_tag = 'CLUSTER ' + str(cluster) + ' CORE EVENT;'
        redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text=event_tag)])])

    # Otherwise it is a cluster event
    else:
        event_tag = 'CLUSTER ' + str(cluster) + ' EVENT;'
        redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text=event_tag)])])

    # Check if event is part of pre-existing catalog, using a inter-event tolerance defined by user
    if min(abs(PEC_event_times - detection_time)) < max_dt:

        # Find the closest event in time from the pre-existing catalog
        PEC_index = np.argmin(abs(PEC_event_times - detection_time))
        PEC_event = PEC_events[PEC_index]

        # Use the closest event's information to fill up event object
        redpy_event.picks = PEC_event.picks
        redpy_event.magnitudes = PEC_event.magnitudes
        redpy_event.origins[0].longitude = PEC_event.origins[0].longitude
        redpy_event.origins[0].latitude = PEC_event.origins[0].latitude
        redpy_event.origins[0].depth = PEC_event.origins[0].depth
        redpy_event.origins[0].arrivals = PEC_event.origins[0].arrivals
        redpy_event.origins[0].comments[0].text += ' IN PEC, DT = ' + '%.2f' % min(abs(PEC_event_times - detection_time))

        # Add cluster number to valid clusters list
        associated_cluster_list.append(cluster)

        # Remove avo index from unmatched list if it still exists
        if PEC_index in unmatched_indices:
            unmatched_indices.remove(PEC_index)

    # If event is not part of the PEC, we tag it as so
    else:
        redpy_event.origins[0].comments[0].text += ' NOT IN PEC'

    # Finally we append the event
    redpy_catalog.append(redpy_event)

# Write the redpy catalog to an xml file
writer(convert_redpy_output_dir + 'redpy_catalog.xml', redpy_catalog)

# Get unique list of AVO-associated clusters and non-associated clusters
associated_clusters = list(np.unique(np.array(associated_cluster_list)))
unassociated_clusters = [cluster for cluster in list(np.unique(np.array(redpy_detections.Cluster)))
                         if cluster not in associated_clusters]

# Write them to text file pickles
with open(convert_redpy_output_dir + 'unassociated_clusters.txt', 'wb') as cluster_pickle:  # Pickling
    pickle.dump(unassociated_clusters, cluster_pickle)
    pickle.dump(associated_clusters, cluster_pickle)
with open(convert_redpy_output_dir + 'unmatched_indices.txt', 'wb') as unmatched_pickle:  # Pickling
    pickle.dump(unmatched_indices, unmatched_pickle)

# Also generate unmatched PEC catalog (that can be used as templates later)
unmatched_PEC_events = Catalog()
for j, PEC_event in enumerate(PEC_events):
    if j in unmatched_indices:
        unmatched_PEC_events += PEC_event

# Write out unmatched PEC catalog to .xml file
writer(convert_redpy_output_dir + 'unmatched_PEC_events.xml', unmatched_PEC_events)

# Conclude process
time_stop = time.time()
print('Catalog object created, processing time: %.2f s' % (time_stop-time_start))

#%% Pull core catalog from full catalog (picks not added)

# Call pull_cores function
core_catalog = pull_cores(redpy_catalog)

# Write core catalog to .xml file
writer(convert_redpy_output_dir + 'core_catalog.xml', core_catalog)

#%% Generate orphan catalog (picks not added)

# Commence process
print('\nNow processing orphans...')
time_start = time.time()

# Initialize empty catalog
orphan_catalog = Catalog()

# Loop through redpy orphan list and populate orphan catalog
for orphan_time in orphan_event_times:

    # Create orphan event
    event_tag = 'ORPHAN EVENT;'
    orphan_event = Event(origins=[Origin(time=orphan_time, comments=[Comment(text=event_tag)])])

    # Check if event is part of PEC, using a inter-event tolerance defined by user
    if min(abs(PEC_event_times - orphan_time)) < max_dt:

        # Find the closest event in time within PEC
        PEC_index = np.argmin(abs(PEC_event_times - orphan_time))
        PEC_event = PEC_events[PEC_index]

        # Use closest event information to fill up event object
        orphan_event.picks = PEC_event.picks
        orphan_event.magnitudes = PEC_event.magnitudes
        orphan_event.origins[0].longitude = PEC_event.origins[0].longitude
        orphan_event.origins[0].latitude = PEC_event.origins[0].latitude
        orphan_event.origins[0].depth = PEC_event.origins[0].depth
        orphan_event.origins[0].arrivals = PEC_event.origins[0].arrivals
        orphan_event.origins[0].comments[0].text += ' IN PEC, DT = ' + '%.2f' % min(abs(PEC_event_times - orphan_time))

    # If event is not part of the PEC, we tag it otherwise
    else:
        orphan_event.origins[0].comments[0].text += ' NOT IN AVO'

    # Finally we append the event
    orphan_catalog.append(orphan_event)

# Conclude process
time_stop = time.time()
print('Orphan catalog created, processing time: %.2f s' % (time_stop-time_start))

#%% Revisit redpy catalog to add picks using single channel STA/LTA

# Commence process
print('\nMaking picks for non-associated cluster events...')
time_start = time.time()

# Define client if using data from server
if not local:
    client = Client(client_name)

# Define some coincidence_trigger arguments

# Loop through unassociated clusters
for unassociated_cluster in unassociated_clusters:

    # Find all detections within unassociated cluster
    contribution_list = list(np.where(np.array(redpy_detections.Cluster) == unassociated_cluster)[0])

    # Loop through these unassociated detections to add picks
    for contribution_index in contribution_list:

        # Determine time difference to offset pick times
        contribution_time = redpy_catalog[contribution_index].origins[0].time

        # Retrieve time limits for data fetch (+/- 12s window)
        starttime = contribution_time - 12
        endtime = contribution_time + 12

        # Remove some stations if necessary
        redpy_stations = redpy_station_list
        # # REMOVE RSO AFTER FIRST EXPLOSION [HARD CODED]
        # if UTCDateTime(2009, 3, 24, 0, 0, 0) < starttime < UTCDateTime(2009, 4, 16, 0, 0, 0):
        #     redpy_stations = ['RDN', 'REF']

        # Gather data from local machine in a +/- 12s window, filter, and taper
        stream = Stream()
        for k, redpy_station in enumerate(redpy_stations):
            for redpy_channel in redpy_channel_list[k]:
                if local:
                    station_tr = read_trace(data_dir=data_dir, station=redpy_station, channel=redpy_channel,
                                            starttime=starttime, endtime=endtime, tolerance=tolerance)
                else:
                    station_tr = client.get_waveforms(redpy_network_list[k], redpy_station, redpy_location_list[k], redpy_channel, starttime, endtime)
                stream = stream + station_tr
        stream = stream.split()
        stream = stream.filter('bandpass',freqmin=1.0, freqmax=10.0, corners=2, zerophase=True)
        stream = stream.taper(0.05, type='hann', max_length=(0.75*1024/100))  # [HARD CODED]
        stream = stream.merge(method=1, fill_value=0)
        stream = stream.trim(starttime=starttime,endtime=endtime)

        # Use coincidence trigger to get a pick time estimate
        try:
            coin_trigger = []
            coin_trigger = coincidence_trigger('classicstalta', 3, 2, stream, 2, sta=0.7, lta=8, details=True)
        # If it fails (usually due to data being too gappy), move to next event
        except:
            continue

        # If there are no coincidence triggers, move to next event
        if not coin_trigger:
            continue

        # Otherwise, continue to single channel STA/LTA pick making
        else:

            # Extract coincidence trigger time
            coin_trigger_time = coin_trigger[0]["time"]

            # For each channel,
            for tr in stream:

                # Calculate the value of the characteristic function
                sampling_rate = tr.stats.sampling_rate
                cft = classic_sta_lta(tr.data, int(0.7 * sampling_rate), int(8 * sampling_rate))
                # To plot function, use: plot_trigger(tr, cft, 3, 2)

                # Obtain trigger limits
                trigger_limits = np.array(trigger_onset(cft, 3, 2))

                # If there exists some trigger limits
                if trigger_limits.size != 0:

                    # Convert to UTCDateTime and find the trigger-on time closest to the coincidence trigger
                    trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
                    pick_time = trigger_on[np.argmin(abs(trigger_on-coin_trigger_time))]

                    # Craft waveform stream ID
                    tr_id = tr.id.split('.')
                    waveform_id = WaveformStreamID(tr_id[0],tr_id[1],tr_id[2],tr_id[3])

                    # Create ObsPy pick and ObsPy arrival objects and add to redpy_catalog
                    add_pick = Pick(time=pick_time,waveform_id=waveform_id,phase_hint='P')
                    add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id),phase='P',time_weight=0.1)
                    redpy_catalog[contribution_index].picks.append(add_pick)
                    redpy_catalog[contribution_index].origins[0].arrivals.append(add_arrival)

# Write the fully picked redpy catalog to an xml file
writer(convert_redpy_output_dir + 'redpy_catalog_picked.xml', redpy_catalog)

# Conclude process
time_stop = time.time()
print('Pick making complete, processing time: %.2f s' % (time_stop-time_start))

#%% Pull the fully picked core catalog from picked redpy catalog

# Call pull_cores function
core_catalog = pull_cores(redpy_catalog)

# Write picked core catalog to .xml file
writer(convert_redpy_output_dir + 'core_catalog_picked.xml', core_catalog)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

# Import all dependencies
import time
import obspy
from obspy import UTCDateTime, Catalog, Stream
from obspy.clients.fdsn import Client
from eqcorrscan.core.match_filter.tribe import Tribe
from toolbox import get_local_stations, prepare_catalog_stream, reader, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = None
output_dir = main_dir + 'output/greatsitkin2/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
sitelist_dir = main_dir + 'data/avo/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
channel_convention = True  # strict compliance for P/S picks on vertical/horizontal components
resampling_frequency = 50  # for resampling traces prior to final merge
tolerance = 4e4            # tolerance for boxcar removal from data (as a factor to median)
lowcut = 1.0               # filter lowcut (Hz), follows Jeremy's recommendation
highcut = 10.0             # filter highcut (Hz), follows Jeremy's recommendation
samp_rate = 50.0           # new sampling rate (Hz)
length = 8.0               # template length (s), Wech et al. (2018) chose 30s
filt_order = 4             # number of corners for filter
prepick = 1.0              # pre-pick time (s), Wech et al. (2018) chose 5s
process_len = 86400        # length to process data in (s)
min_snr = 2                # minimum SNR, Jeremy's recommendation was 5.0 (same as EQcorrscan tutorial)
local_volcano = 'great sitkin'  # for get_local_stations function, since we only take picks from stations near Redoubt
local_radius = 25          # for get_local_stations function; radius around volcano to accept stations
local = False               # if set to True, use data from local machine
client_name = 'IRIS'        # client name for back-up or non-local data query

#%% Define functions

# Nil

#%% Prepare desired catalog to undergo tribe creation

# Read in, and combine, the picked core catalog and unmatched PEC catalog
core_catalog_picked = reader(convert_redpy_output_dir + 'core_catalog_picked.xml')
unmatched_PEC_events = reader(convert_redpy_output_dir + 'unmatched_PEC_events.xml')
catalog = core_catalog_picked + unmatched_PEC_events

# Clean catalog to only include picks from local stations
local_stations = get_local_stations(sitelist_dir,local_volcano,local_radius)
catalog_out = Catalog()
for event in catalog:
    catalog_out.append(event.copy())
    catalog_out[-1].picks = []
    for pick in event.picks:
        if pick.waveform_id.station_code in local_stations:
            catalog_out[-1].picks.append(pick)
catalog = catalog_out

#%% Loop over catalog events to create templates, populating tribe

# Initialize tribe and tracker for valid events
tribe = Tribe()
valid_event = []
time_start = time.time()

# Loop through events in catalog
print('\nCommencing tribe creation. Reminder: index is -1 from event count.')
for k in range(len(catalog)):

    # Isolate event in loop
    print('\nNow at event %d out of %d' % (k+1, len(catalog)))
    event = catalog[k:k+1]

    # Prepare catalog stream, trim to the start and end of the day to enforce daylong processing
    if local:
        stream = prepare_catalog_stream(data_dir,event,resampling_frequency,tolerance)
        day_start = UTCDateTime(event[0].origins[0].time.date)
        day_end = day_start + 86400
        stream = stream.trim(starttime=day_start, endtime=day_end)

        # If the stream has traces, check the length of each trace and remove those that are too short for processing
        # (EQcorrscan produces an error if data is <80% of process_len)
        if stream is not None:
            for trace in stream:
                trace_length = trace.stats.npts / trace.stats.sampling_rate
                if trace_length < (0.8 * process_len):
                    print('%s.%s got removed due to insufficient length.' % (trace.stats.station,trace.stats.channel))
                    stream.remove(trace)

    # Start with an empty template
    template = None

    # Try to construct tribe using local data file
    if local:
        try:
            template = Tribe().construct(
                method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
                filt_order=filt_order, prepick=prepick, meta_file=event, st=stream, process=True,
                process_len=process_len, min_snr=min_snr, parallel=True)
        except:
            print('WARNING: local data failed to produce template, using client method instead.')

            # If local data files fail (e.g. too gappy), we try to construct the tribe using IRIS client downloads
            try:
                client = Client(client_name)
                template = Tribe().construct(
                    method="from_client", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
                    filt_order=filt_order, prepick=prepick, client_id=client, catalog=event, process=True,
                    process_len=process_len, min_snr=min_snr, parallel=True)
            except:
                print('WARNING: data not available on client either, skipping.')

    # If we are using downloaded data, go straight into querying the client
    else:
        try:
            client = Client(client_name)
            template = Tribe().construct(
                method="from_client", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
                filt_order=filt_order, prepick=prepick, client_id=client, catalog=event, process=True,
                process_len=process_len, min_snr=min_snr, parallel=True)
        except:
            print('WARNING: data not available on client, skipping.')


    # Check if template creation was unsuccessful:
    if template is None or len(template) == 0:
        print('Event %d failed.' % (k + 1))
        valid_event.append(0)

    # Append template to tribe if successful
    else:
        tribe = tribe + template
        time_current = time.time()
        print('Event %d passed.' % (k+1))
        valid_event.append(1)

# Conclude process
time_end = time.time()
print('\nTemplate creation complete. Time taken: %.2f s' % (time_end-time_start))
print('%d out of %d events in the catalog were converted to templates.' % (sum(valid_event),len(catalog)))

# Write tribe out into output directory
writer(create_tribe_output_dir + tribe_filename, tribe)

#%% Option to execute bulk EQcorrscan instead (problematic as of Mar 29)

# # Prepare stream object that spans catalog (4 months)
# stream = prepare_catalog_stream(data_dir,catalog,resampling_frequency,tolerance)
# # Now build templates
# print('Constructing templates ...')
# run_time = UTCDateTime()
# tribe = Tribe().construct(
#     method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
#     filt_order=filt_order, prepick=prepick, meta_file=catalog, st=stream, process=True,
#     process_len=process_len, min_snr=min_snr, parallel=True)
# print('Time Elapsed: ',UTCDateTime()-run_time)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

# Import all dependencies
import pickle
import time
import glob
import numpy as np
from eqcorrscan import Party, Tribe
from obspy import UTCDateTime, Stream, Catalog, read
from obspy.clients.fdsn import Client
from toolbox import remove_boxcars, reader, writer

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
data_dir = None
output_dir = main_dir + 'output/greatsitkin2/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
party_filename = 'party.tgz'
catalog_filename = 'party_catalog.xml'
repicked_catalog_filename = None
min_stations = 3                               # to remove templates that are anchored by too little stations
min_picks = 0                                  # to remove templates that are anchored by too little picks
start_time = UTCDateTime(2021, 7, 31, 0, 0, 0)  # start: UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2021, 8, 10, 0, 0, 0)    # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50                                 # to resample streams to match tribes
tolerance = 4e4                                # for boxcar removal
threshold_type = 'av_chan_corr'                # also consider 'MAD', Jeremy uses 12
threshold = 0.60                               # used to be 0.74 from sensitivity test. Was too high
trig_int = 8                                   # trigger interval for each template. Also used to remove repeats (decluster)
parallel_process = 'True'                      # parallel process for speed
generate_repicked_catalog = False              # option to use lag_calc to do catalog repicking
local = False                                  # if set to True, use data from local machine
client_name = 'IRIS'                           # client name for non-local data query

#%% Define functions

# Nil

#%% Read in tribe and clean using min_stations and min_picks

# Read in tribe
tribe = reader(create_tribe_output_dir + tribe_filename)

# Clean tribe using min_stations and min_picks
print('\nBefore cleaning:')
print(tribe)

# Remove based on min_stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= min_stations]
print('\nAfter removing templates with < %d stations...' % min_stations)
print(tribe)

# Remove based on min_picks
new_tribe = Tribe()
for template in tribe:
    if len(template.st) >= min_picks:
        new_tribe += template
tribe = new_tribe
print('\nAfter removing templates with < %d valid picks...' % min_picks)
print(tribe)

#%% If we are generating a repicked catalog, we need to exclude templates without a location
if generate_repicked_catalog:

    # Unpickle the list of unassociated template numbers for comparison
    with open(convert_redpy_output_dir + 'unassociated_clusters.txt', 'rb') as cluster_pickle:
        unassociated_clusters = pickle.load(cluster_pickle)
        #associated_clusters = pickle.load(cluster_pickle)

    # Initialize a master repicked catalog for appending in main loop later
    master_repicked_catalog = Catalog()

    # Initialize placeholder empty tribe
    new_tribe = Tribe()

    # Loop through templates
    for template in tribe:

        # Extract cluster number to compare with unassociated cluster list
        cluster_number = int(template.event.origins[0].comments[0].text.split(' ')[1])

        # Only carry over templates that are associated to the pre-existing catalog
        if cluster_number not in unassociated_clusters:
            new_tribe += template

    # Replace tribe in main workflow for next steps (the templates now all point to events with location)
    tribe = new_tribe

#%% If using local data, get a unique list of all template station-channel combinations for data fetching later
if local:

    # Initialize list for data fileheads
    data_fileheads = []

    # Loop through trace information in every template's stream
    for template in tribe:
        for trace in template.st:

            # Extract station and channel names from trace
            sta = trace.stats.station
            chan = trace.stats.channel

            # Craft filehead and append
            data_filehead = data_dir + sta + '.' + chan + '.'
            data_fileheads.append(data_filehead)

    # Now compile unique and sort
    data_fileheads = list(set(data_fileheads))
    data_fileheads.sort()

#%% Brute force scan: load each day's data and scan day-by-day

# Initialize master party and calculate number of days to scan
party_all = Party()
num_days = int(np.floor((end_time - start_time) / 86400))
time_start = time.time()

# Commence loop over days
for i in range(num_days):

    # Define temporal boundaries of our day's scan
    t1 = start_time + (i * 86400)
    t2 = start_time + ((i + 1) * 86400)
    print('\nNow at %s...' % str(t1))

    # If using local data, read in data and process stream before creating party
    if local:

        # Initialize stream object and fetch data from data_dir
        stream = Stream()

        # Loop over data fileheads, using glob.glob to match data files we want in data_dir
        for data_filehead in data_fileheads:
            data_filename = data_filehead + str(t1.year) + ':' + f'{t1.julday :03}' + ':*'
            matching_filenames = (glob.glob(data_filename))

            # Try to read the miniseed data file and add it to the stream (occasionally fails)
            for matching_filename in matching_filenames:
                try:
                    stream_contribution = read(matching_filename)
                    stream = stream + stream_contribution
                except:
                    continue

        # Process stream (remove spikes, downsample to match tribe, detrend)
        stream = remove_boxcars(stream,tolerance)
        stream = stream.resample(sampling_rate=samp_rate)
        stream = stream.detrend("simple")
        stream = stream.merge()
        print('Stream despiked, resampled and merged. Getting party of detections...')

        # Attempt to scan the current day
        try:
            party = tribe.detect(
                stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, daylong=True, overlap=None, parallel_process=True)

            # Append party to party_all
            party_all = party_all + party
            time_current = time.time()
            print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))

            # If generating a repicked catalog, we prep the stream differently for lag_calc
            if generate_repicked_catalog:

                # Copy the party and stream, then re-merge the stream using method 1
                party_calc = party.copy()
                stream_gappy = stream.copy()
                stream_gappy = stream_gappy.split()
                stream_gappy = stream_gappy.merge(method=1)

                # Use lag_calc to produce the repicked catalog section
                repicked_catalog = party_calc.lag_calc(stream_gappy, pre_processed=False, shift_len=3, min_cc=0.7,
                                                       export_cc=True, cc_dir='/home/ptan/attempt_eqcorrscan/output/cc_data')

                # Add the repicked catalog section to the master repicked catalog
                master_repicked_catalog += repicked_catalog
                time_current = time.time()
                print('Repicked catalog generated. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))

        # If the scan for the day fails, print a notification
        except:
            time_current = time.time()
            print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
            continue

    # If not using local data, directly query client
    else:

        # Define client
        client = Client(client_name)

        # Attempt scan via client downloads
        try:
            party = tribe.client_detect(client=client, starttime=t1, endtime=t2, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int)

            # Append party to party_all
            party_all = party_all + party
            time_current = time.time()
            print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))

        # If the scan for the day fails, print a notification
        except:
            time_current = time.time()
            print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
            continue

# Conclude process
time_end = time.time()
print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))
if generate_repicked_catalog:
    writer(scan_data_output_dir + repicked_catalog_filename, master_repicked_catalog)
    print('Number of relocated events with picks: %d out of %d total' % (len([1 for event in master_repicked_catalog if (event.picks != [])]), len(master_repicked_catalog)))

#%% Clean the party off of repeats (different templates that detect the "same" event)
party_declustered = party_all.decluster(trig_int=trig_int)

#%% Write our results in party form and in catalog form

# Write out party
writer(scan_data_output_dir + party_filename, party_declustered)

# Write out catalog
detected_catalog = party_declustered.get_catalog()
writer(scan_data_output_dir + catalog_filename, detected_catalog)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

#%% SENSITIVITY TEST

# This script loads in a party and rethresholds it over a range of thresholds to be explored.
# (This range of thresholds must be lower than the original party threshold.)
# The user can then choose to plot a detection value histogram, a cumulative event count plot,
# and template-detection waveform comparisons either separately or together.

# Import all dependencies
import numpy as np
import matplotlib.pyplot as plt
from eqcorrscan.utils.plotting import detection_multiplot
from obspy import UTCDateTime
from toolbox import get_detection, reader
from phase_processing.read_hypoddpha import read_hypoddpha

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
party_dir = main_dir + 'output/greatsitkin2/scan_data/'
party_filename = 'party.tgz'
PEC_dir = main_dir + 'data/avo/'
hypoi_file = 'greatsitkin_20210731_20210810_hypoi.txt'
hypoddpha_file = 'greatsitkin_20210731_20210810_hypoddpha.txt'
plot_hist = False
plot_cumu = True
plot_wave = True
separate_wave = False
thres_min = 0.60
thres_max = 0.74
thres_vec = np.linspace(thres_min,thres_max,8)
num_days = 11

# Define a base time for x-axis
base_time = UTCDateTime(2021, 7, 31, 0, 0, 0)
#base_time = UTCDateTime(2009, 1, 1, 0, 0, 0)

# Define swarm start and swarm end hours for vertical dashed lines
swarm_times = []
# swarm_times = [(UTCDateTime(2009, 2, 26, 6, 0, 0), UTCDateTime(2009, 2, 27, 13, 0, 0)),
#                (UTCDateTime(2009, 3, 20, 12, 0, 0), UTCDateTime(2009, 3, 23, 6, 34, 0)),
#                (UTCDateTime(2009, 3, 27, 0, 0, 0), UTCDateTime(2009, 3, 27, 8, 28, 0)),
#                (UTCDateTime(2009, 3, 29, 7, 50, 0), UTCDateTime(2009, 3, 29, 9, 0, 0)),
#                (UTCDateTime(2009, 4, 2, 19, 0, 0), UTCDateTime(2009, 4, 4, 13, 58, 0))]
               # (UTCDateTime(2009,5,2,21,0,0),UTCDateTime(2009,5,8,1,0,0)) May's swarm

# Define tick marks over span of party
tick_times = [UTCDateTime(2021, 7, 31), UTCDateTime(2021, 8, 2),
              UTCDateTime(2021, 8, 4), UTCDateTime(2021, 8, 6),
              UTCDateTime(2021, 8, 8), UTCDateTime(2021, 8, 10)]
# tick_times = [UTCDateTime(2009, 1, 1), UTCDateTime(2009, 1, 15),
#               UTCDateTime(2009, 2, 1), UTCDateTime(2009, 2, 15),
#               UTCDateTime(2009, 3, 1), UTCDateTime(2009, 3, 15),
#               UTCDateTime(2009, 4, 1), UTCDateTime(2009, 4, 15),
#               UTCDateTime(2009, 5, 1)]

#%% Define functions

# NIL

#%% Load party and retrieve all detect_val, detect_time and no_chans

# Load party
party = reader(party_dir + party_filename)

# Extract detection value, number of channels and detection times
all_detect_val = []
all_no_chans = []
all_detect_time = []

# Loop through families in party
for i,family in enumerate(party):

    # Loop through detections in each family
    for detection in family:

        # Append detect value, number of channels and detection time
        all_detect_val.append(detection.detect_val)
        all_no_chans.append(detection.no_chans)
        all_detect_time.append(detection.detect_time)

# Convert lists to arrays for convenience
all_detect_val = np.abs(np.array(all_detect_val))
all_no_chans = np.array(all_no_chans)
all_detect_time = np.array(all_detect_time)

# If a detection value histogram plot is desired
if plot_hist:

    # Get the average channel correlation of all detections
    all_av_chan_corr = all_detect_val / all_no_chans

    # Plot histogram
    fig, ax = plt.subplots()
    ax.hist(all_av_chan_corr, bins=40, rwidth=0.8, color='brown')  # density=False would make counts
    ax.set_ylim([0,50])
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Average Channel Correlation')
    ax.set_title('Located Detections (N=%d)' % len(all_av_chan_corr))
    ax.grid()
    fig.show()

# If a cumulative event count plot for different thresholds is desired
if plot_cumu:

    # Get AVO catalog information for plotting later
    hypoi_path = PEC_dir + hypoi_file
    hypoddpha_path = PEC_dir + hypoddpha_file
    PEC_events = read_hypoddpha(hypoi_path, hypoddpha_path, channel_convention=True)
    PEC_detect_time = [event.origins[0].time for event in PEC_events]
    PEC_detect_time = np.array(PEC_detect_time)
    PEC_detect_hours = (PEC_detect_time - base_time) / 3600

    # Get cumulative numbers for plotting PEC
    PEC_detect_num = []
    PEC_hours = []
    for i in range(num_days * 24 + 1):
        detect_num = sum(PEC_detect_hours <= i)
        PEC_detect_num.append(detect_num)
        PEC_hours.append(i)

    # Start plot
    fig, ax = plt.subplots(figsize=(9,6))

    # Loop through thresholds
    for thres in thres_vec:

        # Calculate threshold sum
        all_thresholds = thres*all_no_chans

        # Check if detections pass this adjusted threshold, and print number of detections
        valid_bool = (all_detect_val >= all_thresholds)
        print('For thres =',thres,', number of detections =',sum(valid_bool))

        # Retrieve detection times to plot, and convert to hours
        plot_detect_time = all_detect_time[valid_bool]
        plot_detect_hours = (plot_detect_time - base_time) / 3600
        plot_detect_num = []
        plot_hours = []

        # Get cumulative numbers for plotting
        for i in range(num_days*24+1):
            detect_num = sum(plot_detect_hours <= i)
            plot_detect_num.append(detect_num)
            plot_hours.append(i)

        # Plot cumulative trend in a stepwise fashion
        ax.step(plot_hours,plot_detect_num,label=str(thres))

    # Add PEC's cumulative trend in black
    ax.step(PEC_hours,PEC_detect_num,color='k',linestyle='--',linewidth=2,label='Original AVO catalog')

    # Add vertical spans indicating swarms
    for swarm_time in swarm_times:
        swarm_start_hour = (swarm_time[0] - base_time) / 3600
        swarm_end_hour = (swarm_time[1] - base_time) / 3600
        ax.axvspan(swarm_start_hour,swarm_end_hour,color='grey',alpha=0.5)

    # Tidy up plot
    ax.legend(title='Threshold',fontsize=12)
    ax.grid()
    ax.set_xlim([0,24*num_days])
    ax.set_xticks((np.array(tick_times) - base_time) / 3600)
    ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
    ax.set_ylabel('Number of detections',fontsize=15)
    ax.set_xlabel('UTC Time',fontsize=15)
    ax.set_title('Cumulative detection plot for different thresholds',fontsize=15)
    fig.show()

# Conduct waveform review if desired
if plot_wave:

    # Craft matrices of appropriate size
    max_len = 0
    for family in party:
        if len(family) > max_len:
            max_len = len(family)
    detect_val_matrix = np.zeros((len(party.families),max_len))
    detect_val_matrix[:] = np.nan
    no_chans_matrix = np.zeros((len(party.families),max_len))
    detect_time_matrix = np.zeros((len(party.families),max_len))

    # Fill up matrix with detection values, no. of chans and detection time
    for i,family in enumerate(party):
        for j,detection in enumerate(family):
            detect_val_matrix[i][j] = abs(detection.detect_val)
            no_chans_matrix[i][j] = detection.no_chans
            detect_time_matrix[i][j] = detection.detect_time

    # Sequentially compare with the different thresholds in the threshold vector
    for thres in thres_vec:

        # Get a matrix of differences
        diff_matrix = abs(detect_val_matrix - (thres*no_chans_matrix))

        # Get the index of the detection that marginally passes the threshold, and pull the detection out of the party
        marginal_index = np.nanargmin(diff_matrix[np.nonzero(diff_matrix)])
        marginal_position = (int(np.floor(marginal_index/max_len)),marginal_index%max_len)
        marginal_detection = party[marginal_position[0]][marginal_position[1]]

        # Pull out the family from which the detection belongs
        family = party[marginal_position[0]]

        # If plotting the template and detection separately
        if separate_wave:
            family.template.st.plot(equal_scale=False, size=(800, 600))
            _ = get_detection(marginal_detection, data_dir=None, client_name='IRIS', length=10, plot=True)

        # If plotting the template and detection together
        else:
            detect_stream = get_detection(marginal_detection, data_dir=None, client_name='IRIS')
            detection_multiplot(detect_stream,family.template.st.merge(),[marginal_detection.detect_time])


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

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
main_dir = '/home/ptan/enhance_catalog/'
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
                                     interpolate=False, max_workers=None, parallel_process=False)

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

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

# import packages
import pygmt
import numpy as np
import pandas as pd
from obspy import read_events, UTCDateTime

# define fixed params
VOLC_LAT = 52.0765
VOLC_LON = -176.1109
MAX_DEPTH = 15  # km
REGION = [VOLC_LON - 0.29, VOLC_LON + 0.29,
          VOLC_LAT - 0.18, VOLC_LAT + 0.18]
MAIN_WIDTH = 10
REGION_EW = [VOLC_LON - 0.3, VOLC_LON + 0.3, -5, MAX_DEPTH]
REGION_NS = [-5, MAX_DEPTH, VOLC_LAT - 0.18, VOLC_LAT + 0.18]

elev_profile_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/data/dem/'
EW_profile_filename = elev_profile_dir + 'ew_profile_GS.txt'
NS_profile_filename = elev_profile_dir + 'ns_profile_GS.txt'
EW = read_csv(EW_profile_filename)
NS = read_csv(NS_profile_filename)

# (1) Pre-existing catalog
PEC_events_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/data/avo/'
PEC_hypoi = PEC_events_dir + 'greatsitkin_20210731_20210810_hypoi.txt'
PEC_hypoddpha = PEC_events_dir + 'greatsitkin_20210731_20210810_hypoddpha.txt'
ncsn2pha(PEC_hypoi, PEC_hypoddpha, channel_convention=True)
PEC_events = read_hypoddpha(PEC_hypoi, PEC_hypoddpha, channel_convention=True)

# (2) All located templates
located_cores = reader(located_cores_filename)
unmatched_PEC_events = reader(unmatched_PEC_events_filename)
templates = located_cores + unmatched_PEC_events
located_templates = Catalog()
for template in templates:
    if template.origins[0].latitude is not None:
        located_templates += template

# (3) GrowClust relocated catalog
relocated_catalog = reader(relocated_catalog_filename)

#%% Extract information from catalog, filtered by max depth

# Choose catalog
catalog = relocated_catalog

# extract latitudes, longitudes, depths and time
latitudes = np.array([event.origins[0].latitude for event in catalog])
longitudes = np.array([event.origins[0].longitude for event in catalog])
depths = np.array([event.origins[0].depth for event in catalog]) / 1000  # km
utctimes = np.array([event.origins[0].time for event in catalog])
times = np.array([np.datetime64(event.origins[0].time) for event in catalog])
days = ((utctimes - UTCDateTime(2021,4,18,0,0,0)) / 86400).astype(int)
days_to = max(days) - days

# define projection function
def projection(region, width, cm=False):
    if cm:
        unit = 'c'
    else:
        unit = 'i'
    return 'S{}/90/{}{}'.format(np.mean(region[:2]), width, unit)

# hypocenter plot
fig = pygmt.Figure()
# bottom plot: E-W cross section
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "4c"), autolabel="c)"):
    # create basemap with correct dimensions
    fig.basemap(region=REGION_EW, projection="X10c/-4c", frame=["xa0.2f0.05+lLongitude","yaf+lDepth", "WSne"], panel=[0, 0])
    # plot landfill
    for i in range(1, len(EW)):
        delta = abs(EW["Longitude"][i] - EW["Longitude"][i - 1]) / (REGION_EW[1]-REGION_EW[0]) * 10  # inter-measurement width
        height_km = (15 + 0.001 * EW["Elevations"][i])  # depth - negative elevation
        height_cm = height_km / 20 * 3.9  # ratio * figure height
        midpoint = 15 - 0.5 * height_km
        data = [[EW["Longitude"][i], midpoint, delta, height_cm]]
        fig.plot(data=data, style="r", color="gray90", pen="1p,gray90")
    # plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2021-07-31T/2021-08-10T/1d")
    fig.plot(x=longitudes, y=depths, color=times, cmap=True, style="c0.10c", pen="black", transparency=20)
    # add colorbar
    fig.colorbar(position="JTR+o1c/-4c+w4c/0.4c+v",frame="xa1Df")
    # plot elevation profile
    fig.plot(x=EW["Longitude"], y=-0.001*EW["Elevations"], pen="1.5p,black")
# move plot origin to plot top-left plot
fig.shift_origin(yshift="h+1c")
# top-left plot: Top-down view
with fig.subplot(nrows=1, ncols=1, figsize=("10c", "10c"), autolabel="a)"):
    # plot topography
    fig.grdimage("@earth_relief_15s",region=REGION,
                 projection=projection(REGION, MAIN_WIDTH, cm=True),
                 shading=True, t=30, cmap="geo")
    # plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2021-07-31T/2021-08-10T/1d")
    fig.plot(x=longitudes, y=latitudes, style="c0.10c", color=times, cmap=True, pen="black", transparency=20)
    fig.plot(x=VOLC_LON, y=VOLC_LAT, color="red", style="t0.3c", pen="black")
    # add frame and scalebar
    pygmt.config(MAP_FRAME_TYPE="plain")
    pygmt.config(FORMAT_GEO_MAP="ddd.x")
    fig.basemap(frame = ['xa0.2f0.05','ya0.1f0.025','WseN+t"Great Sitkin"'],
                L='n0.3/0.08+w20+c{}/{}+l+f'.format(np.mean(REGION[:2]),np.mean(REGION[2:])))
    #fig.basemap(frame=0)
# move plot origin to plot top-right plot
fig.shift_origin(xshift="11c",yshift="0.04c")
# top-right plot: N-S cross section
with fig.subplot(nrows=1, ncols=1, figsize=("4c", "9.95c"), autolabel="b)"):
    # create basemap with correct dimensions
    fig.basemap(region=REGION_NS, projection="X4c/9.95c", frame=["xaf+lDepth","ya0.1f0.025+lLatitude", "wsNE"], panel=[0, 0])
    # plot landfill
    for i in range(1, len(NS)):
        delta = abs(NS["Latitude"][i] - NS["Latitude"][i - 1]) / (REGION_NS[1]-REGION_NS[0]) * 10 # inter-measurement width
        height_km = (15 + 0.001 * NS["Elevations"][i])  # depth - negative elevation
        height_cm = height_km / 20 * 3.9  # ratio * figure height
        midpoint = 15 - 0.5 * height_km
        data = [[midpoint, NS["Latitude"][i], height_cm, delta]]
        fig.plot(data=data, style="r", color="gray90", pen="1.5p,gray90")
    # plot earthquakes
    pygmt.makecpt(cmap="viridis", series="2021-07-31T/2021-08-10T/1d")
    fig.plot(x=depths, y=latitudes, style="c0.10c", color=times, cmap=True, pen="black", transparency=20)
    # plot elevation profile
    fig.plot(x=-0.001*NS["Elevations"], y=NS["Latitude"], pen="1.5p,black")
fig.show(method='external')
