#%% CONVERT REDPY

# This script aims to convert the output from the REDPy tool into ObsPy catalog objects.
# # In the process, the REDPy output is compared to a pre-existing catalog to generate a list of associated cores.
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
from toolbox import read_trace, reader, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/redoubt/'  # redoubt data directory on local
output_dir = main_dir + 'output/redoubt5/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
redpy_results_dir = main_dir + 'redpy_results/redoubt5/'
PEC_dir = main_dir + 'data/avo/'
PEC_file = PEC_dir + 'redoubt_20080401_20090901.xml'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks
redpy_network_list = ['AV','AV','AV']
redpy_station_list = ['RDN','REF','RSO']
redpy_channel_list = [['EHZ'],['EHZ'],['EHZ']]
redpy_location_list = ['','','']
campaign_network_list = ['AV','AV','AV','AV']
campaign_station_list = ['RD01','RD02','RD03','RDW']
campaign_channel_list = [['BHZ'],['BHZ'],['BHZ'],['BHZ']]
campaign_location_list = ['','','','']
campaign_starttime = UTCDateTime(2009,3,21)
campaign_endtime = UTCDateTime(2009,7,1)
tolerance = 4e4  # tolerance for boxcar removal from data (as a factor to median)
local = True
client_name = 'IRIS' # IRIS/NCEDC
add_redpy_pick_to_associated = True
add_campaign_pick_to_associated = True
freqmin = 1
freqmax = 10

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
    for ii in range(1,len(all_cores)):
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
redpy_results_dir = main_dir + 'redpy_results/redoubt5/'
redpy_detections = pd.read_csv(redpy_results_dir + 'catalog.txt', sep=' ', names=['Cluster', 'DateTime'])
redpy_cores = pd.read_csv(redpy_results_dir + 'cores.txt', sep=' ', names=['Cluster', 'DateTime'])
redpy_orphans = pd.read_csv(redpy_results_dir + 'orphancatalog.txt', names=['DateTime'])

# Prepare core and orphan datetimes as an array
core_event_times = np.array([UTCDateTime(t) for t in redpy_cores.DateTime])
orphan_event_times = np.array([UTCDateTime(t) for t in redpy_orphans.DateTime])

# Prepare AVO catalog with known picks
PEC_events = reader(PEC_file)
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
        print('%s already exists.' % output_subdir)
print('All output subdirectories created.')

# Initialize catalog object and lists before commencing loop
redpy_catalog = Catalog()  # empty catalog to populate
associated_cluster_list = []  # store unique cluster numbers that have associated catalog events
unmatched_indices_redpy = list(range(len(PEC_events)))  # list of avo event indices. matched indices are removed
unmatched_indices_core = list(range(len(PEC_events)))  # list of avo event indices. matched indices are removed

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
        redpy_event = PEC_event.copy()
        full_event_tag = event_tag + ' IN PEC, REDPY DETECTION TIME = %s' % str(detection_time)
        redpy_event.origins[0].comments.append(Comment(text=full_event_tag))

        # Add cluster number to valid clusters list
        associated_cluster_list.append(cluster)

        # Remove avo index from unmatched list if it still exists
        if PEC_index in unmatched_indices_redpy:
            unmatched_indices_redpy.remove(PEC_index)
        if PEC_index in unmatched_indices_core and event_tag.split(' ')[2] == 'CORE':
            unmatched_indices_core.remove(PEC_index)

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
with open(convert_redpy_output_dir + 'unmatched_indices_redpy.txt', 'wb') as unmatched_pickle:  # Pickling
    pickle.dump(unmatched_indices_redpy, unmatched_pickle)
with open(convert_redpy_output_dir + 'unmatched_indices_core.txt', 'wb') as unmatched_pickle:  # Pickling
    pickle.dump(unmatched_indices_core, unmatched_pickle)

# Also generate unmatched PEC catalogs (that can be used as templates later)
unmatched_PEC_events_redpy = Catalog() # PEC events that don't match redpy catalog
unmatched_PEC_events_core = Catalog() # PEC events that don't match redpy cores
for j, PEC_event in enumerate(PEC_events):
    if j in unmatched_indices_redpy:
        unmatched_PEC_events_redpy += PEC_event
    if j in unmatched_indices_core:
        unmatched_PEC_events_core += PEC_event

# Write out unmatched PEC catalog to .xml file
writer(convert_redpy_output_dir + 'unmatched_PEC_events_redpy.xml', unmatched_PEC_events_redpy)
writer(convert_redpy_output_dir + 'unmatched_PEC_events_core.xml', unmatched_PEC_events_core)

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
        orphan_event = PEC_event.copy()
        full_event_tag = event_tag + ' IN PEC, ORPHAN DETECTION TIME = %s' % str(orphan_time)
        orphan_event.origins[0].comments.append(Comment(text=full_event_tag))

    # If event is not part of the PEC, we tag it otherwise
    else:
        orphan_event.origins[0].comments[0].text += ' NOT IN PEC'

    # Finally we append the event
    orphan_catalog.append(orphan_event)

# Write the fully picked redpy catalog to an xml file
writer(convert_redpy_output_dir + 'orphan_catalog.xml', orphan_catalog)

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
        stream = stream.filter('bandpass',freqmin=freqmin, freqmax=freqmax, corners=2, zerophase=True)
        stream = stream.taper(0.05, type='hann', max_length=(0.75*1024/100))  # [HARD CODED]
        stream = stream.merge(method=1, fill_value=0)
        stream = stream.trim(starttime=starttime,endtime=endtime)

        # Use coincidence trigger to get a pick time estimate
        try:
            coin_trigger = []
            coincidence_sum_threshold = 2 # int(round(len(redpy_stations)/2))
            coin_trigger = coincidence_trigger('classicstalta', 3, 2, stream, coincidence_sum_threshold, sta=0.7, lta=8, details=True)
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
                    add_pick = Pick(time=pick_time,waveform_id=waveform_id,phase_hint='P',comments=[Comment(text='DERIVED FROM REDPY STA/LTA')])
                    # add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id),phase='P',time_weight=0.1)
                    redpy_catalog[contribution_index].picks.append(add_pick)
                    # redpy_catalog[contribution_index].origins[0].arrivals.append(add_arrival)

if add_redpy_pick_to_associated:

    # %% For associated clusters, give their contribution list redpy station picks based on their closest event in time

    # Commence process
    print('Translating redpy picks for associated cluster events...')
    time_start = time.time()

    # Loop through associated clusters
    for associated_cluster in associated_clusters:

        # Find all detections within associated cluster
        contribution_list = list(np.where(np.array(redpy_detections.Cluster) == associated_cluster)[0])

        # Find indices of events with picks and without picks within cluster
        contribution_indices_with_picks = []
        contribution_indices_without_picks = []
        for contribution_index in contribution_list:
            if redpy_catalog[contribution_index].picks == []:
                contribution_indices_without_picks.append(contribution_index)
            else:
                contribution_indices_with_picks.append(contribution_index)

        # Loop over event indices that need picks
        for contribution_index in contribution_indices_without_picks:

            # Find the pick sharer event and use its REDPy detection time to determine pick offset
            contribution_time = redpy_catalog[contribution_index].origins[0].time
            contribution_time_diffs = [(UTCDateTime(redpy_catalog[ind].origins[0].comments[0].text.split()[-1])
                                        - contribution_time) for ind in contribution_indices_with_picks]
            contribution_time_diffs_abs = [abs(diff) for diff in contribution_time_diffs]
            pick_sharer_event = redpy_catalog[contribution_indices_with_picks[np.argmin(contribution_time_diffs_abs)]]
            contribution_time_diff = contribution_time_diffs[np.argmin(contribution_time_diffs_abs)]

            # Loop through picks and give picks if it is from a redpy station (translating it in time)
            for pick in pick_sharer_event.picks:
                if pick.waveform_id.station_code in redpy_stations:
                    add_pick = pick.copy()
                    add_pick.time -= contribution_time_diff
                    redpy_catalog[contribution_index].picks.append(add_pick)

if add_campaign_pick_to_associated:

    # Commence process
    print('Adding campaign picks for associated cluster events...')
    time_start = time.time()

    # Define client if using data from server
    if not local:
        client = Client(client_name)

    # Define some coincidence_trigger arguments

    # Loop through associated clusters
    for associated_cluster in associated_clusters:

        # Find all detections within unassociated cluster
        contribution_list = list(np.where(np.array(redpy_detections.Cluster) == associated_cluster)[0])

        # Loop through these unassociated detections to add picks
        for contribution_index in contribution_list:

            # Extract contribution time
            contribution_time = redpy_catalog[contribution_index].origins[0].time

            # Retrieve time limits for data fetch (+/- 12s window)
            starttime = contribution_time - 12
            endtime = contribution_time + 12

            # Remove stations if necessary
            campaign_stations = campaign_station_list

            # if the contribution is from before the campaign was installed, we skip the pick-adding for the event
            if contribution_time < campaign_starttime or contribution_time > campaign_endtime:
                continue

            # Gather data from local machine in a +/- 12s window, filter, and taper
            stream = Stream()
            for k, campaign_station in enumerate(campaign_stations):
                for campaign_channel in campaign_channel_list[k]:
                    if local:
                        station_tr = read_trace(data_dir=data_dir, station=campaign_station, channel=campaign_channel,
                                                starttime=starttime, endtime=endtime, tolerance=tolerance)
                    else:
                        station_tr = client.get_waveforms(campaign_network_list[k], campaign_station, campaign_location_list[k],
                                                          campaign_channel, starttime, endtime)
                    stream = stream + station_tr
            stream = stream.split()
            stream = stream.filter('bandpass', freqmin=1.0, freqmax=10.0, corners=2, zerophase=True)
            stream = stream.taper(0.05, type='hann', max_length=(0.75 * 1024 / 100))  # [HARD CODED]
            stream = stream.merge(method=1, fill_value=0)
            stream = stream.trim(starttime=starttime, endtime=endtime)

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

                # Remove traces that belong to channels that already have been picked
                existing_station_list = [pick.waveform_id.station_code for pick in redpy_catalog[contribution_index].picks]
                for tr in stream:
                    if tr.stats.station in existing_station_list:
                        stream.remove(tr)

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
                        pick_time = trigger_on[np.argmin(abs(trigger_on - coin_trigger_time))]

                        # Craft waveform stream ID
                        tr_id = tr.id.split('.')
                        waveform_id = WaveformStreamID(tr_id[0], tr_id[1], tr_id[2], tr_id[3])

                        # Create ObsPy pick and ObsPy arrival objects and add to redpy_catalog
                        add_pick = Pick(time=pick_time, waveform_id=waveform_id, phase_hint='P',comments=[Comment(text='DERIVED FROM CAMPAIGN STA/LTA')])
                        # add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id), phase='P',
                        #                       time_weight=0.1)
                        redpy_catalog[contribution_index].picks.append(add_pick)
                        # redpy_catalog[contribution_index].origins[0].arrivals.append(add_arrival)

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
