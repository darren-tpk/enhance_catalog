import time
import pandas as pd
import numpy as np
from obspy import Catalog, UTCDateTime, Stream, Trace
from obspy.core.event import Event, Origin, Comment, Pick, WaveformStreamID, Arrival, ResourceIdentifier
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from toolbox import read_trace, reader, writer

# define all variables here
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks
redpy_station_list = ['RDN','REF','RSO']  # note that current code that makes picks look at EHZ only
data_dir = '/home/data/redoubt/'  # redoubt data directory

# load in detection, core and orphan lists using pandas
redpy_det = pd.read_csv('/home/ptan/attempt_eqcorrscan/avo_data/redpy/catalog.txt', sep=' ', names=['Cluster','DateTime'])
redpy_core = pd.read_csv('/home/ptan/attempt_eqcorrscan/avo_data/redpy/cores.txt', sep=' ', names=['Cluster','DateTime'])
redpy_orphan = pd.read_csv('/home/ptan/attempt_eqcorrscan/avo_data/redpy/orphancatalog.txt', names=['DateTime'])

# prepare core and orphan datetimes as an array
core_events_time = np.array([UTCDateTime(t) for t in redpy_core.DateTime])
orphan_events_time = np.array([UTCDateTime(t) for t in redpy_orphan.DateTime])

# prepare AVO catalog with known picks
catalog_dir = '/home/ptan/attempt_eqcorrscan/avo_data/'
hypoi_path = catalog_dir + hypoi_file
hypoddpha_path = catalog_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path, channel_convention=True)
avo_events = read_hypoddpha(hypoi_path, hypoddpha_path, channel_convention=True)
avo_events_time = np.array([avo_event.origins[0].time for avo_event in avo_events])

# now create a new ObsPy catalog for REDPy detections
redpy_catalog = Catalog()

# valid clusters to populate cluster event picks later
cluster_list_AVO = []  # store unique cluster numbers with AVO event within

# loop through redpy detection list and replicate clustered events in ObsPy catalog
print('\nCreating ObsPy catalog for REDPy results...')
time_start = time.time()
for i in range(len(redpy_det)):
    # extract values
    cluster = redpy_det.Cluster[i]
    detection_time = UTCDateTime(redpy_det.DateTime[i])
    # check if event is core
    if min(abs(core_events_time - detection_time)) == 0:
        redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text='CLUSTER ' + str(cluster) + ' CORE EVENT;')])])
    # otherwise it is a cluster event
    else:
        redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text='CLUSTER ' + str(cluster) + ' EVENT;')])])
    # check if event is part of AVO catalog, using a inter-event tolerance defined by user
    if min(abs(avo_events_time - detection_time)) < max_dt:
        # find the closest event in time
        avo_index = np.argmin(abs(avo_events_time - detection_time))
        avo_event = avo_events[avo_index]
        # use AVO's catalog information to fill up event object
        redpy_event.picks = avo_event.picks
        redpy_event.magnitudes = avo_event.magnitudes
        redpy_event.origins[0].longitude = avo_event.origins[0].longitude
        redpy_event.origins[0].latitude = avo_event.origins[0].latitude
        redpy_event.origins[0].depth = avo_event.origins[0].depth
        redpy_event.origins[0].arrivals = avo_event.origins[0].arrivals
        redpy_event.origins[0].comments[0].text += ' IN AVO, DT = ' + '%.2f' % min(abs(avo_events_time - detection_time))
        # add cluster number to valid clusters list
        cluster_list_AVO.append(cluster)
    # if event is not part of the AVO catalog, we use...
    else:
        redpy_event.origins[0].comments[0].text += ' NOT IN AVO'
    # finally we append the event
    redpy_catalog.append(redpy_event)
# get unique list of AVO-associated clusters and non-associated clusters
clusters_AVO = list(np.unique(np.array(cluster_list_AVO)))
clusters_NA = [cluster for cluster in list(np.unique(np.array(redpy_det.Cluster))) if cluster not in clusters_AVO]
time_stop = time.time()
print('ObsPy catalog created, processing time: %.2f s' % (time_stop-time_start))

# go through redpy_catalog and add picks to redpy detections associated with an AVO event
print('\nAdopting picks for associated cluster events...')
time_start = time.time()
for cluster_AVO in clusters_AVO:
    # find all associated detections
    adopt_list = list(np.where(np.array(redpy_det.Cluster) == cluster_AVO)[0])
    # find our reference AVO event (use core event if it has picks, if not use event with most picks in the cluster)
    all_tags = [redpy_catalog[adopt_index].origins[0].comments[0].text.split(' ')[2] for adopt_index in adopt_list]
    if 'CORE' in all_tags and redpy_catalog[adopt_list[all_tags.index('CORE')]].picks != []:
        # find the core event and extract it as a reference event
        ref_detection = redpy_catalog[adopt_list[all_tags.index('CORE')]]
        ref_detection_time = ref_detection.origins[0].time
    else:
        # find the event with the most picks and extract it as a reference event
        pick_totals = [len(redpy_catalog[adopt_index].picks) for adopt_index in adopt_list]
        ref_detection = redpy_catalog[adopt_list[np.argmax(pick_totals)]]
        ref_detection_time = ref_detection.origins[0].time
    # explicitly loop through associated detections
    for adopt_index in adopt_list:
        # determine time difference to offset pick times
        adoptee_time = redpy_catalog[adopt_index].origins[0].time
        time_difference = ref_detection_time - adoptee_time
        # duplicate reference event, remove magnitude, reproduce original comment
        ref_copy = ref_detection.copy()
        ref_copy.magnitudes = []
        ref_copy.origins[0].time = adoptee_time
        ref_copy.origins[0].comments = redpy_catalog[adopt_index].origins[0].comments
        # change pick times (account for difference in event time)
        for pick in ref_copy.picks:
            pick.time = pick.time - time_difference
        # change phase weights according to user's specification (ranges from 0 to 1)
        for arrival in ref_copy.origins[0].arrivals:
            arrival.time_weight = adopt_weight
        # update redpy catalog by replacing detection with reference event copy
        redpy_catalog[adopt_index] = ref_copy
time_stop = time.time()
print('Pick adoption complete, processing time: %.2f s' % (time_stop-time_start))

# extract a core catalog from the redpy catalog
print('\nExtracting core events from redpy catalog...')
time_start = time.time()
# initialize catalog object
all_cores = Catalog()
# loop through redpy catalog to check for cores
for detection in redpy_catalog:
    if detection.origins[0].comments[0].text.split(' ')[2] == 'CORE':
        all_cores.append(detection)
# find repeats
remove_index = []
for ii in range(len(all_cores)):
    if all_cores[ii].origins[0].time == all_cores[ii-1].origins[0].time:
        remove_index.append(ii)
# reappend events to core catalog
core_catalog = Catalog()
for ii, core in enumerate(all_cores):
    if ii not in remove_index:
        core_catalog.append(core)
time_stop = time.time()
print('Core extraction complete, processing time: %.2f s' % (time_stop-time_start))

# create a new ObsPy catalog for REDPy orphans
print('\nNow processing orphans...')
time_start = time.time()
# initialize catalog object
orphan_catalog = Catalog()
# loop through redpy orphan list and populate orphan catalog
for orphan_time in orphan_events_time:
    # create orphan event
    orphan_event = Event(origins=[Origin(time=orphan_time, comments=[Comment(text='ORPHAN EVENT;')])])
    # check if event is part of AVO catalog, using a inter-event tolerance defined by user
    if min(abs(avo_events_time - orphan_time)) < max_dt:
        # find the closest event in time
        avo_index = np.argmin(abs(avo_events_time - detection_time))
        avo_event = avo_events[avo_index]
        # use AVO's catalog information to fill up event object
        orphan_event.picks = avo_event.picks
        orphan_event.magnitudes = avo_event.magnitudes
        orphan_event.origins[0].longitude = avo_event.origins[0].longitude
        orphan_event.origins[0].latitude = avo_event.origins[0].latitude
        orphan_event.origins[0].depth = avo_event.origins[0].depth
        orphan_event.origins[0].arrivals = avo_event.origins[0].arrivals
        orphan_event.origins[0].comments[0].text += ' IN AVO, DT = ' + '%.2f' % min(abs(avo_events_time - detection_time))
    # if event is not part of the AVO catalog, we use...
    else:
        orphan_event.origins[0].comments[0].text += ' NOT IN AVO'
    # finally we append the event
    orphan_catalog.append(orphan_event)
time_stop = time.time()
print('Orphan catalog created, processing time: %.2f s' % (time_stop-time_start))

# save catalog files
catalog_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(catalog_outpath+'redpy_catalog.xml', redpy_catalog)
writer(catalog_outpath+'core_catalog.xml', core_catalog)
writer(catalog_outpath+'orphan_catalog.xml', orphan_catalog)

# # read catalog files
# catalog_outpath = '/home/ptan/attempt_eqcorrscan/output/'
# redpy_catalog = reader(catalog_outpath+'redpy_catalog.xml')
# core_catalog = reader(catalog_outpath+'core_catalog.xml')
# orphan_catalog = reader(catalog_outpath+'orphan_catalog.xml')

# catalog_outpath = '/home/ptan/attempt_eqcorrscan/output/'
# redpy_catalog = reader(catalog_outpath+'redpy_catalog.xml')
### MAKE PICKS
from obspy.signal.trigger import coincidence_trigger, classic_sta_lta, trigger_onset
# go through redpy_catalog and add picks to redpy detections associated with an AVO event
print('\nMaking picks for non-associated cluster events...')
time_start = time.time()
for cluster_NA in clusters_NA:
    # find all associated detections
    contribution_list = list(np.where(np.array(redpy_det.Cluster) == cluster_NA)[0])
    # explicitly loop through associated detections
    for contribution_index in contribution_list:
        # determine time difference to offset pick times
        contribution_time = redpy_catalog[contribution_index].origins[0].time
        # load data from local machine in a +/- 12s window
        starttime = contribution_time - 12
        endtime = contribution_time + 12
        # remove RSO after first explosion [HARD CODED]
        redpy_stations = redpy_station_list
        if starttime > UTCDateTime(2009,3,24,0,0,0):
            redpy_stations = ['RDN', 'REF']
        # gather waveforms and filter
        stream = Stream()
        for redpy_station in redpy_stations:
            station_tr = read_trace(data_dir=data_dir, station=redpy_station, channel='EHZ',
                                    starttime=starttime , endtime=endtime, tolerance=5e4)
            stream = stream + station_tr
        stream = stream.split()
        stream = stream.filter('bandpass',freqmin=1.0, freqmax=10.0, corners=2, zerophase=True)
        stream = stream.merge()
        # use coincidence trigger to get pick time estimate
        coin_trigger = coincidence_trigger('classicstalta', 3, 2, stream, 2,
                                      sta=0.7, lta=8, details=True)
        # if there are no coincidence triggers, move to next event
        if coin_trigger == []:
            continue
        # otherwise, continue to single channel STA LTA
        else:
            # extract coincidence trigger time
            coin_trigger_time = coin_trigger[0]["time"]
            # for each channel
            for tr in stream:
                # calculate the value of the characteristic function
                sampling_rate = tr.stats.sampling_rate
                cft = classic_sta_lta(tr.data, int(0.7 * sampling_rate), int(8 * sampling_rate))
                #plot_trigger(tr, cft, 3, 2)
                # obtain trigger limits
                trigger_limits = np.array(trigger_onset(cft, 3, 2))
                # if there exists some trigger limits
                if trigger_limits.size != 0:
                    # convert to UTCDateTime and find the trigger on time closest to the coincidence trigger
                    trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
                    pick_time = trigger_on[np.argmin(abs(trigger_on-coin_trigger_time))]
                    # craft waveform stream ID
                    tr_id = tr.id.split('.')
                    waveform_id = WaveformStreamID(tr_id[0],tr_id[1],tr_id[2],tr_id[3])
                    # create pick and arrival objects and add to redpy_catalog
                    add_pick = Pick(time=pick_time,waveform_id=waveform_id,phase_hint='P')
                    add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id),phase='P',time_weight=0.1)
                    redpy_catalog[contribution_index].picks.append(add_pick)
                    redpy_catalog[contribution_index].origins[0].arrivals.append(add_arrival)
time_stop = time.time()
print('Pick making complete, processing time: %.2f s' % (time_stop-time_start))

# # save catalog files
# writer(catalog_outpath+'redpy_catalog_picked.xml', redpy_catalog_picked)
# writer(catalog_outpath+'core_catalog_picked.xml', core_catalog_picked)

# # read catalog files
# redpy_catalog_picked = reader(catalog_outpath+'redpy_catalog_picked.xml')
# core_catalog_picked = reader(catalog_outpath+'core_catalog_picked.xml')

## TEMPLATES FROM STACKS?
# initialize catalog of stacks
stack_catalog = Catalog()
# for each unassociated cluster
for cluster_NA in clusters_NA:
    redpy_stations = redpy_station_list
    if cluster_NA > 560:
        redpy_stations = ['RDN','REF']
    # initialize list of empty streams
    trace_stacker = [Stream() for redpy_station in redpy_stations]
    # get a list of contributing redpy detections
    contribution_list = list(np.where(np.array(redpy_det.Cluster) == cluster_NA)[0])
    # loop through associated detections
    for contribution_index in contribution_list:
        # use a 10s window to extract data
        contribution_time = redpy_catalog[contribution_index].origins[0].time
        starttime = contribution_time - 10
        endtime = contribution_time + 10
        # gather waveforms and add to list of streams
        stream = Stream()
        for j, redpy_station in enumerate(redpy_stations):
            station_tr = read_trace(data_dir=data_dir, station=redpy_station, channel='EHZ',
                                    starttime=starttime , endtime=endtime, tolerance=5e4)
            trace_stacker[j] = trace_stacker[j] + station_tr
    # initialize stacked stream object
    stacked_stream = Stream()
    # for each station's traces
    for st in trace_stacker:
        # filter and stack
        st = st.split()
        st = st.filter('bandpass', freqmin=1.0, freqmax=10.0, corners=2, zerophase=True)
        st = st.merge()
        st = st.stack()
        stacked_stream = stacked_stream + st
    # use coincidence trigger to get pick time estimate
    coin_trigger = coincidence_trigger('classicstalta', 3, 2, stacked_stream, 2,
                                       sta=0.7, lta=8, details=True)
    # extract coincidence trigger time
    coin_trigger_time = coin_trigger[0]["time"]
    # create stacked event
    stack_event = Event(origins=[Origin(time=coin_trigger_time, comments=[Comment(text='STACKED EVENT;')])])
    # loop through each station's stacked trace
    for tr in stacked_stream:
        # calculate the value of the characteristic function
        sampling_rate = tr.stats.sampling_rate
        cft = classic_sta_lta(tr.data, int(0.7 * sampling_rate), int(8 * sampling_rate))
        #plot_trigger(tr, cft, 3, 2)
        # obtain trigger limits
        trigger_limits = np.array(trigger_onset(cft, 3, 2))
        # if there exists some trigger limits
        if trigger_limits.size != 0:
            # convert to UTCDateTime and find the trigger on time closest to the coincidence trigger
            trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
            pick_time = trigger_on[np.argmin(abs(trigger_on - coin_trigger_time))]
            # craft waveform stream ID
            tr_id = tr.id.split('.')
            waveform_id = WaveformStreamID(tr_id[0], tr_id[1], tr_id[2], tr_id[3])
            # create pick and arrival objects and add to stack event
            add_pick = Pick(time=pick_time, waveform_id=waveform_id, phase_hint='P')
            add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id), phase='P', time_weight=0.1)
            stack_event.picks.append(add_pick)
            stack_event.origins[0].arrivals.append(add_arrival)
    # append stack event
    stack_catalog.append(stack_event)
