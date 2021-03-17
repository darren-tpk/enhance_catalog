import pandas as pd
import numpy as np
from obspy import Catalog, UTCDateTime
from obspy.core.event import Event, Origin, Comment
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from toolbox import writer

# define all variables here
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks

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
cluster_list = []  # store unique cluster numbers with AVO event within

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
        cluster_list.append(cluster)
    # if event is not part of the AVO catalog, we use...
    else:
        redpy_event.origins[0].comments[0].text += ' NOT IN AVO'
    # finally we append the event
    redpy_catalog.append(redpy_event)
# get unique list of AVO-associated clusters
valid_clusters = list(np.unique(np.array(cluster_list)))
time_stop = time.time()
print('ObsPy catalog created, processing time: %.2f s' % (time_stop-time_start))

# go through redpy_catalog and add picks to redpy detections associated with an AVO event
print('\nAdopting picks for associated cluster events...')
time_start = time.time()
for valid_cluster in valid_clusters:
    # find all associated detections
    adopt_list = list(np.where(np.array(redpy_det.Cluster) == valid_cluster)[0])
    # find our reference AVO event (event with most picks in the cluster)
    pick_totals = [len(redpy_catalog[adopt_index].picks) for adopt_index in adopt_list]
    ref_detection = redpy_catalog[adopt_list[np.argmax(pick_totals)]]
    ref_detection_time = ref_detection.origins[0].time
    # explicitly loop through associated detections
    for adopt_index in adopt_list:
        # determine time difference to offset pick times
        adoptee_time = redpy_catalog[adopt_index].origins[0].time
        time_difference = ref_detection_time - adoptee_time
        # duplicate reference event and remove magnitude
        ref_copy = ref_detection.copy()
        ref_copy.origins[0].time = adoptee_time
        ref_copy.magnitudes = []
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
# initialize catalog objects
core_catalog = Catalog()  # all cores
core_loc_catalog = Catalog()  # all cores that have picks
# loop through redpy catalog to check for cores
for detection in redpy_catalog:
    if detection.origins[0].comments[0].text.split(' ')[2] == 'CORE':
        if detection.picks != []:
            core_loc_catalog.append(detection)
        core_catalog.append(detection)
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
writer(catalog_outpath+'core_loc_catalog.xml', core_loc_catalog)
writer(catalog_outpath+'orphan_catalog.xml', orphan_catalog)
