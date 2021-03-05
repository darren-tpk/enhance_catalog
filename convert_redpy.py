import pandas as pd
import numpy as np
from obspy import Catalog, UTCDateTime
from obspy.core.event import Event, Origin, Comment
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha

# define all variables here
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds

# load in detection list and core list using pandas
redpy_det = pd.read_csv('/home/ptan/project/avo_data/redpy/catalog.txt', sep=' ', names=['Cluster','DateTime'])
redpy_core = pd.read_csv('/home/ptan/project/avo_data/redpy/cores.txt', sep=' ', names=['Cluster','DateTime'])
redpy_orphan = pd.read_csv('/home/ptan/project/avo_data/redpy/orphancatalog.txt', names=['DateTime'])
# prepare core and orphan datetimes as an array
core_events_time = np.array([UTCDateTime(t) for t in redpy_core.DateTime])
orphan_events_time = np.array([UTCDateTime(t) for t in redpy_orphan.DateTime])
# prepare AVO catalog with known picks
catalog_dir = '/home/ptan/project/avo_data/'
hypoi_path = catalog_dir + hypoi_file
hypoddpha_path = catalog_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path)
avo_events = read_hypoddpha(hypoi_path, hypoddpha_path)
avo_events_time = np.array([avo_event.origins[0].time for avo_event in avo_events])
# now create a new ObsPy catalog for REDPy detections and cores
redpy_catalog = Catalog()
core_catalog = Catalog()
for i in range(len(redpy_det)):
    # set core boolean to false
    core = False
    # extract values
    cluster = redpy_det.Cluster[i]
    det_time = UTCDateTime(redpy_det.DateTime[i])
    # check if event is core
    if min(abs(core_events_time - det_time)) == 0:
        redpy_event = Event(origins=[Origin(time=det_time,comments=[Comment(text='CORE EVENT;')])])
        core = True # to append to core catalog as well
    # check if event is an orphan
    elif min(abs(orphan_events_time - det_time)) == 0:
        redpy_event = Event(origins=[Origin(time=det_time, comments=[Comment(text='ORPHAN EVENT;')])])
    # otherwise it is a cluster event
    else:
        redpy_event = Event(origins=[Origin(time=det_time, comments=[Comment(text='CLUSTER ' + str(cluster) + ' EVENT;')])])
    # check if event is part of AVO catalog, using a inter-event tolerance of 10.5s
    if min(abs(avo_events_time - det_time)) < max_dt:
        # find the closest event in time
        avo_index = np.argmin(abs(avo_events_time - det_time))
        avo_event = avo_events[avo_index]
        # use AVO's catalog information to fill up event object
        redpy_event.picks = avo_event.picks
        redpy_event.magnitudes = avo_event.magnitudes
        redpy_event.origins[0].longitude = avo_event.origins[0].longitude
        redpy_event.origins[0].latitude = avo_event.origins[0].latitude
        redpy_event.origins[0].depth = avo_event.origins[0].depth
        redpy_event.origins[0].arrivals = avo_event.origins[0].arrivals
        redpy_event.origins[0].comments[0].text += ' IN AVO, DT = ' + '%.2f' % min(abs(avo_events_time - det_time))
    # if event is not part of the AVO catalog, we use...
    else:
        redpy_event.origins[0].comments[0].text += ' NOT IN AVO'
    # finally we append the event
    redpy_catalog.append(redpy_event)
    if core:
        core_catalog.append(redpy_event)
