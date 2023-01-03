"""
Functions to run REDPy and convert REDPy output files into desired input for EQcorrscan.

Specifically, it compares the REDPy output with the pre-existing catalog (PEC),
and gives the user the option of producing either
(a) a consolidated catalog including all REDPy cores and all PEC events
(b) a truncated catalog of REDPy cores and all *unincluded* PEC events, where cores are prioritized
(c) same as (b), but the highest magnitude PEC event in the REDPy cluster is prioritized
"""

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

def execute_redpy():
    pass

def convert_redpy(pec_filepath,
                  redpy_catalog_filepath,
                  redpy_cores_filepath,
                  clusters_to_remove,
                  max_dt,
                  make_stalta_picks,
                  data_source):
    """

    :param redpy_catalog_filepath:
    :param redpy_core_filepath:
    :return:
    """
    # Read pec filepath as a catalog object
    pec_events = reader(pec_filepath)
    # Read redpy text files using pandas
    redpy_detections = pd.read_csv(redpy_catalog_filepath', sep=' ', names=['Cluster', 'DateTime'])
    redpy_cores = pd.read_csv(redpy_cores_filepath, sep=' ', names=['Cluster', 'DateTime'])
    # Remove bad clusters after manual check (highly encouraged)
    cluster_index_to_remove = redpy_cores.index[redpy_cores.Cluster.isin(clusters_to_remove)].tolist()
    redpy_cores = redpy_cores.drop(cluster_index_to_remove)
    redpy_cores = redpy_cores.reset_index(drop=True)
    detection_index_to_remove = redpy_detections.index[redpy_detections.Cluster.isin(clusters_to_remove)].tolist()
    redpy_detections = redpy_detections.drop(detection_index_to_remove)
    redpy_detections = redpy_detections.reset_index(drop=True)
    # Create arrays storing pec event times and core event times
    core_event_times = np.array([UTCDateTime(t) for t in redpy_cores.DateTime])
    pec_event_times = np.array([pec_event.origins[0].time for pec_event in pec_events])
    # Initialize catalog object and lists before commencing loop
    template_catalog = Catalog() # empty catalog to populate
    associated_cluster_list = [] # store unique cluster numbers that have associated catalog events
    unmatched_indices_redpy = list(range(len(PEC_events))) # list of pec event indices, matched indices are removed later
    unmatched_indices_core = list(range(len(PEC_events)))  # list of pec event indices, matched indices are removed later
    # Loop through redpy detection list and construct catalog
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
        if min(abs(pec_event_times - detection_time)) < max_dt:
            # Find the closest event in time from the pre-existing catalog
            pec_index = np.argmin(abs(pec_event_times - detection_time))
            pec_event = pec_events[pec_index]
            # Use the closest event's information to fill up event object
            redpy_event.picks = pec_event.picks
            redpy_event.magnitudes = pec_event.magnitudes
            redpy_event.origins[0].longitude = pec_event.origins[0].longitude
            redpy_event.origins[0].latitude = pec_event.origins[0].latitude
            redpy_event.origins[0].depth = pec_event.origins[0].depth
            redpy_event.origins[0].arrivals = pec_event.origins[0].arrivals
            redpy_event.origins[0].comments[0].text += ' IN PEC, DT = ' + '%.2f' % min(abs(PEC_event_times - detection_time))
            # Add cluster number to valid clusters list
            associated_cluster_list.append(cluster)
            # Remove avo index from unmatched list if it still exists
            if pec_index in unmatched_indices_redpy:
                unmatched_indices_redpy.remove(pec_index)
            if pec_index in unmatched_indices_core and event_tag.split(' ')[2] == 'CORE':
                unmatched_indices_core.remove(pec_index)
        # If event is not part of the PEC, we tag it as so
        else:
            redpy_event.origins[0].comments[0].text += ' NOT IN PEC'
        # Finally we append the event
        redpy_catalog.append(redpy_event)
    # Get unique list of AVO-associated clusters and non-associated clusters
    associated_clusters = list(np.unique(np.array(associated_cluster_list)))
    unassociated_clusters = [cluster for cluster in list(np.unique(np.array(redpy_detections.Cluster)))
                             if cluster not in associated_clusters]
    # Also generate unmatched pec catalogs (that can be used as templates later)
    unmatched_pec_events_redpy = Catalog()  # pec events that don't match redpy catalog
    unmatched_pec_events_core = Catalog()  # pec events that don't match redpy cores
    for j, pec_event in enumerate(pec_events):
        if j in unmatched_indices_redpy:
            unmatched_pec_events_redpy += pec_event
        if j in unmatched_indices_core:
            unmatched_pec_events_core += pec_event
    # Pull core catalog from redpy catalog
    core_catalog = _pull_cores(redpy_catalog)
    # Add STA/LTA picks for non-associated cluster events
    if make_stalta_picks:
        pass

def _pull_cores(full_catalog):
    """
    Extracts a catalog of cores included within a full catalog through comment checking
    :param full_catalog: input Catalog object
    :return: core_catalog: output Catalog object
    """
    # Initialize core catalog object
    core_catalog = Catalog()
    # Loop through full catalog to check for cores
    for detection in full_catalog:
        if detection.origins[0].comments[0].text.split(' ')[2] == 'CORE':
            core_catalog.append(detection)
    # Find repeats and remove them
    for i in reversed(range(1,len(core_catalog))):
        if core_catalog[i].origins[0].time == core_catalog[i-1].origins[0].time:
            core_catalog.events.remove(core_catalog[i])
    # Return core catalog
    return core_catalog

