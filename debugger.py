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
import pandas as pd
import subprocess

from eqcorrscan.core.lag_calc import _concatenate_and_correlate, _xcorr_interp
from eqcorrscan.utils.catalog_to_dd import write_correlations, _compute_dt_correlations, _EventPair, _DTObs
from obspy import Catalog, UTCDateTime
from obspy.core.event import Origin, Event
from toolbox import reader, prepare_stream_dict, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/ptan/enhance_catalog/data/mammoth/'
output_dir = main_dir + 'output/mammoth/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
catalog_filename = 'party_catalog.xml'
relocate_catalog_output_dir = output_dir + 'relocate_catalog/'
relocated_catalog_filename = 'relocated_catalog.xml'
raw_station_list_dir = main_dir + 'data/stations/'
raw_station_list_filename = 'mammoth_station_list.csv'
raw_vzmodel_dir = main_dir + 'data/vz/'
raw_vzmodel_filename = 'mammoth_vzmodel.txt'
ratio_provided = False
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
min_cc = 0.7            # minimum cc to be considered pick

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
stream_dict = prepare_stream_dict(catalog,pre_pick=pre_pick_excess,length=length_excess,local=True,data_dir=data_dir)

catalog=catalog
stream_dict=stream_dict
extract_len=length_actual
pre_pick=pre_pick_actual
shift_len=shift_len
lowcut=lowcut
highcut=highcut
max_sep=max_sep
min_link=min_link
min_cc=min_cc
interpolate=False
max_workers=None
parallel_process=False

def _filter_stream(event_id, st, lowcut, highcut):
    if lowcut is not None and highcut is not None:
        st_out = st.copy().detrend().filter(
            "bandpass", freqmin=lowcut, freqmax=highcut, corners=4,
            zerophase=True)
    elif lowcut is None and highcut is not None:
        st_out = st.copy().detrend().filter(
            "lowpass", freq=highcut, corners=4, zerophase=True)
    elif lowcut is not None and highcut is None:
        st_out = st.copy().detrend().filter(
            "highpass", freq=lowcut, corners=4, zerophase=True)
    else:
        st_out = st  # Don't need to copy if we aren't doing anything.
    return {event_id: st_out}

def _meta_filter_stream(event_id, stream_dict, lowcut, highcut):
    return _filter_stream(
        event_id=event_id, st=stream_dict[event_id], lowcut=lowcut,
        highcut=highcut)

# Process the streams
processed_stream_dict = dict()
for key in stream_dict.keys():
    processed_stream_dict.update(_meta_filter_stream(
        stream_dict=stream_dict, lowcut=lowcut, highcut=highcut,
        event_id=key))


include_master = False
correlation_kwargs = dict(
    min_cc=min_cc, stream_dict=stream_dict, extract_len=extract_len,
    pre_pick=pre_pick, shift_len=shift_len, interpolate=interpolate,
    max_workers=max_workers)

correlation = True

if correlation:
    for arg, name in correlation_kwargs.items():
        assert arg is not None, "{0} is required for correlation".format(
            name)


def _generate_event_id_mapper(catalog, event_id_mapper=None):
    event_id_mapper = event_id_mapper or dict()
    try:
        largest_event_id = max(event_id_mapper.values())
    except ValueError:
        largest_event_id = 0
    for event in catalog:
        if event.resource_id.id not in event_id_mapper.keys():
            event_id = largest_event_id + 1
            largest_event_id = event_id
            event_id_mapper.update({event.resource_id.id: event_id})
    return event_id_mapper

# Ensure all events have locations and picks.
event_id_mapper = _generate_event_id_mapper(
    catalog=catalog, event_id_mapper=None)


from eqcorrscan.utils.clustering import dist_mat_km
distances = dist_mat_km(catalog)
import numpy as np
distance_filter = distances <= max_sep
if not include_master:
    np.fill_diagonal(distance_filter, 0)


additional_args = dict(min_link=min_link, event_id_mapper=event_id_mapper)

differential_times = {}
additional_args.update(correlation_kwargs)
n = len(catalog)

from collections import namedtuple, defaultdict, Counter
SeedPickID = namedtuple("SeedPickID", ["seed_id", "phase_hint"])
from obspy import Stream

def _prepare_stream(stream, event, extract_len, pre_pick, seed_pick_ids=None):
    seed_pick_ids = seed_pick_ids or {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in event.picks if pick.phase_hint.startswith(("P", "S"))}
    stream_sliced = defaultdict(lambda: Stream())
    for seed_pick_id in seed_pick_ids:
        pick = [pick for pick in event.picks
                if pick.waveform_id.get_seed_string() == seed_pick_id.seed_id
                and pick.phase_hint[0] == seed_pick_id.phase_hint]
        if len(pick) > 1:
            pick = sorted(pick, key=lambda p: p.time)
        elif len(pick) == 0:
            continue
        pick = pick[0]
        tr = stream.select(id=seed_pick_id.seed_id).merge()
        if len(tr) == 0:
            continue
        else:
            tr = tr[0]
        tr = stream.select(id=seed_pick_id.seed_id).slice(
            starttime=pick.time - pre_pick,
            endtime=(pick.time - pre_pick) + extract_len).merge()
        if len(tr) == 0:
            continue
        if len(tr) > 1:
            continue
        tr = tr[0]
        stream_sliced.update(
            {seed_pick_id.phase_hint:
             stream_sliced[seed_pick_id.phase_hint] + tr})
    return stream_sliced

for i, master in enumerate(catalog):

    master_id = master.resource_id.id
    sub_catalog = [ev for j, ev in enumerate(catalog)
                   if distance_filter[i][j]]

    differential_times.update({
        master_id: _compute_dt_correlations(
            sub_catalog, master, **additional_args)})

    ### continue here
    max_workers = 1

    differential_times_dict = dict()
    master_stream = _prepare_stream(
        stream=stream_dict[master.resource_id.id], event=master,
        extract_len=extract_len, pre_pick=pre_pick)

    available_seed_ids = {tr.id for st in master_stream.values() for tr in st}


    master_seed_ids = {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in master.picks if
        pick.phase_hint[0] in "PS" and
        pick.waveform_id.get_seed_string() in available_seed_ids}

    # Dictionary of travel-times for master keyed by {station}_{phase_hint}
    master_tts = dict()
    master_origin_time = (master.preferred_origin() or master.origins[0]).time
    for pick in master.picks:
        if pick.phase_hint[0] not in "PS":
            continue
        tt1 = pick.time - master_origin_time
        master_tts.update({
            "{0}_{1}".format(
                pick.waveform_id.station_code, pick.phase_hint[0]): tt1})

    matched_length = extract_len + (2 * shift_len)
    matched_pre_pick = pre_pick + shift_len
    # We will use this to maintain order
    event_dict = {event.resource_id.id: event for event in sub_catalog}
    event_ids = set(event_dict.keys())
    # Check for overlap
    _stream_event_ids = set(stream_dict.keys())

    matched_streams = {
        event_id: _prepare_stream(
            stream=stream_dict[event_id], event=event_dict[event_id],
            extract_len=matched_length, pre_pick=matched_pre_pick,
            seed_pick_ids=master_seed_ids)
        for event_id in event_ids}

    sampling_rates = {tr.stats.sampling_rate for st in master_stream.values()
                      for tr in st}
    for phase_hint in master_stream.keys():  # Loop over P and S separately
        for sampling_rate in sampling_rates:  # Loop over separate samp rates
            delta = 1.0 / sampling_rate
            _master_stream = master_stream[phase_hint].select(
                sampling_rate=sampling_rate)
            _matched_streams = dict()
            for key, value in matched_streams.items():
                _st = value[phase_hint].select(sampling_rate=sampling_rate)
                if len(_st) > 0:
                    _matched_streams.update({key: _st})

            # Check lengths
            master_length = [tr.stats.npts for tr in _master_stream]

            master_length = Counter(master_length).most_common(1)[0][0]
            _master_stream = _master_stream.select(npts=master_length)
            matched_length = Counter(
                (tr.stats.npts for st in _matched_streams.values()
                 for tr in st))

            matched_length = matched_length.most_common(1)[0][0]

            # Remove empty streams and generate an ordered list of event_ids
            used_event_ids, used_matched_streams = [], []
            for event_id, _matched_stream in _matched_streams.items():
                _matched_stream = _matched_stream.select(npts=matched_length)
                if len(_matched_stream) > 0:
                    used_event_ids.append(event_id)
                    used_matched_streams.append(_matched_stream)
            # Check that there are matching seed ids.
            master_seed_ids = set(tr.id for tr in _master_stream)
            matched_seed_ids = set(
                tr.id for st in used_matched_streams for tr in st)


            ccc_out, used_chans = _concatenate_and_correlate(
                template=_master_stream, streams=used_matched_streams,
                cores=max_workers)
            # Convert ccc_out to pick-time
            for i, used_event_id in enumerate(used_event_ids):
                for j, chan in enumerate(used_chans[i]):
                    if not chan.used:
                        continue
                    correlation = ccc_out[i][j]
                    if interpolate:
                        shift, cc_max = _xcorr_interp(correlation, dt=delta)
                    else:
                        cc_max = np.amax(correlation)
                        shift = np.argmax(correlation) * delta
                    if cc_max < min_cc:
                        continue
                    shift -= shift_len
                    pick = [p for p in event_dict[used_event_id].picks
                            if p.phase_hint[0] == phase_hint
                            and p.waveform_id.station_code == chan.channel[0]
                            and p.waveform_id.channel_code == chan.channel[1]]
                    pick = sorted(pick, key=lambda p: p.time)[0]
                    tt2 = pick.time - (
                            event_dict[used_event_id].preferred_origin() or
                            event_dict[used_event_id].origins[0]).time
                    tt2 += shift
                    diff_time = differential_times_dict.get(
                        used_event_id, None)
                    if diff_time is None:
                        diff_time = _EventPair(
                            event_id_1=event_id_mapper[master.resource_id.id],
                            event_id_2=event_id_mapper[used_event_id])
                    diff_time.obs.append(
                        _DTObs(station=chan.channel[0],
                               tt1=master_tts["{0}_{1}".format(
                                   chan.channel[0], phase_hint)],
                               tt2=tt2, weight=cc_max ** 2,
                               phase=phase_hint[0]))
                    differential_times_dict.update({used_event_id: diff_time})
    # Threshold on min_link
    differential_times = [dt for dt in differential_times_dict.values()
                          if len(dt.obs) >= min_link]

# Remove Nones
for key, value in differential_times.items():
    differential_times.update({key: [v for v in value if v is not None]})

correlation_times = differential_times


with open("dt2.cc", "w") as f:
    for master_id, linked_events in correlation_times.items():
        for linked_event in linked_events:
            f.write(linked_event.cc_string)
            f.write("\n")
return event_id_mapper