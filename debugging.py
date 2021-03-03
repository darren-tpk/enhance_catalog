# Scan data to make detections
# Data is from 20km radius around Augustine

# Import packages and functions that we need
import pickle
from collections import namedtuple

import numpy as np
from eqcorrscan import Party, get_stream_xcorr
from eqcorrscan.core.match_filter.helpers import get_waveform_client
from obspy import Catalog, Stream, Trace
from eqcorrscan.core.match_filter.template import group_templates

# Load one party's data


save_dir = '/home/ptan/project/output/'
party_file = save_dir + 'sample_party.pickle'
party_in = open(party_file,'rb')
party = pickle.load(party_in)
st = pickle.load(party_in)
print(party)
print(st)

party = Party() + party[3]

import logging
Logger = logging.getLogger(__name__)

#  party.lag_calc params
stream = st.merge(method=1)
shift_len = 0.5
min_cc = 0.4
horizontal_chans = ['E','N','1','2','Z']
vertical_chans = ['E','N','1','2','Z']
cores = 1
interpolate = False
plot = False
plotdir = False
pre_processed = False
process_cores = None
parallel = True
ignore_bad_data = False
ignore_length = False

# party.lag_calc

process_cores = process_cores or cores
template_groups = group_templates(
    [_f.template for _f in party.families
     if len(_f) > 0])  # Fix for #341
catalog = Catalog()

template_group = template_groups[0]
family = [_f for _f in party.families
          if _f.template == template_group[0]][0]
group_seed_ids = {tr.id for template in template_group
                  for tr in template.st}

template_stream = Stream()
for seed_id in group_seed_ids:
    net, sta, loc, chan = seed_id.split('.')
    template_stream += stream.select(
        network=net, station=sta, location=loc, channel=chan)

processed_stream = family._process_streams(
    stream=template_stream, pre_processed=pre_processed,
    process_cores=process_cores, parallel=parallel,
    ignore_bad_data=ignore_bad_data, ignore_length=ignore_length,
    select_used_chans=False)

template = template_group[0]
family = [_f for _f in party.families
          if _f.template == template][0]
catalog += family.lag_calc(
    stream=processed_stream, pre_processed=True,
    shift_len=shift_len, min_cc=min_cc,
    horizontal_chans=horizontal_chans,
    vertical_chans=vertical_chans, cores=cores,
    interpolate=interpolate, plot=plot, plotdir=plotdir,
    parallel=parallel, process_cores=process_cores,
    ignore_bad_data=ignore_bad_data,
    ignore_length=ignore_length, **kwargs)

# family.lag_calc
from eqcorrscan.core.lag_calc import xcorr_pick_family, _prepare_data, _concatenate_and_correlate

processed_stream = family._process_streams(
    stream=stream, pre_processed=pre_processed,
    process_cores=process_cores, parallel=parallel,
    ignore_bad_data=ignore_bad_data, ignore_length=ignore_length)

# xcorr_pick_family
# picked_dict = xcorr_pick_family(
#     family=family, stream=processed_stream, shift_len=shift_len,
#     min_cc=min_cc, horizontal_chans=horizontal_chans,
#     vertical_chans=vertical_chans, cores=cores,
#     interpolate=interpolate, plot=plot, plotdir=plotdir)

picked_dict = {}
delta = family.template.st[0].stats.delta
detect_streams_dict = _prepare_data(
    family=family, detect_data=stream, shift_len=shift_len)
detection_ids = list(detect_streams_dict.keys())
detect_streams = [detect_streams_dict[detection_id]
                  for detection_id in detection_ids]

if len(detect_streams) == 0:
    Logger.warning("No appropriate data found, check your family and "
                   "detections - make sure seed ids match")
if len(detect_streams) != len(family):
    Logger.warning("Not all detections have matching data. "
                   "Proceeding anyway. HINT: Make sure SEED IDs match")

# concatenate and correlate
ccc, chans = _concatenate_and_correlate(streams=detect_streams, template=family.template.st, cores=cores)


### NEW ATTEMPT
for family in party:
    for detection in family:
        detection._calculate_event(template=family.template)
party.lag_calc(stream, pre_processed=False, shift_len=0.5, min_cc=0.4,horizontal_chans = ['E','N','1','2','Z'], vertical_chans = ['E','N','1','2','Z'])


### Comparing client approach with raw data approach
from obspy.clients.fdsn import Client
client = Client('IRIS')
start_time = UTCDateTime(2009, 2, 26, 0, 0, 0)
i = 0
t1 = start_time + (i * 86400)
t2 = start_time + ((i + 1) * 86400)
starttime = t1
endtime = t2
threshold = 20
threshold_type = "MAD"
trig_int = 30
plot = False
plotdir = None
min_gap = None
daylong = False
parallel_process = True
xcorr_func = None
concurrency = None
cores = None
ignore_length = False
ignore_bad_data = False
group_size = None
return_stream = True
full_peaks = False
save_progress = False
process_cores = None
retries = 3

from obspy.clients.fdsn.client import FDSNException
if not hasattr(client, "get_waveforms_bulk"):
    assert hasattr(client, "get_waveforms"), (
        f"client {client} must have at least a get_waveforms method")
    Logger.info(f"Client {client} does not have a get_waveforms_bulk "
                "method, monkey-patching this")
    client = get_waveform_client(client)

party = Party()
buff = 300
data_length = max([t.process_length for t in self.templates])
pad = 0
starttime = UTCDateTime(2009, 2, 26, 0, 0, 0)
endtime = start_time + (24 * 60 * 60)
return_stream = True

for template in self.templates:
    max_delay = (template.st.sort(['starttime'])[-1].stats.starttime -
                 template.st.sort(['starttime'])[0].stats.starttime)
    if max_delay > pad:
        pad = max_delay # important
download_groups = int(endtime - starttime) / data_length
template_channel_ids = []
for template in self.templates:
    for tr in template.st:
        if tr.stats.network not in [None, '']:
            chan_id = (tr.stats.network,)
        else:
            chan_id = ('*',)
        if tr.stats.station not in [None, '']:
            chan_id += (tr.stats.station,)
        else:
            chan_id += ('*',)
        if tr.stats.location not in [None, '']:
            chan_id += (tr.stats.location,)
        else:
            chan_id += ('*',)
        if tr.stats.channel not in [None, '']:
            if len(tr.stats.channel) == 2:
                chan_id += (tr.stats.channel[0] + '?' +
                            tr.stats.channel[-1],)
            else:
                chan_id += (tr.stats.channel,)
        else:
            chan_id += ('*',)
        template_channel_ids.append(chan_id)

        template_channel_ids = list(set(template_channel_ids))
        if return_stream:
            stream = Stream()
        if int(download_groups) < download_groups:
            download_groups = int(download_groups) + 1
        else:
            download_groups = int(download_groups)
         ##########
        for i in range(download_groups):
            bulk_info = []
            for chan_id in template_channel_ids:
                bulk_info.append((
                    chan_id[0], chan_id[1], chan_id[2], chan_id[3],
                    starttime + (i * data_length) - (pad + buff),
                    starttime + ((i + 1) * data_length) + (pad + buff)))
            for retry_attempt in range(retries):
                try:
                    st = client.get_waveforms_bulk(bulk_info)
                    break
                except FDSNException as e:
                    if "Split the request in smaller" in " ".join(e.args):
                        st = Stream()
                        for _bulk in bulk_info:
                            try:
                                st += client.get_waveforms_bulk([_bulk])
                            except Exception as e:
                                continue
                        break
                except Exception as e:
                    continue
            else:
                raise MatchFilterError(
                    "Could not download data after {0} attempts".format(
                        retries))
            # Get gaps and remove traces as necessary
            if min_gap:
                gaps = st.get_gaps(min_gap=min_gap)
                if len(gaps) > 0:
                    Logger.warning("Large gaps in downloaded data")
                    st.merge()
                    gappy_channels = list(
                        set([(gap[0], gap[1], gap[2], gap[3])
                             for gap in gaps]))
                    _st = Stream()
                    for tr in st:
                        tr_stats = (tr.stats.network, tr.stats.station,
                                    tr.stats.location, tr.stats.channel)
                        if tr_stats in gappy_channels:
                            Logger.warning(
                                "Removing gappy channel: {0}".format(tr))
                        else:
                            _st += tr
                    st = _st
                    st.split()
            st.detrend("simple").merge()
            st.trim(starttime=starttime + (i * data_length) - pad,
                    endtime=starttime + ((i + 1) * data_length) + pad)
            for tr in st:
                if not _check_daylong(tr):
                    st.remove(tr)
                    Logger.warning(
                        "{0} contains more zeros than non-zero, "
                        "removed".format(tr.id))
            for tr in st:
                if tr.stats.endtime - tr.stats.starttime < \
                   0.8 * data_length:
                    st.remove(tr)
                    Logger.warning(
                        "{0} is less than 80% of the required length"
                        ", removed".format(tr.id))
            if return_stream:
                stream += st
            try:
                party += self.detect(
                    stream=st, threshold=threshold,
                    threshold_type=threshold_type, trig_int=trig_int,
                    plot=plot, plotdir=plotdir, daylong=daylong,
                    parallel_process=parallel_process, xcorr_func=xcorr_func,
                    concurrency=concurrency, cores=cores,
                    ignore_length=ignore_length,
                    ignore_bad_data=ignore_bad_data, group_size=group_size,
                    overlap=None, full_peaks=full_peaks,
                    process_cores=process_cores, **kwargs)
                if save_progress:
                    party.write("eqcorrscan_temporary_party")
            except Exception as e:
                Logger.critical(
                    'Error, routine incomplete, returning incomplete Party')
                Logger.error('Error: {0}'.format(e))
                if return_stream:
                    return party, stream
                else:
                    return party
        for family in party:
            if family is not None:
                family.detections = family._uniq().detections
        if return_stream:
            return party, stream
        else:
            return party
