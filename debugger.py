import numpy as np
import logging
from obspy.clients.fdsn import Client
from obspy import Catalog
from obspy import Stream
from obspy import read_events
from obspy import UTCDateTime
from eqcorrscan.core import template_gen
from eqcorrscan.core.template_gen import _group_events
from eqcorrscan.core.template_gen import _download_from_client

class TemplateGenError(Exception):
    def __init__(self, value):
        self.value = value
    def __repr__(self):
        return self.value
    def __str__(self):
        return 'TemplateGenError: ' + self.value

from waveform_collection import gather_waveforms
st = gather_waveforms(source='IRIS', network='AV', station='AUL',
                           location='', channel='*N', starttime=UTCDateTime(2015, 1, 15, 5, 27, 45, 560000),
                           endtime=UTCDateTime(2015, 1, 16, 5, 27, 45, 560000))

swin = 'all'
lowcut = 4
highcut = 15
samp_rate = 50
length = 6
filt_order = 4
prepick = 0.5
client_id = Client('IRIS')
catalog = catalog
data_pad = 20
process_len = 86400
min_snr = None
all_horiz = False
delayed = True
plot = False
plotdir = None
parallel = False
num_cores = False
skip_short_chans = False
save_progress = False
url = None

# Sub-sample our catalog
catalog_raw = read_events(output_file, "HYPODDPHA")
catalog = Catalog()
for i in range(0,10):
    num_picks = len(catalog_raw[i].picks)
    for j in range(num_picks):
        catalog_raw[i].picks[j].waveform_id.network_code = catalog_raw[i].picks[j].waveform_id.station_code[0:2]
        catalog_raw[i].picks[j].waveform_id.station_code = catalog_raw[i].picks[j].waveform_id.station_code[2:]
        catalog_raw[i].picks[j].waveform_id.channel_code  = "*" + catalog_raw[i].picks[j].waveform_id.channel_code
    if num_picks != 0:
        catalog.append(catalog_raw[i])

# now try template_gen
client_map = {'from_client': 'fdsn', 'from_seishub': 'seishub'}
assert method in ('from_client', 'from_seishub', 'from_meta_file',
                  'from_sac')
if not isinstance(swin, list):
    swin = [swin]
process = True
if method in ['from_client', 'from_seishub']:
    catalog = catalog
    data_pad = data_pad
    # Group catalog into days and only download the data once per day
    sub_catalogs = _group_events(
        catalog=catalog, process_len=process_len, template_length=length,
        data_pad=data_pad)
    if method == 'from_client':
        if isinstance(client_id, str):
            client = FDSNClient(client_id)
        else:
            client = client_id
        available_stations = []
    else:
        client = SeisHubClient(url, timeout=10)
        available_stations = client.waveform.get_station_ids()
elif method == 'from_meta_file':
    if isinstance(meta_file, Catalog):
        catalog = meta_file
    elif meta_file:
        catalog = read_events(meta_file)
    else:
        catalog = catalog
    sub_catalogs = [catalog]
    st = Stream()
    process = process
elif method == 'from_sac':
    sac_files = sac_files
    if isinstance(sac_files, list):
        if isinstance(sac_files[0], (Stream, Trace)):
            # This is a list of streams...
            st = Stream(sac_files[0])
            for sac_file in sac_files[1:]:
                st += sac_file
        else:
            sac_files = [read(sac_file)[0] for sac_file in sac_files]
            st = Stream(sac_files)
    else:
        st = sac_files
    # Make an event object...
    catalog = Catalog([sactoevent(st)])
    sub_catalogs = [catalog]

temp_list = []
process_lengths = []
catalog_out = Catalog()

if "P_all" in swin or "S_all" in swin or all_horiz:
    all_channels = True
else:
    all_channels = False
for sub_catalog in sub_catalogs:
    if method in ['from_seishub', 'from_client']:
        Logger.info("Downloading data")
        client_type = client_map[method]
        st = Stream()
        sub_catalog = Catalog(sorted(sub_catalog, key=lambda e: e.origins[0].time))
        all_waveform_info = []
        for event in sub_catalog:
            for pick in event.picks:
                if not pick.waveform_id:
                    Logger.warning(
                        "Pick not associated with waveforms, will not use:"
                        " {0}".format(pick))
                    continue
                if all_channels:
                    channel_code = pick.waveform_id.channel_code[0:2] + "?"
                else:
                    channel_code = pick.waveform_id.channel_code
                all_waveform_info.append((
                    pick.waveform_id.network_code, pick.waveform_id.station_code,
                    channel_code, pick.waveform_id.location_code))
        starttime = UTCDateTime(
            sub_catalog[0].origins[0].time - data_pad)
        endtime = starttime + process_len
        # Check that endtime is after the last event
        if not endtime > sub_catalog[-1].origins[0].time + data_pad:
            raise TemplateGenError(
                'Events do not fit in processing window')
        all_waveform_info = sorted(list(set(all_waveform_info)))
        dropped_pick_stations = 0
        for waveform_info in all_waveform_info:
            net, sta, chan, loc = waveform_info
            if client_type == 'seishub' and sta not in available_stations:
                Logger.error("Station not found in SeisHub DB")
                dropped_pick_stations += 1
                continue
            Logger.info('Downloading for start-time: {0} end-time: {1}'.format(
                starttime, endtime))
            Logger.debug('.'.join([net, sta, loc, chan]))
            query_params = dict(
                network=net, station=sta, location=loc, channel=chan,
                starttime=starttime, endtime=endtime)
            try:
                st += client.get_waveforms(**query_params)
            except Exception as e:
                Logger.error(e)
                Logger.error('Found no data for this station: {0}'.format(
                    query_params))
                dropped_pick_stations += 1
        if not st and dropped_pick_stations == len(event.picks):
            raise Exception('No data available, is the server down?')
        st.merge()
        # clients download chunks, we need to check that the data are
        # the desired length
        final_channels = []
        for tr in st:
            tr.trim(starttime, endtime)
            if len(tr.data) == (process_len * tr.stats.sampling_rate) + 1:
                tr.data = tr.data[1:len(tr.data)]
            if tr.stats.endtime - tr.stats.starttime < 0.8 * process_len:
                Logger.warning(
                    "Data for {0}.{1} is {2} hours long, which is less than 80 "
                    "percent of the desired length, will not use".format(
                        tr.stats.station, tr.stats.channel,
                        (tr.stats.endtime - tr.stats.starttime) / 3600))
            elif not pre_processing._check_daylong(tr):
                Logger.warning(
                    "Data are mostly zeros, removing trace: {0}".format(tr.id))
            else:
                final_channels.append(tr)
        st.traces = final_channels
    Logger.info('Pre-processing data')
    st.merge()
    if len(st) == 0:
        Logger.info("No data")
        continue
    if process:
        data_len = max([len(tr.data) / tr.stats.sampling_rate
                        for tr in st])
        if 80000 < data_len < 90000:
            daylong = True
            starttime = min([tr.stats.starttime for tr in st])
            min_delta = min([tr.stats.delta for tr in st])
            # Cope with the common starttime less than 1 sample before the
            #  start of day.
            if (starttime + min_delta).date > starttime.date:
                starttime = (starttime + min_delta)
            # Check if this is stupid:
            if abs(starttime - UTCDateTime(starttime.date)) > 600:
                daylong = False
            starttime = starttime.date
        else:
            daylong = False
        # Check if the required amount of data have been downloaded - skip
        # channels if arg set.
        for tr in st:
            if np.ma.is_masked(tr.data):
                _len = np.ma.count(tr.data) * tr.stats.delta
            else:
                _len = tr.stats.npts * tr.stats.delta
            if _len < process_len * .8:
                Logger.info(
                    "Data for {0} are too short, skipping".format(
                        tr.id))
                if skip_short_chans:
                    continue
            # Trim to enforce process-len
            tr.data = tr.data[0:int(process_len * tr.stats.sampling_rate)]
        if len(st) == 0:
            Logger.info("No data")
            continue
        if daylong:
            st = pre_processing.dayproc(
                st=st, lowcut=lowcut, highcut=highcut,
                filt_order=filt_order, samp_rate=samp_rate,
                parallel=parallel, starttime=UTCDateTime(starttime),
                num_cores=num_cores)
        else:
            st = pre_processing.shortproc(
                st=st, lowcut=lowcut, highcut=highcut,
                filt_order=filt_order, parallel=parallel,
                samp_rate=samp_rate, num_cores=num_cores)
    data_start = min([tr.stats.starttime for tr in st])
    data_end = max([tr.stats.endtime for tr in st])

    for event in sub_catalog:
        stations, channels, st_stachans = ([], [], [])
        if len(event.picks) == 0:
            Logger.warning(
                'No picks for event {0}'.format(event.resource_id))
            continue
        use_event = True
        # Check that the event is within the data
        for pick in event.picks:
            if not data_start < pick.time < data_end:
                Logger.warning(
                    "Pick outside of data span: Pick time {0} Start "
                    "time {1} End time: {2}".format(
                        str(pick.time), str(data_start), str(data_end)))
                use_event = False
        if not use_event:
            Logger.error('Event is not within data time-span')
            continue
        # Read in pick info
        Logger.debug("I have found the following picks")
        for pick in event.picks:
            if not pick.waveform_id:
                Logger.warning(
                    'Pick not associated with waveforms, will not use:'
                    ' {0}'.format(pick))
                continue
            Logger.debug(pick)
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
        # Check to see if all picks have a corresponding waveform
        for tr in st:
            st_stachans.append('.'.join([tr.stats.station,
                                         tr.stats.channel]))
        # Cut and extract the templates
        template = _template_gen(
            event.picks, st, length, swin, prepick=prepick, plot=plot,
            all_horiz=all_horiz, delayed=delayed, min_snr=min_snr,
            plotdir=plotdir)
        process_lengths.append(len(st[0].data) / samp_rate)
        temp_list.append(template)
        catalog_out += event
    if save_progress:
        if not os.path.isdir("eqcorrscan_temporary_templates"):
            os.makedirs("eqcorrscan_temporary_templates")
        for template in temp_list:
            template.write(
                "eqcorrscan_temporary_templates{0}{1}.ms".format(
                    os.path.sep, template[0].stats.starttime.strftime(
                        "%Y-%m-%dT%H%M%S")),
                format="MSEED")
    del st