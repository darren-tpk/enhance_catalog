#%% CREATE TRIBE

# This script creates a tribe based on an input catalog, taking in miniseed files from a local directory or pulling in data from IRIS
# The tribe creation process is done event-by-event, to track progress.
# Tribe is then saved in the user-defined output_dir.

# Import all dependencies
import time
import obspy
import numpy as np
from itertools import compress
from obspy import UTCDateTime, Catalog, Stream
from obspy.clients.fdsn import Client
from eqcorrscan.core.match_filter.tribe import Tribe
from toolbox import download_catalog, get_local_stations, prepare_catalog_stream, reader, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/augustine/'
PEC_dir = main_dir + 'data/avo/'
PEC_file = PEC_dir + 'augustine_20050401_20060501.xml'
output_dir = main_dir + 'output/augustine5/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
sitelist_dir = main_dir + 'data/avo/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe2.tgz'
channel_convention = True   # strict compliance for P/S picks on vertical/horizontal components
use_all_PEC_events = False  # If set to True, we will include PEC events that are non-core RedPy detections
max_dt = 4  # maximum time difference between current template catalog and AVO events allowed, in seconds
resampling_frequency = 50   # for resampling traces prior to final merge
tolerance = 4e4             # tolerance for boxcar removal from data (as a factor to median)
lowcut = 1.0                # filter lowcut (Hz)
highcut = 10.0              # filter highcut (Hz)
samp_rate = 50.0            # new sampling rate (Hz)
length = 8.0                # template length (s), Wech et al. (2018) chose 30s
filt_order = 4              # number of corners for filter
prepick = 1.0               # pre-pick time (s), Wech et al. (2018) chose 5s
process_len = 86400         # length to process data in (s)
min_snr = 1                 # minimum SNR, Jeremy's recommendation was 5.0 (same as EQcorrscan tutorial)
local_volcano = 'augustine' # for get_local_stations function, since we only take picks from stations near Augustine
local_radius = 25           # for get_local_stations function; radius around volcano to accept stations
local = True                # if set to True, use data from local machine
client_name = 'IRIS'        # client name for back-up or non-local data query

#%% Define functions

# Nil

#%% Prepare desired catalog to undergo tribe creation

# Read in, and combine, the picked core catalog and unmatched PEC catalog
core_catalog_picked = reader(convert_redpy_output_dir + 'core_catalog_picked.xml')
if use_all_PEC_events:
    unmatched_PEC_events = reader(convert_redpy_output_dir + 'unmatched_PEC_events_core.xml')
else:
    unmatched_PEC_events = reader(convert_redpy_output_dir + 'unmatched_PEC_events_redpy.xml')
catalog = core_catalog_picked + unmatched_PEC_events

# We only include picks from our hand-picked stations
# mammoth_station_list = ['MINS','MDPB','OMMB','MRD','MDC','MCM','MMP','MLC','MCV','MB01','MB02','MB03','MB05','MB06','MB07','MB08','MB09','MB10','MB11','MQ1P','MMS','MCY','MDY','MLI','MGPB','MLK','MEM','MDR','MMLB','MCB','MLAC']
# redoubt_station_list = ['DFR','NCT','RDDF','RDDR','RDE','RDJH','RDN','RDSO','RDT','RDW','RDWB','RED','REF','RSO','RD01','RD02','RD03']
augustine_station_list = ['AU11','AU12','AU13','AU14','AU15','AU20','AUE','AUH','AUI','AUL','AUNW','AUP','AUR','AUS','AUSE','AUW']
# davidof_station_list = ['LSNW','LSSA','LSPA','LSSE']
# greatsitkin_station_list = ['GSCK','GSIG','GSMY','GSSP','GSTR']

# Clean catalog to only include picks from our station list
for i, event in enumerate(catalog):
    for pick in reversed(event.picks):
        if pick.waveform_id.station_code not in augustine_station_list:
            print('Removed ' + pick.waveform_id.station_code + ' pick from combined catalog index ' + str(i))
            event.picks.remove(pick)

# # Clean catalog to only include picks from local stations
# local_stations = get_local_stations(sitelist_dir,local_volcano,local_radius)
# catalog_out = Catalog()
# for event in catalog:
#     catalog_out.append(event.copy())
#     catalog_out[-1].picks = []
#     for pick in event.picks:
#         if pick.waveform_id.station_code in local_stations:
#             catalog_out[-1].picks.append(pick)
# catalog = catalog_out

#%% Loop over catalog events to create templates, populating tribe

# Initialize tribe and tracker for valid events
tribe = Tribe()
valid_event = []
time_start = time.time()

# Create a boolean list to keep track of which events have been attempted
tracker = np.array([False for i in range(len(catalog))])

# Initialize array of all event times
catalog_times = np.array([event.origins[0].time for event in catalog])

# Initialize a catalog to create tribes with the client method if local data doesn't work
client_catalog = Catalog()

# Loop through events in catalog
print('\nCommencing tribe creation...')
for k, event in enumerate(catalog):

    print(k)

    # Check if event has had its template creation attempted. If True, skip to next
    if tracker[k]:
        continue
    # Otherwise check for other events occuring on the same day so that we load the local data at one go
    else:
        event_daystart = UTCDateTime(event.origins[0].time.date)
        event_dayend = event_daystart + 86400
        sub_catalog_index = (np.array(catalog_times) > event_daystart) & (np.array(catalog_times) < event_dayend)
        sub_catalog = Catalog(list(compress(catalog, sub_catalog_index)))
        tracker[np.where(sub_catalog_index==True)] = True

    # Prepare catalog stream, trim to the start and end of the day to enforce daylong processing
    stream = prepare_catalog_stream(data_dir, sub_catalog, resampling_frequency, tolerance)
    stream = stream.trim(starttime=event_daystart, endtime=event_dayend, pad=True)

    # Do stream checks
    for trace in stream:
        # Check if trace is masked with too many zeros (more than half of samples are masked)
        if hasattr(trace.data, 'mask') and (np.sum(trace.data.mask) / len(trace.data.mask)) > 0.5:
            print('%s.%s got removed due to overly large data gaps.' % (trace.stats.station, trace.stats.channel))
            stream.remove(trace)

    # Construct sub tribe out of sub catalog
    try:
        sub_tribe = Tribe().construct(
            method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
            filt_order=filt_order, prepick=prepick, meta_file=sub_catalog, st=stream, process=True,
            process_len=process_len, min_snr=min_snr, parallel=True)
        if sub_tribe is not None or len(sub_tribe) > 0:
            tribe += sub_tribe
    # If tribe fails then we add the current sub catalog to the client catalog to try the client method later
    except:
        client_catalog += sub_catalog

# Conclude process
time_end = time.time()
print('\nTemplate creation attempts with local data complete. Time elapsed: %.2f s' % (time_end-time_start))

# Now attempt to tackle leftovers with client data
print('\nNow attempting leftovers with client data.')
try:
    client = Client(client_name)
    client_tribe = Tribe().construct(
        method="from_client", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
        filt_order=filt_order, prepick=prepick, client_id=client, catalog=client_catalog, process=True,
        process_len=process_len, min_snr=min_snr, parallel=True)
    tribe += client_tribe
    time_end = time.time()
    print('\nTemplate creation attempts with client data complete. Time elapsed: %.2f s' % (time_end-time_start))
except:
    time_end = time.time()
    print('\nTemplate creation attempts with client data failed. Time elapsed: %.2f s' % (time_end-time_start))

print('%d out of %d events in the catalog were converted to templates.' % (len(tribe),len(catalog)))
writer(create_tribe_output_dir + tribe_filename, tribe)
