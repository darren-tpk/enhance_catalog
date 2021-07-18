#%% CREATE TRIBE

# This script creates a tribe based on an input catalog, taking in miniseed files from a local directory or pulling in data from IRIS
# The tribe creation process is done event-by-event, to track progress.
# Tribe is then saved in the user-defined output_dir.

# Import all dependencies
import time
from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client
from eqcorrscan.core.match_filter.tribe import Tribe
from toolbox import get_local_stations, prepare_catalog_stream, reader, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/redoubt/'  # redoubt data directory
output_dir ='/home/ptan/enhance_catalog/output/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
sitelist_dir = main_dir + 'data/avo/'
create_tribe_output_dir = main_dir + 'output/create_tribe/'
tribe_filename = 'tribe_len8_snr2.tgz'
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
local_volcano = 'redoubt'  # for get_local_stations function, since we only take picks from stations near Redoubt
local_radius = 25          # for get_local_stations function; radius around volcano to accept stations

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

    # Prepare catalog stream, trim to the start and end of the day to enforce daylong processing
    print('\nNow at event %d out of %d' % (k+1, len(catalog)))
    event = catalog[k:k+1]
    stream = prepare_catalog_stream(data_dir,event,resampling_frequency,tolerance)
    day_start = UTCDateTime(event[0].origins[0].time.date)
    day_end = day_start + 86400
    stream = stream.trim(starttime=day_start, endtime=day_end)

    # Start with an empty template
    template = None

    # Try to construct tribe using local data file
    try:
        template = Tribe().construct(
            method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
            filt_order=filt_order, prepick=prepick, meta_file=event, st=stream, process=True,
            process_len=process_len, min_snr=min_snr, parallel=True)
    except AttributeError:
        print('WARNING: local data failed to produce template, using client method instead.')

    # If local data files fail (e.g. too gappy), we construct the tribe using IRIS client downloads
    try:
        client = Client('IRIS')
        template = Tribe().construct(
            method="from_client", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
            filt_order=filt_order, prepick=prepick, client_id=client, catalog=event, process=True,
            process_len=process_len, min_snr=min_snr, parallel=True)
    except Exception:
        print('WARNING: data not available on client either, skipping.')

    # Append template to tribe if successful
    if template is not None:
        tribe = tribe + template
        time_current = time.time()
        print('Event %d passed.' % (k+1))
        valid_event.append(1)

    # Skip if unsuccessful
    else:
        print('Event %d failed.' % (k+1))
        valid_event.append(0)

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
