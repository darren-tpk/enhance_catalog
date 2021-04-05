# create tribes from input catalog and miniseed files

# import packages that we need
import time
from obspy import UTCDateTime, Catalog
from eqcorrscan.core.match_filter.tribe import Tribe
# import toolbox functions
from toolbox import get_local_stations, prepare_catalog_stream, reader, writer

# define all variables here
data_dir = '/home/data/redoubt/'  # redoubt data directory
channel_convention = True  # strict compliance for P/S picks on vertical/horizontal components
start_time = UTCDateTime(2009, 2, 26, 0, 0, 0)
end_time = start_time + (24 * 60 * 60)
resampling_frequency = 50  # for resampling traces prior to final merge
tolerance = 4e4            # for boxcar removal
lowcut = 1.0               # filter lowcut (Hz), follows Jeremy's recommendation
highcut = 10.0             # filter highcut (Hz), follows Jeremy's recommendation
samp_rate = 50.0           # new sampling rate (Hz)
length = 20.0              # template length (s), follows Wech et al. (2018)
filt_order = 4             # number of corners for filter
prepick = 5.0              # pre-pick time (s), follows Wech et al. (2018)
process_len = 86400        # length to process (s)
min_snr = 5.0              # minimum SNR, follows Jeremy's recommendation
local_volcano = 'redoubt'  # for get_local_stations function, since we only take picks from stations near Redoubt
local_radius = 25          # for get_local_stations function; radius around volcano to accept stations

# convert hypoi phase data to hypoddpha form
catalog_inpath = '/home/ptan/attempt_eqcorrscan/output/'
catalog = reader(catalog_inpath+'core_catalog_picked.xml')

# clean catalog to only include picks from local stations
local_stations = get_local_stations(local_volcano,local_radius)
catalog_out = Catalog()
for event in catalog:
    catalog_out.append(event.copy())
    catalog_out[-1].picks = []
    for pick in event.picks:
        if pick.waveform_id.station_code in local_stations:
            catalog_out[-1].picks.append(pick)
catalog = catalog_out

# create tribes event-by-event
# initialize tribe and tracker for valid events
tribe = Tribe()
valid_event = []
time_start = time.time()
# loop through events in catalog
for k in range(len(catalog)):
    # prepare catalog stream, trim to the start and end of the day to enforce daylong processing
    print('\nNow at k = %d' % k)
    event = catalog[k:k+1]
    stream = prepare_catalog_stream(data_dir,event,resampling_frequency,tolerance)
    day_start = UTCDateTime(event[0].origins[0].time.date)
    day_end = day_start + 86400
    stream.trim(starttime=day_start, endtime=day_end)
    # construct template
    try:
        # construct using local data file
        template = Tribe().construct(
            method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
            filt_order=filt_order, prepick=prepick, meta_file=event, st=stream, process=True,
            process_len=process_len, min_snr=min_snr, parallel=True)
    except:
        # if local data files fail, use client download
        print('WARNING: local data failed to produce template, using client method instead.')
        from obspy.clients.fdsn import Client
        client = Client('IRIS')
        template = Tribe().construct(
            method="from_client", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
            filt_order=filt_order, prepick=prepick, client_id=client, catalog=event, process=True,
            process_len=process_len, min_snr=min_snr, parallel=True)
    # append if template is created
    if len(template) != 0:
        tribe = tribe + template
        time_current = time.time()
        print('k = %d passed.' % k)
        valid_event.append(1)
    else:
        print('k = %d failed.' % k)
        valid_event.append(0)
time_end = time.time()
print('\nTemplate creation complete. Time taken: %.2f s' % (time_end-time_start))

# bulk EQcorrscan (problematic as of Mar 29)
# stream = prepare_catalog_stream(data_dir,catalog,resampling_frequency,tolerance)
# # now build templates
# print('Constructing templates ...')
# run_time = UTCDateTime()
# tribe = Tribe().construct(
#     method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
#     filt_order=filt_order, prepick=prepick, meta_file=catalog, st=stream, process=True,
#     process_len=process_len, min_snr=min_snr, parallel=True)
# print('Time Elapsed: ',UTCDateTime()-run_time)

# write tribe into output
tribe_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(tribe_outpath+'tribe.tgz', tribe)