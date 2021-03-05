# create tribes from input catalog and miniseed files

# import packages that we need
import os
import glob
from obspy import UTCDateTime, Stream, Catalog, read
from eqcorrscan.core.match_filter.tribe import Tribe
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
# import toolbox functions
from toolbox import remove_boxcars

# define all variables here
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
start_time = UTCDateTime(2009, 2, 26, 0, 0, 0)
end_time = start_time + (24 * 60 * 60)
tolerance = 5e4     # for boxcar removal
lowcut = 1.0        # filter lowcut (Hz), follows Jeremy's recommendation
highcut = 10.0      # filter highcut (Hz), follows Jeremy's recommendation
samp_rate = 50.0    # new sampling rate (Hz)
length = 30.0       # template length (s), follows Wech et al. (2018)
filt_order = 4      # number of corners for filter
prepick = 5.0       # pre-pick time (s), follows Wech et al. (2018)
process_len = 86400 # length to process (s)
min_snr = 5.0       # minimum SNR, follows Jeremy's recommendation
tribe_out = 'tribe' # name of tribe tar file saved in output

# redoubt data directory
data_dir = '/home/data/redoubt/'
data_filenames_complete = os.listdir(data_dir)

# convert hypoi phase data to hypoddpha form
catalog_dir = '/home/ptan/project/avo_data/'
hypoi_path = catalog_dir + hypoi_file
hypoddpha_path = catalog_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path)

# read hypoddpha file into a python catalog
catalog_raw = read_hypoddpha(hypoi_path, hypoddpha_path)

# sub-sample catalog, picking only events on swarm day for testing
catalog_sample = Catalog()
for i in range(len(catalog_raw)):
    time = catalog_raw[i].origins[0].time
    if time >= start_time and time <= end_time:
        catalog_sample.append(catalog_raw[i])
catalog = catalog_sample

# now prepare stream
# get unique list of all station channel data needed for the day's events
data_filenames = []
for event in catalog:
    for pick in event.picks:
        # extract key information from pick
        sta = pick.waveform_id.station_code
        chan = pick.waveform_id.channel_code
        pick_year = pick.time.year
        pick_julday = pick.time.julday
        # craft file string and append
        data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{pick_julday:03}' + ':*'
        data_filenames.append(data_filename)
        # add next day if pick occurs in the first 15 minutes of the day
        if pick.time.hour == 0 and pick.time.minute < 15:
            if pick.time.julday == 1: # special case for first day of the year
                data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year-1) + ':365:*'
                data_filenames.append(data_filename)
            else:
                data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{(pick_julday-1):03}' + ':*'
                data_filenames.append(data_filename)
        # add previous day if pick occurs in the last 15 minutes of the day
        if pick.time.hour == 23 and pick.time.minute > 45:
            if pick.time.julday == 365: # special case for last day of the year
                data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year+1) + ':001:*'
                data_filenames.append(data_filename)
            else:
                data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{(pick_julday+1):03}' + ':*'
                data_filenames.append(data_filename)
# now compile unique and sort
data_filenames = list(set(data_filenames))
data_filenames.sort()

# read in all streams we need
stream = Stream()
for data_filename in data_filenames:
    real_data_filename = (glob.glob(data_filename)).pop()
    stream_contribution = read(real_data_filename)
    stream = stream + stream_contribution
stream = remove_boxcars(stream,tolerance)
stream.detrend("simple").merge()

# now build templates
print('Constructing templates ...')
run_time = UTCDateTime()
tribe = Tribe().construct(
    method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
    filt_order=filt_order, prepick=prepick, meta_file=catalog, st=stream, process=True,
    process_len=process_len, min_snr=min_snr, parallel=True)
print('Time Elapsed: ',UTCDateTime()-run_time)

# write tribe into output
tribe_outpath = '/home/ptan/project/output/' + tribe_out
if os.path.exists(tribe_outpath+'.tgz'):
    os.remove(tribe_outpath+'.tgz')
    tribe.write(tribe_outpath)
else:
    tribe.write(tribe_outpath)