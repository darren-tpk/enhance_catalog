# create tribes from input catalog and miniseed files

# import packages that we need
import os
from obspy import UTCDateTime, Catalog
from eqcorrscan.core.match_filter.tribe import Tribe
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
# import toolbox functions
from toolbox import prepare_catalog_stream, writer

# define all variables here
data_dir = '/home/data/redoubt/'  # redoubt data directory
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
channel_convention = True  # strict compliance for P/S picks on vertical/horizontal components
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

# convert hypoi phase data to hypoddpha form
catalog_dir = '/home/ptan/attempt_eqcorrscan/avo_data/'
hypoi_path = catalog_dir + hypoi_file
hypoddpha_path = catalog_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path, channel_convention)

# read hypoddpha file into a python catalog
catalog_raw = read_hypoddpha(hypoi_path, hypoddpha_path, channel_convention)

# sub-sample catalog, picking only events on swarm day for testing
catalog_sample = Catalog()
for i in range(len(catalog_raw)):
    time = catalog_raw[i].origins[0].time
    if time >= start_time and time <= end_time:
        catalog_sample.append(catalog_raw[i])
catalog = catalog_sample

stream = prepare_catalog_stream(data_dir,catalog,tolerance)

# now build templates
print('Constructing templates ...')
run_time = UTCDateTime()
tribe = Tribe().construct(
    method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, length=length,
    filt_order=filt_order, prepick=prepick, meta_file=catalog, st=stream, process=True,
    process_len=process_len, min_snr=min_snr, parallel=True)
print('Time Elapsed: ',UTCDateTime()-run_time)

# write tribe into output
tribe_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(tribe_outpath+'tribe.tgz', tribe)