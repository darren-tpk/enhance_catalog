# Create tribes from input catalog

# Import packages and functions that we need
import pickle
import numpy as np
import matplotlib.pyplot as plt
from obspy import Catalog
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from eqcorrscan import Tribe
from eqcorrscan.utils.catalog_utils import filter_picks
from phase_processing.phase_processing.ncsn2pha import ncsn2pha
from phase_processing.phase_processing.read_hypoddpha import read_hypoddpha

# Convert hypoi phase data to hypoddpha form
main_dir = '/Users/darrentpk/Desktop/avo_data/'
hypoi_file = main_dir + 'redoubt_20090215_20090316_hypoi.txt'
hypoddpha_file = main_dir + 'redoubt_20090215_20090316_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)

# read hypoddpha file into a python catalog
catalog_raw = read_hypoddpha(hypoi_file, hypoddpha_file)

# sub-sample catalog, picking only events with mag >0.9
catalog_sample = Catalog()
fetch_time_start = UTCDateTime(2009,2,15,21,0,0)
fetch_time_end = UTCDateTime(2009,3,15,21,0,0)
for i in range(len(catalog_raw)):
    mag = catalog_raw[i].magnitudes[0].mag
    time = catalog_raw[i].origins[0].time
    # if mag >1:
    #     catalog_sample.append(catalog_raw[i])
    if time >= fetch_time_start and time <= fetch_time_end and mag >0.9:
        catalog_sample.append(catalog_raw[i])
catalog = catalog_sample

# construct tribe
run_time = UTCDateTime()
client = Client("IRIS")
tribe = Tribe().construct(
    method="from_client", lowcut=1.0, highcut=10.0, samp_rate=50.0, length=30.0,
    filt_order=4, prepick=10.0, client_id=client, catalog=catalog, data_pad=20,
    process_len=86400, min_snr=5.0, parallel=False) # process_len = 21600
print(tribe)
print('Time Elapsed: ',UTCDateTime()-run_time)

# save tribe
save_dir = '/Volumes/NorthStar/data/'
tribe_file = save_dir + 'tribe.pickle'
tribe_out = open(tribe_file,'wb')
pickle.dump(tribe,tribe_out)
tribe_out.close()
