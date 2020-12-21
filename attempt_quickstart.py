# Data is from 30km radius around Augustine

# Read catalog file

import pandas as pd

catalog_dir = '/Users/darrentpk/Desktop/Research/PREEVENTS/Data/avo_internal/detailed_readable_text.txt'
catalog = pd.read_csv(catalog_dir,header=0,skiprows=list([1]),skipfooter=1,delim_whitespace=True,engine='python')
catalog

# Parse phase info

from prepare_catalog.parser import printer
from prepare_catalog.parser import construct_phase_catalog

phase_dir = '/Users/darrentpk/Desktop/Research/PREEVENTS/Data/avo_internal/hypoinverse.txt'
all_Event = construct_phase_catalog(phase_dir)

print('Event Sample:\n----------')
printer(all_Event[0])

print('\nPhase Sample:\n----------')
printer(all_Event[0].all_Phase[0])

# Convert hypoi phase data to hypoddpha form
from prepare_catalog.ncsn2pha import ncsn2pha

main_dir = '/Users/darrentpk/opt/anaconda3/envs/eqcorrscan/darren/avo_data/'
input_file = main_dir + 'augustine2_hypoi.txt'
output_file = main_dir + 'augustine2_hypoddpha.txt'
ncsn2pha(input_file, output_file)

# read hypoddpha file into a python catalog
from obspy import read_events
catalog = read_events(output_file, "HYPODDPHA")

# quick loop to fix network and station entries, and remove EQs without picks
has_picks = []
for i in range(len(catalog)):
    num_picks = len(catalog[i].picks)
    if num_picks != 0:
        has_picks.append(i)
    for j in range(num_picks):
        catalog[i].picks[j].waveform_id.network_code = catalog[i].picks[j].waveform_id.station_code[0:2]
        catalog[i].picks[j].waveform_id.station_code = catalog[i].picks[j].waveform_id.station_code[2:]

# plot catalog
fig = catalog.plot(projection="local",resolution="l")

# attempt to filter catalog
from eqcorrscan.utils.catalog_utils import filter_picks
catalog_filtered = filter_picks(catalog=catalog, evaluation_mode="manual", top_n_picks=5)

# construct tribe
from obspy.clients.fdsn import Client
client = Client("IRIS")
from eqcorrscan import Tribe
tribe = Tribe().construct(
    method="from_client", lowcut=4.0, highcut=15.0, samp_rate=50.0, length=6.0,
    filt_order=4, prepick=0.5, client_id=client, catalog=catalog, data_pad=20.,
    process_len=21600, min_snr=5.0, parallel=True)
print(tribe)
