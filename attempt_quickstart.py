# Data is from 30km radius around Augustine

# # Read catalog file
#
# import pandas as pd
#
# catalog_dir = '/Users/darrentpk/Desktop/Research/PREEVENTS/Data/avo_internal/detailed_readable_text.txt'
# catalog = pd.read_csv(catalog_dir,header=0,skiprows=list([1]),skipfooter=1,delim_whitespace=True,engine='python')
# catalog

# # Parse phase info
#
# from phase_processing.hypoi_parser import printer
# from phase_processing.hypoi_parser import construct_phase_catalog
#
# phase_dir = '/Users/darrentpk/Desktop/Research/PREEVENTS/Data/avo_internal/hypoinverse.txt'
# all_Event = construct_phase_catalog(phase_dir)
#
# print('Event Sample:\n----------')
# printer(all_Event[0])
#
# print('\nPhase Sample:\n----------')
# printer(all_Event[0].all_Phase[0])

# Convert hypoi phase data to hypoddpha form
from phase_processing.phase_processing.ncsn2pha import ncsn2pha

main_dir = '/Users/darrentpk/Desktop/avo_data/'
hypoi_file = main_dir + 'augustine2_hypoi.txt'
hypoddpha_file = main_dir + 'augustine2_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)

# read hypoddpha file into a python catalog
from phase_processing.phase_processing.ncsn2pha import read_hypoddpha
catalog = read_hypoddpha(hypoi_file, hypoddpha_file)

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
    process_len=3600, min_snr=5.0, parallel=False)
print(tribe)

####################################
# From EQcorrscan tutorial

import logging

logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from eqcorrscan.utils.catalog_utils import filter_picks

client = Client("NCEDC")
t1 = UTCDateTime(2004, 9, 28)
t2 = t1 + 86400
catalog = client.get_events(
    starttime=t1, endtime=t2, minmagnitude=2.5, minlatitude=35.7, maxlatitude=36.1,
    minlongitude=-120.6, maxlongitude=-120.2, includearrivals=True)
fig = catalog.plot(projection="local", resolution="l")
catalog = filter_picks(catalog=catalog, evaluation_mode="manual", top_n_picks=20)

from eqcorrscan import Tribe

tribe = Tribe().construct(
    method="from_client", lowcut=4.0, highcut=15.0, samp_rate=50.0, length=6.0,
    filt_order=4, prepick=0.5, client_id=client, catalog=catalog, data_pad=20.,
    process_len=21600, min_snr=5.0, parallel=False)
print(tribe)