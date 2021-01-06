# Data is from 25km radius around Augustine

# Convert hypoi phase data to hypoddpha form
from phase_processing.phase_processing.ncsn2pha import ncsn2pha

main_dir = '/Users/darrentpk/Desktop/avo_data/'
hypoi_file = main_dir + 'augustine2015_2017_hypoi.txt'
hypoddpha_file = main_dir + 'augustine2015_2017_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)

# read hypoddpha file into a python catalog
from phase_processing.phase_processing.read_hypoddpha import read_hypoddpha
catalog = read_hypoddpha(hypoi_file, hypoddpha_file)

# # plot catalog
# fig = catalog.plot(projection="local",resolution="l")

# attempt to filter catalog
#from eqcorrscan.utils.catalog_utils import filter_picks
#catalog_filtered = filter_picks(catalog=catalog, evaluation_mode="manual", top_n_picks=5)

# # extract magnitudes
# import numpy as np
# mag_list = []
# for i in range(len(catalog)):
#     mag_list.append(catalog[i].magnitudes[0].mag)
# all_mag = np.array(mag_list)

# sub-sample catalog, picking only events with mag >1
from obspy import Catalog
catalog_sample = Catalog()
for i in range(len(catalog)):
    mag = catalog[i].magnitudes[0].mag
    if mag > 1:
        catalog_sample.append(catalog[i])
catalog = catalog_sample

# construct tribe
from obspy.clients.fdsn import Client
client = Client("IRIS")
from eqcorrscan import Tribe
tribe = Tribe().construct(
    method="from_client", lowcut=4.0, highcut=15.0, samp_rate=50.0, length=6.0,
    filt_order=4, prepick=0.5, client_id=client, catalog=catalog, data_pad=20.,
    process_len=86400, min_snr=5.0, parallel=True) # process_len = 21600
print(tribe)

# # to look at one template
# print(tribe[2])
# fig = tribe[2].st.plot(equal_scale=False, size=(800, 600))

# detect events
from obspy import UTCDateTime
t1 = UTCDateTime(2016,8,26,0,0,0)
# only use tribes that have 5 or more stations...
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 5] # 5
print(tribe)
party, st = tribe.client_detect(
    client=client, starttime=t1, endtime=t1 + (86400 * 2), threshold=9.,
    threshold_type="MAD", trig_int=2.0, plot=False, return_stream=True)

# view most productive family
family = sorted(party.families, key=lambda f: len(f))[-1]
print(family)
fig = family.template.st.plot(equal_scale=False, size=(800, 600))

# get a dictionary of streams for each detection in a Family and look at those
streams = family.extract_streams(stream=st, length=10, prepick=2.5)
print(family.detections[0])
fig = streams[family.detections[0].id].plot(equal_scale=False, size=(800, 600))

# # remove response for streams
sample_st = streams[family.detections[0].id]
inv = client.get_stations(network = "AV", station = "*", level = 'response')
sample_st.remove_response(inventory=inv, output="VEL")
sample_st.filter('highpass',freq=1)
fig = sample_st.plot(equal_scale=False, size=(800, 600))

st_merged = st.merge(method=1)
repicked_catalog = party.lag_calc(st_merged, pre_processed=False, shift_len=0.5, min_cc=0.4)
print(repicked_catalog[0].picks[0])
print(repicked_catalog[0].picks[0].comments[0])
