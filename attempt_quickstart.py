# Data is from 25km radius around Augustine

# Import packages and functions that we need
import numpy as np
from obspy import Catalog
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from eqcorrscan import Tribe
from eqcorrscan.utils.catalog_utils import filter_picks
from phase_processing.phase_processing.ncsn2pha import ncsn2pha
from phase_processing.phase_processing.read_hypoddpha import read_hypoddpha

# Convert hypoi phase data to hypoddpha form
main_dir = '/Users/darrentpk/Desktop/avo_data/'
hypoi_file = main_dir + 'augustine2015_2017_hypoi.txt'
hypoddpha_file = main_dir + 'augustine2015_2017_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)

# read hypoddpha file into a python catalog
catalog_raw = read_hypoddpha(hypoi_file, hypoddpha_file)

# # plot catalog
# fig = catalog.plot(projection="local",resolution="l")

# attempt to filter catalog
#catalog_filtered = filter_picks(catalog=catalog, evaluation_mode="manual", top_n_picks=5)

# sub-sample catalog, picking only events with mag >1
catalog_sample = Catalog()
for i in range(len(catalog_raw)):
    mag = catalog_raw[i].magnitudes[0].mag
    if mag > 1:
        catalog_sample.append(catalog_raw[i])
catalog = catalog_sample

# construct tribe
client = Client("IRIS")
tribe = Tribe().construct(
    method="from_client", lowcut=4.0, highcut=15.0, samp_rate=50.0, length=6.0,
    filt_order=4, prepick=0.5, client_id=client, catalog=catalog, data_pad=20.,
    process_len=86400, min_snr=5.0, parallel=True) # process_len = 21600
print(tribe)

# detect events
t1 = UTCDateTime(2016,8,26,0,0,0)
t2 = t1 + (86400 * 3)
# only use tribes that have 5 or more stations...
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 5] # 5
print(tribe)
party, st = tribe.client_detect(
    client=client, starttime=t1, endtime=t2, threshold=9.,
    threshold_type="MAD", trig_int=2.0, plot=False, return_stream=True)

# # view most productive family
# family = sorted(party.families, key=lambda f: len(f))[-1]
# print(family)
# fig = family.template.st.plot(equal_scale=False, size=(800, 600))

# # get a dictionary of streams for each detection in a Family and look at those
# streams = family.extract_streams(stream=st, length=10, prepick=2.5)
# print(family.detections[0])
# #fig = streams[family.detections[0].id].plot(equal_scale=False, size=(800, 600))
#
# # # remove response for streams
# sample_st = streams[family.detections[0].id]
# inv = client.get_stations(network = "AV", station = "*", level = 'response')
# sample_st.remove_response(inventory=inv, output="VEL")
# sample_st.filter('highpass',freq=1)
# fig = sample_st.plot(equal_scale=False, size=(800, 600))

inv = client.get_stations(network = "AV", station = "*", level = 'response')
# commence a loop to store figures of each template and their detections
for j in range(len(party.families)):
    # extract family
    family = party[j]
    # plot template and save fig
    outfile = '/Users/darrentpk/Desktop/examples/temp' + str(j) + '.png'
    fig = family.template.st.plot(equal_scale=False, size=(800, 600), outfile=outfile)
    # loop through detections related to this template
    for k in range(len(family.detections)):
        # remove response and filter
        streams = family.extract_streams(stream=st, length=10, prepick=2.5)
        stream = streams[family.detections[k].id]
        stream.remove_response(inventory=inv, output="VEL")
        stream.filter('highpass',freq=1)
        # plot streams and save fig
        outfile = '/Users/darrentpk/Desktop/examples/temp' + str(j) + '_det' + str(k) + '.png'
        fig = stream.plot(equal_scale=False, size=(800, 600), outfile=outfile)

# merge streams and create repicked catalog
st_merged = st.merge(method=1)
repicked_catalog = party.lag_calc(st_merged, pre_processed=False, shift_len=0.5, min_cc=0.4)

# print repicked catalog events and their pick cross correlations
for a in range(len(repicked_catalog)-1):
    print('Event ' + str(a))
    for b in range(len(repicked_catalog[a].picks)):
        print(repicked_catalog[a].picks[b].comments)

# print UTCDatetime of repicked catalog events
all_time = []
for i in range(len(repicked_catalog)):
    time_str = str(repicked_catalog[i].resource_id)[20:]
    time = UTCDateTime(int(time_str[0:4]),int(time_str[4:6]),int(time_str[6:8]),int(time_str[9:11]),int(time_str[11:13]),int(time_str[13:15]),int(time_str[15:]))
    all_time.append(time)
    print('Event ' + str(i+1) + ':  ' + str(time))
    for j in range(len(repicked_catalog[i].picks)):
        print(repicked_catalog[i].picks[j].comments)
all_time = list(np.sort(all_time))

# cross-check with original catalog
catalog_sample = Catalog()
for i in range(len(catalog_raw)):
    time = catalog_raw[i].origins[0].time
    if time > t1 and time < t2:
        catalog_sample.append(catalog_raw[i])
catalog_original = catalog_sample