# Data is from 30km radius around Augustine

# Import packages and functions that we need
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
save_dir = '/Volumes/NorthStar/examples/'
hypoi_file = main_dir + 'redoubt_20090215_20090316_hypoi.txt'
hypoddpha_file = main_dir + 'redoubt_20090215_20090316_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)

# read hypoddpha file into a python catalog
catalog_raw = read_hypoddpha(hypoi_file, hypoddpha_file)

# # plot catalog
# fig = catalog_raw.plot(projection="local",resolution="l")

# attempt to filter catalog
#catalog_filtered = filter_picks(catalog=catalog, evaluation_mode="manual", top_n_picks=5)

# # extract station list
# netstachan_list = []
# for i in range(len(catalog_raw)):
#     num_picks = len(catalog_raw[i].picks)
#     for j in range(num_picks):
#         network = catalog_raw[i].picks[j].waveform_id.network_code
#         station = catalog_raw[i].picks[j].waveform_id.station_code
#         channel = catalog_raw[i].picks[j].waveform_id.channel_code
#         netstachan = network + ' ' + station + ' ' + channel
#         netstachan_list.append(netstachan)
# unique_netstachan = np.unique(netstachan_list)
# netstachan_file = '/Users/darrentpk/Desktop/avo_data/combined_list.txt'
# output = open(netstachan_file, 'w')
# for k in range(len(unique_netstachan)):
#     output.write(unique_netstachan[k] + '\n')
# output.close()

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
client = Client("IRIS")
tribe = Tribe().construct(
    method="from_client", lowcut=1.0, highcut=10.0, samp_rate=50.0, length=30.0,
    filt_order=4, prepick=10.0, client_id=client, catalog=catalog, data_pad=20,
    process_len=86400, min_snr=5.0, parallel=False) # process_len = 21600
print(tribe)

# print station numbers of each template before filter
print('Before filter:')
for i in range(len(tribe)):
     t = tribe[i]
     num_stations = len({tr.stats.station for tr in t.st})
     print('Template ' + str(i+1) + ': ' + str(num_stations) + ' stations')
# filter away templates with < 3 stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 5]
print('After removing templates with < 5 stations...')
# print station numbers of each template again
for i in range(len(tribe)):
     t = tribe[i]
     num_stations = len({tr.stats.station for tr in t.st})
     print('Template ' + str(i+1) + ': ' + str(num_stations) + ' stations')

# detect events across desired duration, in time steps to cope with memory limitations
scan_time_start = fetch_time_start
scan_time_end = fetch_time_end
num_time_steps = 10
step_duration = (scan_time_end - scan_time_start)/num_time_steps
# Initialize detection counter for outfile numbering
detection_counter = np.zeros(len(tribe))
# Initialize lists for detection and threshold values
detection_list = []
threshold_list = []
# commence loop for time steps
print('Total number of time steps to be scanned:' + str(num_time_steps))
for i in range(num_time_steps):
    # print statement to keep tabs on time step
    print('Scanning time step ' + str(i+1) + '...')
    # define t1 and t2 for current time step
    t1 = scan_time_start + (i*step_duration)
    t2 = scan_time_start + ((i+1)*step_duration)
    # detect events in step duration
    party, st = tribe.client_detect(
        client=client, starttime=t1, endtime=t2, threshold=20.,
        threshold_type="MAD", trig_int=30.0, plot=False, return_stream=True)
    # check if length of tribe matches number of families
    if len(tribe) != len(party.families):
        print('WARNING: Length of tribe does not match number of families. i = ',str(i))
    # commence a nested loop to extract detection value and store templates and detections
    print('Extracting templates and detections...')
    for j in range(len(party.families)):
        # extract family
        family = party[j]
        # plot template and save fig
        outfile = save_dir + 'temp' + str(j+1) + '.png'
        fig = family.template.st.plot(equal_scale=False, size=(800, 600), outfile=outfile)
        # extract streams for detections
        streams = family.extract_streams(stream=st, length=30.0, prepick=10.0)
        # loop through detections related to this template to save figures
        for k in range(len(family.detections)):
            # add one to detection counter
            detection_counter[j] = detection_counter[j] + 1
            # filter stream and save detection figure
            try:
                stream = streams[family.detections[k].id]
                stream.filter('bandpass', freqmin=1, freqmax=10)
                # plot streams and save fig
                outfile = save_dir + 'temp' + str(j+1) + '_det' + str(detection_counter[j]) + '.png'
                fig = stream.plot(equal_scale=False, size=(800, 600), outfile=outfile)
            except NotImplementedError:
                continue
            # append detection & threshold values for every detection in the current family
            try:
                detection_value = abs(family[k].detect_val)
                detection_list.append(detection_value)
                threshold_value = family[k].threshold
                threshold_list.append(threshold_value)
            except:
                continue
    # before terminating the loop, empty party and st to free memory
    party = None
    st = None

# plot distribution of detections as histogram and save
detection_array = np.array(detection_list)
detection_floor = np.floor(min(detection_array))
detection_ceil = np.ceil(max(detection_array))
threshold_array = np.unique(threshold_list)
fig,ax = plt.subplots()
ax.grid(True)
ax.hist(detection_list, bins=np.arange(detection_floor,detection_ceil,0.1),color='teal',edgecolor='black')
for j in range(len(threshold_array)):
    ax.axvline(x=threshold_array[j],color='red',linewidth=2)
ax.set_xlabel('Detection Value');
ax.set_ylabel('Frequency');
ax.set_title('Histogram of Detection Values');
plt.show()
outfile = save_dir + 'det_hist.png'
plt.savefig(outfile)

# # view most productive family
# family = sorted(party.families, key=lambda f: len(f))[-1]
# print(family)
# fig = family.template.st.plot(equal_scale=False, size=(800, 600))

# # get a dictionary of streams for each detection in a Family and look at those
# streams = family.extract_streams(stream=st, length=30.0, prepick=5.0)
# print(family.detections[0])
# stream = streams[family.detections[0].id]
# stream.filter('bandpass',freqmin=1,freqmax=10)
# fig = stream.plot(equal_scale=False, size=(800, 600))

# # # remove response for streams
# sample_st = streams[family.detections[0].id]
# inv = client.get_stations(network = "AV", station = "*", level = 'response')
# sample_st.remove_response(inventory=inv, output="VEL")
# sample_st.filter('highpass',freq=1)
# fig = sample_st.plot(equal_scale=False, size=(800, 600))

# merge streams and create repicked catalog
st_merged = st.merge(method=1)
repicked_catalog = party.lag_calc(st_merged, pre_processed=False, shift_len=0.5, min_cc=0.4)

# # print repicked catalog events and their pick cross correlations
# for a in range(len(repicked_catalog)-1):
#     print('Event ' + str(a))
#     for b in range(len(repicked_catalog[a].picks)):
#         print(repicked_catalog[a].picks[b].comments)

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
catalog_original