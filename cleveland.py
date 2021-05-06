from obspy import Catalog, UTCDateTime, Stream
from obspy.clients.fdsn import Client
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from eqcorrscan import Tribe, Party, Family

hypoi_file = 'cleveland_20210301_20210326_hypoi.txt'
hypoddpha_file = 'cleveland_20210301_20210326_hypoddpha.txt'

catalog_dir = '/home/ptan/attempt_eqcorrscan/avo_data/'
hypoi_path = catalog_dir + hypoi_file
hypoddpha_path = catalog_dir + hypoddpha_file
ncsn2pha(hypoi_path, hypoddpha_path, channel_convention=True)
raw_catalog = read_hypoddpha(hypoi_path, hypoddpha_path, channel_convention=True)

catalog = Catalog()
for raw_event in raw_catalog:
    event = raw_event.copy()
    event.picks.clear()
    clev_index = []
    for i, pick in enumerate(raw_event.picks):
        if pick.waveform_id.station_code == 'CLES' or pick.waveform_id.station_code == 'CLCO':
            clev_index.append(i)
    for j in clev_index:
        event.picks.append(raw_event.picks[j])
    catalog = catalog + event

client = Client('IRIS')
tribe = Tribe().construct(
    method="from_client", lowcut=1.0, highcut=10.0, samp_rate=50.0, length=30.0,
    filt_order=4, prepick=5.0, client_id=client, catalog=catalog,
    process_len=21600, min_snr=5.0, parallel=True)

t1 = UTCDateTime(2021,3,10,12,0,0)
t2 = UTCDateTime(2021,3,26,20,0,0)
party, st = tribe.client_detect(
    client=client, starttime=t1, endtime=t2, threshold=20,
    threshold_type="MAD", trig_int=30, plot=False, return_stream=False)

from toolbox import remove_repeats, gen_detect_hist
time_interval = 30  # number of seconds before and after each detection to check for repeats
party_clean = remove_repeats(party,time_interval)
gen_detect_hist(party_clean)

from waveform_collection import gather_waveforms
party_num = 0
detection_index = 0
detect_time = party_clean[party_num][detection_index].detect_time
all_stachans = party_clean[party_num][detection_index].chans
st = Stream()
for station, channel in all_stachans:
    tr = gather_waveforms(source='IRIS', network='AV', station=station,
                          location='*', channel=channel, starttime=detect_time-5,
                          endtime=detect_time+25)
    st = st + tr
st = st.filter('bandpass',freqmin=1.0, freqmax=10.0, corners=4, zerophase=True)
st = st.taper(0.05, type='hann', max_length=(0.75*1024/50))
st.plot()

all_surplus = []
for i in valid_index:
    detection = party_clean[party_num][i]
    surplus = abs(detection.detect_val)-detection.threshold
    all_surplus.append(surplus)
all_surplus = np.array(all_surplus)

for i,detection in enumerate(party_clean[party_num]):
    print(i, detection.detect_time, abs(detection.detect_val) - detection.threshold)

valid_index = []
for i,detection in enumerate(party_clean[party_num]):
    if (abs(detection.detect_val) - detection.threshold) >0.12:
        valid_index.append(i)



