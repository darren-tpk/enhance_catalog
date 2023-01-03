
from toolbox import reader
from obspy import UTCDateTime
import numpy as np
from obspy import Catalog
import matplotlib.pyplot as plt
# file:///Users/darrentpk/Desktop/GitHub/enhance_catalog/redpy/runs/redoubt2/clusters/1022.html

main_dir = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/'
PEC = reader(main_dir + 'output/redoubt5/scan_data/PEC_FI.xml')
tribe = reader(main_dir + 'output/redoubt5/create_tribe/tribe2.tgz')
catalog = reader(main_dir + 'output/redoubt5/scan_data/party_FI.xml')
cat2 = reader(main_dir + 'output/augustine5/redoubt5/relocatable_catalog_FImag.xml')

all_pec_times = [e.origins[0].time for e in PEC]
all_template_times = [t.event.origins[0].time for t in tribe]
all_detection_times = [UTCDateTime(e.resource_id.id.split('_')[-1]) for e in catalog]

valid_templates = []
for t in tribe.templates:
    if len(t.event.origins[0].comments) != 0:
        txt_split = t.event.origins[0].comments[0].text.split(' ')
        if txt_split[2] == 'CORE' and txt_split[4] != 'NOT' and len(t.event.picks) >5:
            valid_templates.append(t.name)

all_party_temp = ['_'.join(e.resource_id.id.split('/')[1].split('_')[:-1]) for e in catalog]
all_det_count = [all_party_temp.count(vt) for vt in valid_templates]
valid_det_count = np.flatnonzero([dc > 15 for dc in all_det_count])
valid_templates2 = np.array(valid_templates)[valid_det_count]


core_time = UTCDateTime(2006,1,12,6,3,44)
# '2006_01_12t06_03_44'
# '2006_01_14t01_27_05'
# '2006_01_27t00_59_15'
# '2006_01_27t22_54_16'
all_t_diff = [t - core_time for t in all_pec_times]
abs_t_diff = [np.abs(t) for t in all_t_diff]
PEC_i = np.argmin(abs_t_diff)
all_t_diff = [t - core_time for t in all_template_times]
abs_t_diff = [np.abs(t) for t in all_t_diff]
core_i = np.argmin(abs_t_diff)

pec_event = PEC[PEC_i]
core_template = tribe.templates[core_i]

sub_catalog = Catalog()
for event in catalog:
    parent_template_name = event.resource_id.id.split('/')[1][0:19]
    if parent_template_name == core_template.name:
        sub_catalog.events.append(event)
select_detection_times = [UTCDateTime(e.resource_id.id.split('_')[-1]) for e in sub_catalog]
# select_detection_times = list(np.sort(select_detection_times))

# Define client and target network
from obspy.clients.fdsn import Client
from obspy import Stream
client = Client('IRIS')
stream = Stream()
for pick in pec_event.picks:
    net = pick.waveform_id.network_code
    sta = pick.waveform_id.station_code
    chan = pick.waveform_id.channel_code
    starttime = pick.time-20
    endtime = pick.time+20
    try:
        stream_contribution = client.get_waveforms(net, sta, "*", chan, starttime, endtime)
        stream = stream + stream_contribution
    except:
        pass


from toolbox import prepare_catalog_stream
st_main = prepare_catalog_stream(data_dir='/home/data/redoubt/',catalog=Catalog()+pec_event,resampling_frequency=50,tolerance=4e4)
st_main.traces.remove(st_main[11])
st_main.traces.remove(st_main[3])
st_main.traces.remove(st_main[2])
st_main.traces.remove(st_main[0])
st_copy = st_main.copy()
st = st_main.copy()
st.trim(starttime=core_time-5,endtime=core_time+13)
st = st.split()
st.detrend()
st.taper(0.05,type='cosine')
st.filter(type='bandpass',freqmin=1,freqmax=10)

f = plt.figure(figsize=(10,5))
st.plot(fig=f,equal_scale=False,show=False)
for frame in f.axes:
    frame.axes.get_yaxis().set_visible(False)
    text = frame.get_children()[1].get_text()
    text = text[0:7] + text[8:]
    frame.get_children()[1].set_text(text)
    frame.get_children()[1].set_x(0.02)
    frame.get_children()[1].set_y(0.7)
    frame.get_children()[1].set_fontsize(20)
frame = f.axes[-1]
frame.tick_params(axis='x',which='major',labelsize=18)
plt.savefig('waveforms4.png',bbox_inches = 'tight')
plt.show(bbox_inches='tight')


from eqcorrscan.utils.plotting import detection_multiplot

st = st_main.copy()
st.trim(starttime=core_time-10,endtime=core_time+10)
st = st.split()
st.detrend()
st.taper(0.05,type='cosine')
st.filter(type='bandpass',freqmin=1,freqmax=10)
st.trim(starttime=core_time-1,endtime=core_time+7)
# f = plt.figure(figsize=(4,6))
# st.plot(fig=f,equal_scale=False)
# plt.show()

chosen_detection_time = select_detection_times[5]
detection = sub_catalog[5]
st_plot = st_main.copy()
st_plot.trim(starttime=chosen_detection_time-5,endtime=chosen_detection_time+15)
st_plot = st_plot.split()
st_plot.detrend()
st_plot.taper(0.05,type='cosine')
st_plot.filter(type='bandpass',freqmin=1,freqmax=10)
f = detection_multiplot(st_plot, st, [chosen_detection_time+2.1], streamcolour='gray', templatecolour='k',size=(11,7),return_figure=True)
for frame in f.axes:
    ylabel = frame.get_ylabel()
    ylabel = 'AV.' + ylabel
    frame.text(0.02,0.7,ylabel, fontsize=20, verticalalignment='top', transform=frame.transAxes, bbox=dict(facecolor='white',alpha=0.5, edgecolor='black', boxstyle='round'))
    frame.axes.get_yaxis().set_visible(False)
frame = f.axes[-1]
frame.tick_params(axis='x',which='major',labelsize=18,labelrotation=0)
frame.set_xlabel(None)
# f.show()
f.savefig('waveforms2.png',bbox_inches = 'tight')

chosen_detection_time = select_detection_times[6]
detection = sub_catalog[6]
st_plot = st_main.copy()
st_plot.trim(starttime=chosen_detection_time-5,endtime=chosen_detection_time+15)
st_plot = st_plot.split()
st_plot.detrend()
st_plot.taper(0.05,type='cosine')
st_plot.filter(type='bandpass',freqmin=1,freqmax=10)
f = detection_multiplot(st_plot, st, [chosen_detection_time+2.1], streamcolour='gray', templatecolour='k',size=(11,7),return_figure=True)
for frame in f.axes:
    ylabel = frame.get_ylabel()
    ylabel = 'AV.' + ylabel
    frame.text(0.02,0.7,ylabel, fontsize=20, verticalalignment='top', transform=frame.transAxes, bbox=dict(facecolor='white',alpha=0.5, edgecolor='black', boxstyle='round'))
    frame.axes.get_yaxis().set_visible(False)
frame = f.axes[-1]
frame.tick_params(axis='x',which='major',labelsize=18,labelrotation=0)
frame.set_xlabel(None)
# f.show()
f.savefig('waveforms3.png',bbox_inches = 'tight')
