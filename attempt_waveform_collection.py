from obspy import UTCDateTime
from waveform_collection import gather_waveforms, gather_waveforms_bulk

# Seismic Data
STARTTIME = UTCDateTime(2005, 7, 3, 10, 18, 22, 110000)
ENDTIME   = UTCDateTime(2005, 7, 3, 16, 18, 22, 110000)

st = gather_waveforms(source='IRIS', network='AV', station='AUL',
                      location='*', channel='BHN', starttime=STARTTIME,
                      endtime=ENDTIME)

## Infrasound
#STARTTIME = UTCDateTime(2020, 10, 7, 9, 41)
#ENDTIME   = UTCDateTime(2020, 10, 7, 9, 50)
#
#st = gather_waveforms(source='IRIS', network='AV', station='KENI',
#                      location='01', channel='HDF', starttime=STARTTIME,
#                      endtime=ENDTIME)

# Plot and get stats
st.plot(color='black')
st[0].stats

# Response removal
st2 = st
st2.remove_response(output = 'ACC', plot = False)

# Filter
st3 = st2
st3.filter('bandpass',freqmin=1,freqmax=20)
st3.plot(color='black')

# Trim
st4 = st3
STARTTIME2 = STARTTIME + 11.999*60*60
ENDTIME2 = ENDTIME - 11.999*60*60
st4.trim(starttime=STARTTIME2,endtime=ENDTIME2)
st4.plot(color='black')

# Decimate
st5 = st4
st5.decimate(10)
st5.plot(color='black')

# ROSES tutorial
import obspy
from obspy.clients.fdsn import Client

# Edit client to use your data center of interest
client = Client("IRIS")

# Edit this to request metadata from your favorite station(s)
t1 = obspy.UTCDateTime("2020-07-01")
inv = client.get_stations(network="IW", station="PLID", channel="BHZ", level="response", starttime=t1)
inv += client.get_stations(network="GS", station="PR01", channel="HHZ", level="response", starttime=t1)
# may get a warning about StationXML1.1 -- OK to ignore it
