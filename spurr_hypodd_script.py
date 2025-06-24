### HypoDD Tutorial for Spurr Seismicity (December 2024 to January 2025)

# Import dependencies
import requests
from obspy import UTCDateTime
from obspy.geodetics import kilometers2degrees
from obspy.core.event.catalog import read_events, Catalog
from obspy.core.event.base import Comment
from io import BytesIO
import numpy as np
import csv
from toolbox import reader, writer
from functions import initialize_run, download_data

### Initialize run (creates sub-directories)
subdir_name = 'spurr'
initialize_run(subdir_name)

### Query catalog from FDSN

# Query parameters
starttime = UTCDateTime(2024, 12, 1)
endtime = UTCDateTime(2025, 2, 1)
volcano_latitude = 61.2989  # spurr
volcano_longitude = -152.2539  # spurr
radius_km = 15
only_finalized = True

# Set up session
session = requests.Session()
session.auth = ('username', 'password')
catalog_params = {'format': 'csv',
                 'starttime': starttime.strftime('%Y-%m-%dT%H:%M:%S.000'),
                 'endtime': endtime.strftime('%Y-%m-%dT%H:%M:%S.000'),
                 'latitude': volcano_latitude,
                 'longitude': volcano_longitude,
                 'maxradius': kilometers2degrees(radius_km),
                 'eventtype': 'eq,lp'}
catalog_request = session.get('https://avoaqmspp2.avo.alaska.edu/fdsnws/event/1/query', params=catalog_params)
catalog_request.raise_for_status()
csv_dict = csv.DictReader(catalog_request.iter_lines(decode_unicode=True))

# Store event ids, types, and statuses
evids = []
types = []
status = []
count = 0

for event_dict in csv_dict:
    count += 1
    evids.append(event_dict['id'])
    types.append(event_dict['type'])
    status.append(event_dict['status'])
    print('Event ID: %s, Time: %s, Magnitude: %s, Type: %s, Status: %s, RMS: %s' %
          (event_dict['id'], event_dict['time'], event_dict['mag'], event_dict['type'], event_dict['status'], event_dict['rms']))
print('There are %d events in the requested catalog.\n' % count)

# Initialize catalog and populate with events
catalog = Catalog()
for i, evid in enumerate(evids):
    if only_finalized and status[i] != 'F':
        print('Event ID %s skipped, event not finalized.' % (evid))
        continue
    event_params = {'eventid': evid}
    event_request = session.get('https://avoaqmspp2.avo.alaska.edu/fdsnws/event/1/query', params=event_params)
    event_request.raise_for_status()
    if len(read_events(BytesIO(event_request.content), 'QUAKEML').events) > 0:
        event = read_events(BytesIO(event_request.content), 'QUAKEML').events[0]
        if event.comments:
            print('Event ID %s pull sucessful: %s from %s' % (evid, event.event_type, event.comments[0].text))
            event.comments[-1] = Comment(text='type:%s' % types[i])
            catalog.append(event)
    else:
        print('Event ID %s pull failed' % evid)
print('Catalog populated with %d events.\n' % len(catalog))

# Remove events with no magnitude
catalog_withmag = Catalog([e for e in catalog if e.magnitudes])
print('Removed %d events with no magnitude.' % (len(catalog) - len(catalog_withmag)))
catalog = catalog_withmag
print('Catalog now contains %d events.\n' % len(catalog))

# Write catalog
writer('./spurr_20241201_20250201.xml', catalog)

### Filter catalog by phase picks, RMSE, and gap
# catalog = reader('./spurr_20241201_20250201.xml')

# Filter parameters
min_P = 4
min_S = 2
max_rmse = 0.30
max_gap = 150

# Filter catalog
print('Initial catalog size: %d' % len(catalog))
filtered_catalog = Catalog()
for e in catalog:
    try:
        phases = [a.phase for a in e.preferred_origin().arrivals]
        num_P = phases.count('P')
        num_S = phases.count('S')
        rmse = e.preferred_origin().quality.standard_error
        gap = e.preferred_origin().quality.azimuthal_gap
        if num_P >= min_P and num_S >= min_S and rmse <= max_rmse and gap <= max_gap:
            filtered_catalog.events.append(e)
    except:
        continue
print('Filtered catalog size: %d' % len(filtered_catalog))

# Transfer phase hint to filtered catalog
for e in filtered_catalog:
    arrivals = e.preferred_origin().arrivals
    pick_ids = [p.resource_id for p in e.picks]
    for arrival in arrivals:
        try:
            pick_index = pick_ids.index(arrival.pick_id)
            e.picks[pick_index].phase_hint = arrival.phase
        except:
            continue
    P_to_keep = [p for p in e.picks if p.phase_hint == 'P' and p.waveform_id.channel_code[-1] == 'Z']
    S_to_keep = [p for p in e.picks if p.phase_hint == 'S' and p.waveform_id.channel_code[-1] != 'Z']
    e.picks = P_to_keep + S_to_keep

# Additionally filter picks by retaining only stations with > 50 event observations
all_stations = []
all_stachans = []
for e in filtered_catalog:
    stachans = [p.waveform_id.station_code + '.' + p.waveform_id.channel_code for p in e.picks]
    all_stachans = all_stachans + stachans
unique_stachans = list(np.sort(np.unique(all_stachans)))
unique_stachan_count = [all_stachans.count(stachan) for stachan in unique_stachans]
stachans_keep = [unique_stachans[i] for i in range(len(unique_stachans)) if unique_stachan_count[i] >= 50]

relocatable_catalog = Catalog()
for e in filtered_catalog:
    picks_to_keep = [p for p in e.picks if (p.waveform_id.station_code + '.' + p.waveform_id.channel_code) in stachans_keep]
    e.picks = picks_to_keep
    relocatable_catalog.events.append(e)

### Download data to local directory

# Parse station and channel list
networks = ['AV'] * len(stachans_keep)
stations = [s.split('.')[0] for s in stachans_keep]
channels = [s.split('.')[1] for s in stachans_keep]
locations = ['--'] * len(stations)
stations_str = ','.join(stations)
channels_str = ','.join(channels)

# Define data download parameters
data_destination = './data/spurr/'
starttime = UTCDateTime(2025,1,1)
endtime = UTCDateTime(2025,2,1)
client = 'IRIS'
network = ','.join(networks)
station = ','.join(stations)
channel = ','.join(channels)
location = ','.join(locations)

download_data(data_destination=data_destination,
              starttime=starttime,
              endtime=endtime,
              client=client,
              network=network,
              station=station,
              channel=channel,
              location=location)

#############################################################################
import pandas as pd
from toolbox import reader
from functions import generate_dtcc, run_hypoDD, plot_hypoDD_results

# relocatable_catalog = reader('./spurr_20241201_20250201_relocatable_catalog.xml')
relocate_catalog_output_dir = './output/spurr/relocate_catalog/'
data_path = './data/spurr/'
pre_pick_actual = 0.15
pre_pick_excess = 1
length_actual = 2.66
length_excess = 4
shift_len = 0.4
resampling_frequency = 100
min_link = 3
min_cc = 0.4
fmin = 1
fmax = 10

generate_dtcc(catalog=relocatable_catalog,
              relocate_catalog_output_dir=relocate_catalog_output_dir,
              data_path=data_path,
              pre_pick_actual=pre_pick_actual,
              pre_pick_excess=pre_pick_excess,
              length_actual=length_actual,
              length_excess=length_excess,
              shift_len=shift_len,
              lowcut=fmin,
              highcut=fmax,
              resampling_frequency=resampling_frequency,
              min_cc=min_cc,
              min_link=min_link)

## (8) Relocate catalog candidates
D = 160
ph2dt_inc_dict = {'MEV': 12000,
                  'MSTA': 2600,
                  'MOBS': 12000}
ph2dt_inp_dict = {'MINWGHT': 0,
                  'MAXDIST': 200,
                  'MAXSEP': 10,
                  'MAXNGH': 10,
                  'MINLNK': 3,
                  'MINOBS': 3,
                  'MAXOBS': 50}
hypoDD_inc_dict = {'MAXEVE': 12000,
                   'MAXDATA': 5000000,
                   'MAXEVE0': 2,
                   'MAXDATA0': 1,
                   'MAXLAY': 30,
                   'MAXSTA': 100,
                   'MAXCL': 200}
hypoDD_inp_dict = {'IDAT': 3,
                   'IPHA': 3,
                   'DIST': 3000,
                   'OBSCC': 0,
                   'OBSCT': 0,
                   'MINDS': -999,
                   'MAXDS': -999,
                   'MAXGAP': -999,
                   'ISTART': 2,
                   'ISOLVE': 2,
                   'IAQ': 0,
                   'NSET': 5,
                   'NITER': [10, 10, 10, 10, 10],
                   'WTCCP': [0.01, 0.01, 0.10, 0.50, 1.00],
                   'WTCCS': [0.01, 0.01, 0.10, 0.50, 1.00],
                   'WRCC': [10, 10, 10, 6, 6],
                   'WDCC': [4, 4, 4, 2, 2],
                   'WTCTP': [1.00, 1.00, 1.00, 0.10, 0.01],
                   'WTCTS': [1.00, 1.00, 1.00, 0.10, 0.01],
                   'WRCT': [12, 6, 6, 6, 6],
                   'WDCT': [10, 5, 2.5, 2.5, 2.5],
                   'DAMP': [D, D, D, D, D]}

# relocatable_catalog = reader('./spurr_20241201_20250201_relocatable_catalog.xml')
stations = 'RDT,SPBG,SPBL,SPCG,SPCL,SPCN,SPCP,SPNN,SPU,SPWE,STLK'
# stations = ','.join(np.unique(stations))
avosta_df = pd.read_csv('./avo_stations.csv')
stalats = ','.join([str(avosta_df.Latitude[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
stalons = ','.join([str(avosta_df.Longitude[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
staelevs = ','.join([str(avosta_df.Elevation[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
vzmodel_path = './vz/spurr_vzmodel.txt' # './vz/spurr_vzmodel2.txt'
has_ps_ratio = True  # If the input vz model has p/s ratios in the 3rd column
correct_depths = True

# ph2dt is bugged
hypoDD_loc, hypoDD_reloc = run_hypoDD(catalog=relocatable_catalog,
                                      relocate_catalog_output_dir=relocate_catalog_output_dir,
                                      stations=stations,
                                      stalats=stalats,
                                      stalons=stalons,
                                      staelevs=staelevs,
                                      vzmodel_path=vzmodel_path,
                                      has_ps_ratio=has_ps_ratio,
                                      correct_depths=correct_depths,
                                      ph2dt_inc_dict=ph2dt_inc_dict,
                                      ph2dt_inp_dict=ph2dt_inp_dict,
                                      hypoDD_inc_dict=hypoDD_inc_dict,
                                      hypoDD_inp_dict=hypoDD_inp_dict)

## (9) Plot relocated earthquakes
# relocatable_catalog = reader('./spurr_20241201_20250201_relocatable_catalog.xml')
hypoDD_loc = reader(relocate_catalog_output_dir + 'hypoDD_loc.xml')
hypoDD_reloc = reader(relocate_catalog_output_dir + 'hypoDD_reloc.xml')
lat_lims = [61.25, 61.375]
lon_lims = [-152.5, -152.15]
dep_lims = [-5, 10]
markersize = 2
figsize = (10,7)

plot_hypoDD_results(hypoDD_in=relocatable_catalog,
                    hypoDD_loc=hypoDD_loc,
                    hypoDD_reloc=hypoDD_reloc,
                    lat_lims=lat_lims,
                    lon_lims=lon_lims,
                    dep_lims=dep_lims,
                    markersize=markersize,
                    figsize=figsize,
                    legend_loc='lower left')
