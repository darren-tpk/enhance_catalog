#%% SCAN DATA

# This script loads in a tribe object and scans a user-specified duration of seismic data.
# If desired, a repicked catalog can be produced during the scanning process.
# (Picks are made at times of peak channel correlation between template and detection)
# The output party and catalogs are saved in the user-defined output_dir.

# Import all dependencies
import pickle
import time
import glob
import numpy as np
from eqcorrscan import Party, Tribe, Family
from obspy import UTCDateTime, Stream, Catalog, read
from obspy.clients.fdsn import Client
from obspy.core.event import Origin
from toolbox import remove_boxcars, remove_bad_traces, calculate_catalog_FI, calculate_relative_magnitudes, reader, writer

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/mammoth/'
output_dir = main_dir + 'output/mammoth3/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
party_filename = 'party.tgz'
detected_catalog_filename = 'party.xml'
relocatable_catalog_filename = 'relocatable_catalog.xml'
repicked_catalog_filename = 'repicked_catalog.xml'
min_stations = 3                               # to remove templates that are anchored by too little stations
min_picks = 0                                  # to remove templates that are anchored by too little picks
start_time = UTCDateTime(2012, 10, 1, 0, 0, 0)  # start: UTCDateTime(2008, 5, 1, 0, 0, 0)
end_time = UTCDateTime(2013, 2, 1, 0, 0, 0)    # goal: UTCDateTime(2009, 9, 1, 0, 0, 0)
max_zeros = 100                                # maximum number of zeros allowed in daylong merged trace
npts_threshold = 100                           # drop data chunks that are overly short / have too little samples
samp_rate = 50                                 # to resample streams to match tribes
tolerance = 4e4                                # for boxcar removal
threshold_type = 'av_chan_corr'                # also consider 'MAD', Jeremy uses 12
threshold = 0.7                                # used to be 0.74 from sensitivity test. Was too high
# rethreshold = 0.70                           # filter party using revised threshold from sensitivity test
trig_int = 8                                   # trigger interval for each template. Also used to remove repeats (decluster)
parallel_process = True                      # parallel process for speed
local = True                                   # if set to True, use data from local machine
client_name = 'IRIS'                           # client name for non-local data query

# Settings for FI calculation
reference_stations = ['MINS','MDPB','OMMB','MRD','MDC','MCM','MMP','MLC','MCV','MB01','MB02','MB03','MB05','MB06','MB07','MB08','MB09','MB10','MB11']
reference_channels = ['HHZ','HHZ','HHZ','EHZ','EHZ','EHZ','EHZ','EHZ','EHZ','HHZ','HHZ','HHZ','HHZ','HHZ','HHZ','HHZ','HHZ','HHZ','HHZ']
min_match = 3  # minimum number of matches with reference channels before FI calculation
resampling_frequency = 50  # resampling frequency before fft
tolerance = 4e4  # factor to median for boxcar removal
prepick = 1.0  # s
length = 8.0  # s
lowcut = 1  # Hz
highcut = 10  # Hz
filt_order = 4  # corners
filomin = 1  # Hz
filomax = 2.5  # Hz
fiupmin = 5   # Hz
fiupmax = 10   # Hz

# Settings for magnitude calculation
noise_window = (-20.0, -prepick)  # s
signal_window = (-prepick, -prepick+length)  # s
shift_len = 1.5  # s
min_cc = 0.6
min_snr = 1

#%% Define functions

# Nil

#%% Read in tribe and clean using min_stations and min_picks

# Read in tribe
tribe = reader(create_tribe_output_dir + tribe_filename)

# Clean tribe using min_stations and min_picks
print('\nBefore cleaning:')
print(tribe)

# First remove template traces that are all nans
for t in tribe:
    for tr in t.st:
        if np.isnan(tr.data).all():
            t.st.traces.remove(tr)

# Remove based on min_stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= min_stations]
print('\nAfter removing templates with < %d stations...' % min_stations)
print(tribe)

# Remove based on min_picks
tribe.templates = [t for t in tribe if len(t.event.picks) >= min_picks]
print('\nAfter removing templates with < %d valid picks...' % min_picks)
print(tribe)

#%% If using local data, get a unique list of all template station-channel combinations for data fetching later
if local:

    # Initialize list for data fileheads
    data_fileheads = []

    # Loop through trace information in every template's stream
    for template in tribe:
        for trace in template.st:

            # Extract station and channel names from trace
            sta = trace.stats.station
            chan = trace.stats.channel

            # Craft filehead and append
            data_filehead = data_dir + sta + '.' + chan + '.'
            data_fileheads.append(data_filehead)

    # Now compile unique and sort
    data_fileheads = list(set(data_fileheads))
    data_fileheads.sort()

#%% Brute force scan: load each day's data and scan day-by-day

# Loop over month long chunks
time_list = [UTCDateTime(2012,10,1),
             UTCDateTime(2012,11,1),
             UTCDateTime(2012,12,1),
             UTCDateTime(2013,1,1),
             UTCDateTime(2013,2,1)]

# Loop over time list and save in monthly chunks
for k in range(len(time_list)-1):

    # Define start and end time rather than at start of script
    start_time = time_list[k]
    end_time = time_list[k+1]
    time_tag = '_%4d%02d%02d_%4d%02d%02d' % (start_time.year,start_time.month,start_time.day,end_time.year,end_time.month,end_time.day)
    print('NOW PROCESSING TIME TAG: ' + time_tag[1:])

    # Initialize master party and calculate number of days to scan
    party_all = Party()
    num_days = int(np.floor((end_time - start_time) / 86400))
    time_start = time.time()

    # Commence loop over days
    for i in range(num_days):

        # Define temporal boundaries of our day's scan
        t1 = start_time + (i * 86400)
        t2 = start_time + ((i + 1) * 86400)
        print('\nNow at %s...' % str(t1))

        # If using local data, read in data and process stream before creating party
        if local:

            # Initialize stream object and fetch data from data_dir
            stream = Stream()

            # Loop over data fileheads, using glob.glob to match data files we want in data_dir
            for data_filehead in data_fileheads:
                data_filename = data_filehead + str(t1.year) + ':' + f'{t1.julday :03}' + ':*'
                matching_filenames = (glob.glob(data_filename))

                # Try to read the miniseed data file and add it to the stream (occasionally fails)
                for matching_filename in matching_filenames:
                    try:
                        stream_contribution = read(matching_filename)
                        stream = stream + stream_contribution
                    except:
                        continue

            # Remove bad traces (too many zeros or too little samples)
            stream = remove_bad_traces(stream, max_zeros=max_zeros, npts_threshold=npts_threshold)

            # Process stream (remove spikes, downsample to match tribe, detrend, merge)
            stream = remove_boxcars(stream,tolerance)
            stream = stream.resample(sampling_rate=samp_rate)
            stream = stream.detrend("simple")
            stream = stream.merge()
            stream = stream.trim(starttime=t1, endtime=t2, pad=True)

            # Check if stream is already empty. If yes, declare failure and skip
            if stream is None:
                print('The stream is already empty.')
                time_current = time.time()
                print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))
                continue

            # Remove traces that are overly masked (i.e. more zeros than actual data points)
            for tr in stream:
                if len(np.nonzero(tr.data)[0]) < 0.5 * len(tr.data):
                    stream.remove(tr)

            print('Stream despiked, resampled, merged and trimmed. Getting party of detections...')

            # Attempt to scan the current day
            try:
                party = tribe.detect(
                    stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, daylong=True, overlap=None, parallel_process=True)

                # Append party to party_all
                party_all = party_all + party
                time_current = time.time()
                print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))

            # If the scan for the day fails, print a notification
            except:
                time_current = time.time()
                print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
                continue

        # If not using local data, directly query client
        else:

            # Define client
            client = Client(client_name)

            # Attempt scan via client downloads
            try:
                party = tribe.client_detect(client=client, starttime=t1, endtime=t2, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int)

                # Append party to party_all
                party_all = party_all + party
                time_current = time.time()
                print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))

            # If the scan for the day fails, print a notification
            except:
                time_current = time.time()
                print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
                continue

    # Conclude process
    time_end = time.time()
    print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))

    #%% Process catalog derived from party, then write out both the party and the catalog

    # Remove detections that are anchored by too little stations
    print('Filtering away detections anchored by too little stations...')
    for family in party_all:
        for detection_index in reversed(range(len(family.detections))):
            detection = family[detection_index]
            if detection.no_chans < min_stations:
                family.detections.pop(detection_index)

    # Clean the party off of repeats (different templates that detect the "same" event)
    print('Declustering party with leftover detections...')
    party_declustered = party_all.decluster(trig_int=trig_int)

    # Extract catalog from declustered party
    detected_catalog = party_declustered.get_catalog()

    # Write out party and catalog with all detections
    writer(scan_data_output_dir + party_filename[0:-4] + time_tag + party_filename[-4:], party_declustered)
    writer(scan_data_output_dir + detected_catalog_filename[0:-4] + time_tag + detected_catalog_filename[-4:], detected_catalog)

    # Commence FI and magnitude calculation for relocatable catalog
    print('\nNow processing FI and magnitudes of relocatable catalog...')
    time_start = time.time()

    # Extract all template names from tribe
    template_names = [template.name for template in tribe]

    # Initialize new catalog that stores detections that can be relocated
    relocatable_catalog = Catalog()

    # Loop over the detected catalog, copying over events if their templates have locations
    # Also give each event pick phase information based on channel convention
    for event in detected_catalog:

        # Get source event
        source_template_name = event.comments[0].text.split(' ')[1]
        source_template_index = template_names.index(source_template_name)
        source_template = tribe[source_template_index]
        source_event = source_template.event

        # Check if the source event has a location. Only continue if it does
        if source_event.origins[0].latitude is not None:

            # We extract the source event's location (lat, lon, dep)
            lat = source_event.origins[0].latitude
            lon = source_event.origins[0].longitude
            dep = source_event.origins[0].depth

            # Reformat the event's detection time to UTC
            time_text = event.resource_id.id.split('_')[-1]
            UTCdt = UTCDateTime(time_text)

            # Reweight arrivals and add phase weight information [TO BE DONE]
            pass

            # Create an ObsPy Origin object and store it in the event
            event.origins = [Origin(time= UTCdt, latitude = lat, longitude = lon, depth = dep)]

            # Loop over picks in event
            for pick in event.picks:

                # Check for non-vertical channel, and give an 'S' hint
                if pick.waveform_id.channel_code[-1] == 'N' or pick.waveform_id.channel_code[-1] == 'E':
                    pick.phase_hint = 'S'

                # Also check for a vertical channel, and give an 'P' hint
                elif pick.waveform_id.channel_code[-1] == 'Z':
                    pick.phase_hint = 'P'

                # Throw out a ValueError if the channel code is not recognized
                else:
                    raise ValueError

            # Append to new catalog
            relocatable_catalog += event

    # Write out catalog with relocatable detections
    writer(scan_data_output_dir + relocatable_catalog_filename[0:-4] + time_tag + relocatable_catalog_filename[-4:], relocatable_catalog)

    # Calculate FI for relocatable catalog
    print('\nCalculating FIs...')
    relocatable_catalog = calculate_catalog_FI(relocatable_catalog, data_dir, reference_stations, reference_channels,
                                               min_match, resampling_frequency, tolerance, prepick, length, lowcut,
                                               highcut, filt_order, filomin, filomax, fiupmin, fiupmax, verbose=True,
                                               histogram=False)
    print('\nFIs calculated. Now calculating magnitudes...')

    # Calculate magnitude for relocatable catalog
    print('\nCalculating relative magnitudes...')
    relocatable_catalog = calculate_relative_magnitudes(relocatable_catalog, tribe, data_dir, noise_window, signal_window,
                                                        min_cc, min_snr, shift_len, tolerance, samp_rate, lowcut,
                                                        highcut, filt_order, verbose=True)
    print('\nMagnitudes calculated. Saving catalog...')

    # Write out catalog with relocatable detections
    writer(scan_data_output_dir + relocatable_catalog_filename[0:-4] + time_tag + relocatable_catalog_filename[-4:], relocatable_catalog)
    time_end = time.time()
    print('\nRelocatable catalog saved. Time taken: %.2f hours' % ((time_end - time_start)/3600))

# Now combine all time tag catalogs to get the full party catalog and relocatable catalog

print('\nCombining all sub-catalogs...')
time_start = time.time()

full_party_catalog = Catalog()
full_relocatable_catalog = Catalog()

# Loop over time list and save in monthly chunks
for k in range(len(time_list)-1):
    start_time = time_list[k]
    end_time = time_list[k+1]
    time_tag = '_%4d%02d%02d_%4d%02d%02d' % (start_time.year,start_time.month,start_time.day,end_time.year,end_time.month,end_time.day)
    sub_party_catalog = reader(scan_data_output_dir + detected_catalog_filename[0:-4] + time_tag + detected_catalog_filename[-4:])
    sub_relocatable_catalog = reader(scan_data_output_dir + relocatable_catalog_filename[0:-4] + time_tag + relocatable_catalog_filename[-4:])
    full_party_catalog += sub_party_catalog
    full_relocatable_catalog += sub_relocatable_catalog

writer(scan_data_output_dir + detected_catalog_filename, full_party_catalog)
writer(scan_data_output_dir + relocatable_catalog_filename, full_relocatable_catalog)
print('\nCombination complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))
