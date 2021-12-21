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
data_dir = '/home/data/redoubt/'
# main_dir = './enhance_catalog/'
# data_dir = './data/redoubt/'
output_dir = main_dir + 'output/redoubt2/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
party_filename = 'party.tgz'
detected_catalog_filename = 'party_catalog.xml'
relocatable_catalog_filename = 'relocatable_catalog.xml'
repicked_catalog_filename = 'repicked_catalog.xml'
min_stations = 3                               # to remove templates that are anchored by too little stations
min_picks = 0                                  # to remove templates that are anchored by too little picks
# start_time = UTCDateTime(2008, 10, 1, 0, 0, 0)  # start: UTCDateTime(2008, 5, 1, 0, 0, 0)
# end_time = UTCDateTime(2009, 5, 1, 0, 0, 0)    # goal: UTCDateTime(2009, 9, 1, 0, 0, 0)
max_zeros = 100                                # maximum number of zeros allowed in daylong merged trace
npts_threshold = 100                           # drop data chunks that are overly short / have too little samples
samp_rate = None                                 # to resample streams to match tribes
tolerance = 4e4                                # for boxcar removal
threshold_type = 'av_chan_corr'                # also consider 'MAD', Jeremy uses 12
threshold = 0.60                               # used to be 0.74 from sensitivity test. Was too high
rethreshold = 0.70                             # filter party using revised threshold from sensitivity test
trig_int = 8                                   # trigger interval for each template. Also used to remove repeats (decluster)
parallel_process = 'True'                      # parallel process for speed
generate_repicked_catalog = False              # option to use lag_calc to do catalog repicking
local = True                                  # if set to True, use data from local machine
client_name = 'IRIS'                           # client name for non-local data query

# Settings for FI calculation
reference_station = 'REF'
reference_channel = 'EHZ'
prepick = 1.0  # s
length = 8.0  # s
filomin = 2  # Hz
filomax = 4  # Hz
fiupmin = 8   # Hz
fiupmax = 25   # Hz

# Settings for magnitude calculation
noise_window = (-20.0, -prepick)  # s
signal_window = (-prepick, -prepick+length)  # s
shift_len = 1.5  # s
min_cc = 0.7
min_snr = 2

#%% Define functions

# Nil

#%% Read in tribe and clean using min_stations and min_picks

# Read in tribe
tribe = reader(create_tribe_output_dir + tribe_filename)

# Clean tribe using min_stations and min_picks
print('\nBefore cleaning:')
print(tribe)

# Remove based on min_stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= min_stations]
print('\nAfter removing templates with < %d stations...' % min_stations)
print(tribe)

# Remove based on min_picks
new_tribe = Tribe()
for template in tribe:
    if len(template.st) >= min_picks:
        new_tribe += template
tribe = new_tribe
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
# time_list = [UTCDateTime(2009, 1, 15, 0, 0, 0),
#              UTCDateTime(2009, 2, 1, 0, 0, 0),
#              UTCDateTime(2009, 2, 15, 0, 0, 0),
#              UTCDateTime(2009, 3, 1, 0, 0, 0),
#              UTCDateTime(2009, 3, 15, 0, 0, 0),
#              UTCDateTime(2009, 4, 1, 0, 0, 0),
#              UTCDateTime(2009, 4, 15, 0, 0, 0),
#              UTCDateTime(2009, 5, 1, 0, 0, 0)]
# time_list = [UTCDateTime(2009, 5, 1, 0, 0, 0),
# 	         UTCDateTime(2009, 5, 15, 0, 0, 0),
#              UTCDateTime(2009, 6, 1, 0, 0, 0),
#              UTCDateTime(2009, 6, 15, 0, 0, 0),
#              UTCDateTime(2009, 7, 1, 0, 0, 0),
#              UTCDateTime(2009, 7, 15, 0, 0, 0),
#              UTCDateTime(2009, 8, 1, 0, 0, 0),
#              UTCDateTime(2009, 8, 15, 0, 0, 0),
#              UTCDateTime(2009, 9, 1, 0, 0, 0)]

# time_list = [UTCDateTime(2009, 3, 15, 0, 0, 0),
#              UTCDateTime(2009, 4, 1, 0, 0, 0),
#              UTCDateTime(2009, 4, 15, 0, 0, 0),
#              UTCDateTime(2009, 5, 1, 0, 0, 0),
# 	         UTCDateTime(2009, 5, 15, 0, 0, 0),
#              UTCDateTime(2009, 6, 1, 0, 0, 0),
#              UTCDateTime(2009, 6, 15, 0, 0, 0),
#              UTCDateTime(2009, 7, 1, 0, 0, 0),
#              UTCDateTime(2009, 7, 15, 0, 0, 0),
#              UTCDateTime(2009, 8, 1, 0, 0, 0),
#              UTCDateTime(2009, 8, 15, 0, 0, 0),
#              UTCDateTime(2009, 9, 1, 0, 0, 0)]

time_list = [UTCDateTime(2009,5,1), UTCDateTime(2009,5,15)]

# Initialize master catalog
master_catalog = Catalog()

# Loop over time list and save in monthly chunks
for k in range(len(time_list)-1):

    # Define start and end time rather than at start of script
    start_time = time_list[k]
    end_time = time_list[k+1]
    time_tag = '_%4d%02d%02d_%4d%02d%02d' % (start_time.year,start_time.month,start_time.day,end_time.year,end_time.month,end_time.day)
    print('NOW PROCESSING TIME TAG: ' + time_tag[1:])

    # Load catalog with current time tag
    add_catalog = reader(scan_data_output_dir + relocatable_catalog_filename[0:-4] + time_tag + relocatable_catalog_filename[-4:])
    print('Time Tag: ' + time_tag + ', Number of events before filtering: ' + str(len(add_catalog)))

    # Implement threshold filter
    for event_index in reversed(range(len(add_catalog.events))):
        event = add_catalog[event_index]
        av_chan_corr = float(event.comments[2].text.split('=')[1]) / event.comments[3].text.count(')')
        if av_chan_corr < min_cc:
            add_catalog.events.pop(event_index)

    # Print number of events after threshold filter
    print('Time Tag: ' + time_tag + ', Number of events after filtering: ' + str(len(add_catalog)))

    # Calculate magnitude for relocatable catalog
    add_catalog = calculate_relative_magnitudes(add_catalog, tribe, data_dir, noise_window, signal_window,
                                                min_cc, min_snr, shift_len, tolerance, samp_rate)
    print('Magnitudes successfully calculated. Adding catalog.')

    # Append to master catalog
    master_catalog = master_catalog + add_catalog

    # Write out master catalog at every step to save progress
    writer(scan_data_output_dir + 'master_catalogM2.xml', master_catalog)

# all_times = [event.origins[0].time for event in master_catalog]