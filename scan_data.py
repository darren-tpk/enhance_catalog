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
from eqcorrscan import Party, Tribe
from obspy import UTCDateTime, Stream, Catalog, read
from toolbox import remove_boxcars, reader, writer

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
data_dir = '/home/data/redoubt/'
output_dir ='/home/ptan/enhance_catalog/output/'
convert_redpy_output_dir = output_dir + 'convert_redpy/'
create_tribe_output_dir = output_dir + 'create_tribe/'
tribe_filename = 'tribe_test.tgz'
scan_data_output_dir = output_dir + 'scan_data/'
party_filename = 'party_test.tgz'
catalog_filename = 'party_catalog_test.xml'
repicked_catalog_filename = None
min_stations = 3                               # to remove templates that are anchored by too little stations
min_picks = 0                                  # to remove templates that are anchored by too little picks
start_time = UTCDateTime(2009, 3, 20, 0, 0, 0) # start: UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 3, 21, 0, 0, 0)    # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50                                 # to resample streams to match tribes
tolerance = 4e4                                # for boxcar removal
threshold_type = 'av_chan_corr'                # also consider 'MAD', Jeremy uses 12
threshold = 0.60                               # used to be 0.74 from sensitivity test. Was too high
trig_int = 10                                  # trigger interval for each template. Also used to remove repeats (decluster)
parallel_process = 'True'                      # parallel process for speed
generate_repicked_catalog = False              # option to use lag_calc to do catalog repicking

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

#%% If we are generating a repicked catalog, we need to exclude templates without a location
if generate_repicked_catalog:

    # Unpickle the list of unassociated template numbers for comparison
    with open(convert_redpy_output_dir + 'unassociated_clusters.txt', 'rb') as cluster_pickle:
        unassociated_clusters = pickle.load(cluster_pickle)
        #associated_clusters = pickle.load(cluster_pickle)

    # Initialize a master repicked catalog for appending in main loop later
    master_repicked_catalog = Catalog()

    # Initialize placeholder empty tribe
    new_tribe = Tribe()

    # Loop through templates
    for template in tribe:

        # Extract cluster number to compare with unassociated cluster list
        cluster_number = int(template.event.origins[0].comments[0].text.split(' ')[1])

        # Only carry over templates that are associated to the pre-existing catalog
        if cluster_number not in unassociated_clusters:
            new_tribe += template

    # Replace tribe in main workflow for next steps (the templates now all point to events with location)
    tribe = new_tribe

#%% Get a unique list of all template station-channel combinations for data fetching later

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

    # Process stream (remove spikes, downsample to match tribe, detrend)
    stream = remove_boxcars(stream,tolerance)
    stream = stream.resample(sampling_rate=samp_rate)
    stream = stream.detrend("simple")
    stream = stream.merge()
    print('Stream despiked, resampled and merged. Getting party of detections...')

    # Attempt to scan the current day
    try:
        party = tribe.detect(
            stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, daylong=True, overlap=None, parallel_process=True)

        # Append party to party_all
        party_all = party_all + party
        time_current = time.time()
        print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))

        # If generating a repicked catalog, we prep the stream differently for lag_calc
        if generate_repicked_catalog:

            # Copy the party and stream, then re-merge the stream using method 1
            party_calc = party.copy()
            stream_gappy = stream.copy()
            stream_gappy = stream_gappy.split()
            stream_gappy = stream_gappy.merge(method=1)

            # Use lag_calc to produce the repicked catalog section
            repicked_catalog = party_calc.lag_calc(stream_gappy, pre_processed=False, shift_len=3, min_cc=0.7,
                                                   export_cc=True, cc_dir='/home/ptan/attempt_eqcorrscan/output/cc_data')

            # Add the repicked catalog section to the master repicked catalog
            master_repicked_catalog += repicked_catalog
            time_current = time.time()
            print('Repicked catalog generated. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))

    # If the scan for the day fails, print a notification
    except:
        time_current = time.time()
        print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
        continue

# Conclude process
time_end = time.time()
print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))
if generate_repicked_catalog:
    writer(scan_data_output_dir + repicked_catalog_filename, master_repicked_catalog)
    print('Number of relocated events with picks: %d out of %d total' % (len([1 for event in master_repicked_catalog if (event.picks != [])]), len(master_repicked_catalog)))

#%% Clean the party off of repeats (different templates that detect the "same" event)
party_all = party_all.decluster(trig_int=trig_int)

#%% Write our results in party form and in catalog form

# Write out party
writer(scan_data_output_dir + party_filename, party_all)

# Write out catalog
detected_catalog = party_all.get_catalog()
writer(scan_data_output_dir + catalog_filename, detected_catalog)
