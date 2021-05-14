# scan data to make detections

# import packages that we need
import pickle
import time
import glob
import numpy as np
from eqcorrscan import Party, Tribe
from obspy import UTCDateTime, Stream, Catalog, read
# import toolbox functions
from toolbox import remove_boxcars, reader, writer

# define all variables here
start_time = UTCDateTime(2009, 3, 18, 0, 0, 0)  # start: UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 3, 27, 0, 0, 0)  # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50             # to resample streams to match tribes
tolerance = 4e4             # for boxcar removal
threshold = 0.65             # Used to be 20 without RDWB and RDJH
threshold_type = 'av_chan_corr'
trig_int = 10
parallel_process = 'True'
generate_repicked_catalog = 'True'

# read tribe data
tribe = reader('/home/ptan/attempt_eqcorrscan/output/tribe_10.tgz')

# filter away templates with < 3 stations
print('\nBefore filter:')
print(tribe)
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 3]
print('\nAfter removing templates with < 3 stations...')
print(tribe)

# # filter away templates with < 3 valid picks (with data stream attached)
# print('\nBefore filter:')
# print(tribe)
# new_tribe = Tribe()
# for template in tribe:
#     if len(template.st) >= 3:
#         new_tribe += template
# tribe = new_tribe
# print('\nAfter removing templates with < 3 valid picks...')
# print(tribe)

# if we are generating the repicked catalog for the relocation track,
# we want to exclude templates that are not associated to AVO
if generate_repicked_catalog:
    # we have to unpickle the list of unassociated template numbers for comparison later
    with open('/home/ptan/attempt_eqcorrscan/avo_data/redpy/clusters_NA_AVO.txt', 'rb') as cluster_pickle:  # Unpickling
        clusters_NA = pickle.load(cluster_pickle)
        #clusters_AVO = pickle.load(cluster_pickle)
    # we also initialize a master repicked catalog for appending in main loop later
    master_repicked_catalog = Catalog()
    # initialize placeholder empty tribe
    new_tribe= Tribe()
    # loop through templates
    for template in tribe:
        # extract cluster number for debugging / reference for relocation track
        cluster_number = int(template.event.origins[0].comments[0].text.split(' ')[1])
        # for relocation track, we append templates that are associated to AVO
        if cluster_number not in clusters_NA:
            new_tribe += template
# replace tribe in main workflow for next steps
tribe = new_tribe

# get unique list of all template station & channel combinations
data_dir = '/home/data/redoubt/'
data_fileheads = []
for template in tribe:
    for trace in template.st:
        # extract key information from trace
        sta = trace.stats.station
        chan = trace.stats.channel
        # craft file string and append
        data_filehead = data_dir + sta + '.' + chan + '.'
        data_fileheads.append(data_filehead)
# now compile unique and sort
data_fileheads = list(set(data_fileheads))
data_fileheads.sort()

# brute force scan: load each day's data and scan day by day
party_all = Party()
num_days = int(np.floor((end_time - start_time) / 86400))
time_start = time.time()
for i in range(num_days):
    # define boundaries of our day's scan
    t1 = start_time + (i * 86400)
    t2 = start_time + ((i + 1) * 86400)
    print('\nNow at %s...' % str(t1))
    # initialize and read in stream
    stream = Stream()
    for data_filehead in data_fileheads:
        data_filename = data_filehead + str(t1.year) + ':' + f'{t1.julday :03}' + ':*'
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                stream_contribution = read(matching_filename)
                stream = stream + stream_contribution
            except:
                continue
    # process stream (remove spikes, downsample to match tribe, detrend)
    stream = remove_boxcars(stream,tolerance)
    stream = stream.resample(sampling_rate=samp_rate)
    stream = stream.detrend("simple")
    stream = stream.merge()
    print('Stream despiked, resampled and merged. Getting party of detections...')
    # scan the current day
    try:
        party = tribe.detect(
            stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, daylong=True, overlap=None, parallel_process='True')
        # append party to party_all
        party_all = party_all + party
        time_current = time.time()
        print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start)/3600))
        # for relocation track
        if generate_repicked_catalog:
            # we copy the party, re-merge the stream, and generate repicked catalog
            party_calc = party.copy()
            stream_gappy = stream.copy()
            stream_gappy = stream_gappy.split()
            stream_gappy = stream_gappy.merge(method=1)
            # use lag_calc to produce repicked catalog
            repicked_catalog = party_calc.lag_calc(stream_gappy, pre_processed=False, shift_len=3, min_cc=0.7,
                                                   export_cc=True, cc_dir='/home/ptan/attempt_eqcorrscan/output/cc_data')
            # add repicked catalog to master repicked catalog
            master_repicked_catalog += repicked_catalog
            print('Repicked catalog generated. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))
    except:
        print('Party failed, skipping.')
        continue
time_end = time.time()
print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))
print('Number of relocated events with picks: %d out of %d total' % (len([1 for event in master_repicked_catalog if (event.picks != [])]), len(master_repicked_catalog)))

# save master repicked catalog
if generate_repicked_catalog:
    catalog_outpath = '/home/ptan/attempt_eqcorrscan/output/'
    writer(catalog_outpath+'repicked_catalog.xml', master_repicked_catalog)

# write the combined party object as a tar file
party_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_outpath+'party_marchswarm.tgz',party_all)


