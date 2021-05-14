# generate data products
# code still needs to be updated

# import packages and functions that we need
import glob

import numpy as np
from eqcorrscan import Party
from obspy import UTCDateTime, Stream, read

from toolbox import remove_repeats, reader, writer, remove_boxcars

# define all variables here
time_interval = 5  # number of seconds before and after each detection to check for repeats

# grab our party
party = party_all
#party = reader('/home/ptan/attempt_eqcorrscan/output/party_marchswarm.tgz')
# unpickle our cluster lists (NA: not associated to AVO event, AVO: associated to AVO event)
with open('/home/ptan/attempt_eqcorrscan/avo_data/redpy/clusters_NA_AVO.txt', 'rb') as cluster_pickle:   # Unpickling
    clusters_NA = pickle.load(cluster_pickle)
    clusters_AVO = pickle.load(cluster_pickle)
# now initialize an empty party, and fill it with families that come from a template assocaited with an AVO event
party_associated = Party()
for family in party:
    cluster_number = int(family.template.event.origins[0].comments[0].text.split(' ')[1])
    if cluster_number not in clusters_NA:
        party_associated = party_associated + family
total_detect_num = sum([len(family) for family in party_associated])
print('Party of associated detections created. Number of detections =',int(total_detect_num))

# now attempt lag_calc
# define some params to download appropriate data
start_time = UTCDateTime(2009, 3, 18, 0, 0, 0)  # start: UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 3, 27, 0, 0, 0)  # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50             # to resample streams to match tribes
tolerance = 4e4             # for boxcar removal
# get unique list of all template station & channel combinations
data_dir = '/home/data/redoubt/'
data_fileheads = []
for family in party_associated:
    for trace in family.template.st:
        # extract key information from trace
        sta = trace.stats.station
        chan = trace.stats.channel
        # craft file string and append
        data_filehead = data_dir + sta + '.' + chan + '.'
        data_fileheads.append(data_filehead)
# now compile unique and sort
data_fileheads = list(set(data_fileheads))
data_fileheads.sort()
# no loop through days and download data
num_days = int(np.floor((end_time - start_time) / 86400))
for i in range(num_days):
    # define boundaries of our day's scan
    t1 = start_time + (i * 86400)
    t2 = start_time + ((i + 1) * 86400)
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
stream = stream.merge(method=1)
print('Stream despiked, resampled and merged. Getting party of detections...')
# now run lag_calc
repicked_catalog = party_associated.lag_calc(stream, pre_processed=False, shift_len=0.5, min_cc=0.4)
# party = reader('/home/ptan/attempt_eqcorrscan/output/party_marchswarm.tgz')

# derive a party that only has associated cores


# clean party off of repeated detections
party_clean = remove_repeats(party,time_interval)



# save cleaned party
party_clean_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_clean_outpath+'party_clean.tgz',party_clean)