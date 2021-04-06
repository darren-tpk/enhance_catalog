# scan data to make detections

# import packages that we need
import time
import glob
import numpy as np
from eqcorrscan import Party, Tribe
from obspy import UTCDateTime, Stream, read
# import toolbox functions
from toolbox import remove_boxcars, reader, writer

# define all variables here
start_time = UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 5, 1, 0, 0, 0)  # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50             # to resample streams to match tribes
threshold = 30             # Used to be 20 without RDWB and RDJH
threshold_type = 'MAD'
trig_int = 30
parallel_process = 'True'
tolerance = 4e4 # for boxcar removal

# read tribe data
tribe = reader('/home/ptan/attempt_eqcorrscan/output/tribe.tgz')

# # filter away templates with < 3 stations
# print('\nBefore filter:')
# print(tribe)
# tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 3]
# print('\nAfter removing templates with < 3 stations...')
# print(tribe)

# filter away templates with < 3 valid picks (with data stream attached)
print('\nBefore filter:')
print(tribe)
new_tribe = Tribe()
for template in tribe:
    if len(template.st) >= 3:
        new_tribe += template
tribe = new_tribe
print('\nAfter removing templates with < 3 valid picks...')
print(tribe)

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
    stream = stream.resample(sampling_rate=50.0)
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
    except:
        print('Party failed, skipping.')
        continue
time_end = time.time()
print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start)/3600))

# write the combined party object as a tar file
party_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_outpath+'party.tgz',party_all)

# # use client detect and compare
# from obspy.clients.fdsn import Client
# party2, stream2 = tribe.client_detect(
#         client=Client('IRIS'), starttime=start_time, endtime=end_time, threshold=20.,
#         threshold_type="MAD", trig_int=30.0, plot=False, return_stream=True)
# party2_out = 'party2'
# party2_outpath = '/home/ptan/attempt_eqcorrscan/output/' + party2_out
# if os.path.exists(party2_outpath+'.tgz'):
#     os.remove(party2_outpath+'.tgz')
#     party2.write(party2_outpath + '.tgz' , format='tar')
# else:
#     party2.write(party2_outpath + '.tgz' , format='tar')

