# scan data to make detections

# import packages that we need
import os
import glob
import numpy as np
from eqcorrscan import Party
from eqcorrscan.core.match_filter.tribe import read_tribe
from obspy import UTCDateTime, Stream, read
# import toolbox functions
from toolbox import remove_boxcars, reader, writer

# define all variables here
start_time = UTCDateTime(2009, 2, 26, 0, 0, 0)
end_time = start_time + (24 * 60 * 60)
threshold = 30  # Used to be 20 without RDWB and RDJH
threshold_type = 'MAD'
trig_int = 30
parallel_process = 'True'
tolerance = 5e4 # for boxcar removal

# read tribe data
tribe = reader('/home/ptan/attempt_eqcorrscan/output/tribe.tgz')

# print station numbers of each template before filter
print('\nBefore filter:')
print(tribe)
# for template in tribe:
#     num_stations = len({tr.stats.station for tr in template.st})
#     print(template.name + ': ' + str(num_stations) + ' stations')
# filter away templates with < 3 stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 3]
print('\nAfter removing templates with < 3 stations...')
print(tribe)
# for template in tribe:
#     num_stations = len({tr.stats.station for tr in template.st})
#     print(template.name + ': ' + str(num_stations) + ' stations')

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

# BRUTE FORCE SCAN: load each day's data and scan day by day
party_all = Party()
num_days = int(np.floor((end_time - start_time) / 86400))
for i in range(num_days):
    # define boundaries of our day's scan
    t1 = start_time + (i * 86400)
    t2 = start_time + ((i + 1) * 86400)
    # initialize and download stream
    stream = Stream()
    for data_filehead in data_fileheads:
        data_filename = data_filehead + str(t1.year) + ':' + f'{t1.julday :03}' + ':*'
        real_data_filename = (glob.glob(data_filename))
        if real_data_filename:
            stream_contribution = read(real_data_filename.pop())
            # for trace_contribution in stream_contribution:
            #    median_filter(trace_contribution,windowlength=despike_window,interp_len=despike_interp)
            stream = stream + stream_contribution
    stream = remove_boxcars(stream,tolerance)
    stream.detrend("simple").merge()
    # scan the current day
    party = tribe.detect(
        stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, overlap=None, parallel_process='True')
    # append party to party_all
    party_all = party_all + party

# write the combined party object as a tar file
party_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_outpath+'party.tgz')

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

