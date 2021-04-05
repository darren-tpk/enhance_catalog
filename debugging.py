# scan data to make detections

# import packages that we need
import os
import glob
import numpy as np
from eqcorrscan import Party
from eqcorrscan.core.match_filter.matched_filter import _group_detect, _group_process, match_filter
from eqcorrscan.core.match_filter.template import group_templates
from eqcorrscan.core.match_filter.tribe import read_tribe
from obspy import UTCDateTime, Stream, read
# import toolbox functions
from toolbox import remove_boxcars, reader, writer

# define all variables here
start_time = UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 5, 1, 0, 0, 0)
threshold = 30  # Used to be 20 without RDWB and RDJH
threshold_type = 'MAD'
trig_int = 30
parallel_process = 'True'
tolerance = 4e4 # for boxcar removal

# read tribe data
tribe = reader('/home/ptan/attempt_eqcorrscan/output/new_tribe.tgz')

# filter away templates with < 3 stations
print('\nBefore filter:')
print(tribe)
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 3]
print('\nAfter removing templates with < 3 stations...')
print(tribe)

# # filter away templates with < 3 picks
# print('\nBefore filter:')
# print(tribe)
# tribe.templates = [t for t in tribe if len(t.event.picks) >= 3]
# print('\nAfter removing templates with < 3 picks...')
# print(tribe)

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
    print('Now at %s...' % str(t1))
    # initialize and download stream
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
    stream = remove_boxcars(stream,tolerance)
    stream = stream.detrend("simple")
    stream = stream.merge()
    # scan the current day
    try:
        party = tribe.detect(
            stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, overlap=None, parallel_process='False')
        # append party to party_all
        party_all = party_all + party
        print('Party created, appending.')
    except:
        print('Party failed, skipping.')
        continue

#### DEBUGGING
# inputs
stream = stream
threshold = threshold
threshold_type = threshold_type
trig_int = trig_int
plot = False
plotdir = None
daylong = True
parallel_process = False
xcorr_func = None
concurrency = None
cores = None
ignore_length = False
ignore_bad_data = False
group_size = None
overlap = None
full_peaks = False
save_progress = False
process_cores = None

party = Party()
template_groups = group_templates(tribe.templates)
# now we can compute the detections for each group
for group in template_groups:
    templates=group
    stream=stream.copy()
    threshold=threshold
    threshold_type=threshold_type
    trig_int=trig_int
    plot=plot
    group_size=group_size
    pre_processed=False
    daylong=daylong
    parallel_process=parallel_process
    xcorr_func=xcorr_func
    concurrency=concurrency
    cores=cores
    ignore_length=ignore_length
    overlap=overlap
    plotdir=plotdir,
    full_peaks=full_peaks
    process_cores=process_cores,
    ignore_bad_data=ignore_bad_data
    arg_check=False

    master = templates[0]
    # Check that they are all processed the same.
    lap = 0.0
    for template in templates:
        starts = [t.stats.starttime for t in template.st.sort(['starttime'])]
        if starts[-1] - starts[0] > lap:
            lap = starts[-1] - starts[0]
        if not template.same_processing(master):
            raise MatchFilterError('Templates must be processed the same.')
    if overlap is None:
        overlap = 0.0
    elif not isinstance(overlap, float) and str(overlap) == str("calculate"):
        overlap = lap
    elif not isinstance(overlap, float):
        raise NotImplementedError(
            "%s is not a recognised overlap type" % str(overlap))
    if overlap >= master.process_length:
        Logger.warning(
            f"Overlap of {overlap} s is greater than process "
            f"length ({master.process_length} s), ignoring overlap")
        overlap = 0
    if not pre_processed:
        if process_cores is None:
            process_cores = cores
        streams = _group_process(
            template_group=templates, parallel=parallel_process,
            cores=process_cores, stream=stream, daylong=daylong,
            ignore_length=ignore_length, ignore_bad_data=ignore_bad_data,
            overlap=overlap)
        # for _st in streams:
        #     Logger.debug(f"Processed stream:\n{_st.__str__(extended=True)}")
    else:
        Logger.warning('Not performing any processing on the continuous data.')
        streams = [stream]
    detections = []

    party = Party()
    if group_size is not None:
        n_groups = int(len(templates) / group_size)
        if n_groups * group_size < len(templates):
            n_groups += 1
    else:
        n_groups = 1

    for st_chunk in streams:
        chunk_start, chunk_end = (min(tr.stats.starttime for tr in st_chunk),
                                  max(tr.stats.endtime for tr in st_chunk))
        st_chunk.trim(starttime=chunk_start, endtime=chunk_end)
        for tr in st_chunk:
            if len(tr) > len(st_chunk[0]):
                tr.data = tr.data[0:len(st_chunk[0])]
        for i in range(n_groups):
            if group_size is not None:
                end_group = (i + 1) * group_size
                start_group = i * group_size
                if i == n_groups:
                    end_group = len(templates)
            else:
                end_group = len(templates)
                start_group = 0
            template_group = [t for t in templates[start_group: end_group]]
            detections += match_filter(
                template_names=[t.name for t in template_group],
                template_list=[t.st for t in template_group], st=st_chunk,
                xcorr_func=xcorr_func, concurrency=concurrency,
                threshold=threshold, threshold_type=threshold_type,
                trig_int=trig_int, plot=plot, plotdir=plotdir, cores=cores,
                full_peaks=full_peaks, peak_cores=process_cores,
                **kwargs)
            for template in template_group:
                family = Family(template=template, detections=[])
                for detection in detections:
                    if detection.template_name == template.name:
                        for pick in detection.event.picks:
                            pick.time += template.prepick
                        for origin in detection.event.origins:
                            origin.time += template.prepick
                        family.detections.append(detection)
                party += family

    party += group_party

if len(party) > 0:
    for family in party:
        if family is not None:
            family.detections = family._uniq().detections
return party