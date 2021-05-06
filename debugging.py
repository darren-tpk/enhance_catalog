
#### DEBUGGING
# inputs
from eqcorrscan.core.match_filter.matched_filter import _group_process, match_filter
from eqcorrscan.core.match_filter.template import group_templates
from eqcorrscan.core.match_filter.detection import Detection
from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.utils.plotting import _match_filter_plot
import logging
from timeit import default_timer

import numpy as np
from obspy import Catalog, UTCDateTime, Stream

from eqcorrscan.core.match_filter.helpers import (
    _spike_test, extract_from_stream)

from eqcorrscan.utils.correlate import get_stream_xcorr, _get_registerd_func, XCOR_FUNCS, _get_array_dicts
from eqcorrscan.utils.findpeaks import multi_find_peaks
from eqcorrscan.utils.pre_processing import (
    dayproc, shortproc, _prep_data_for_correlation)


import logging

from pyfftw.pyfftw import next_fast_len

logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

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
group = template_groups[0]
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
n_groups = 1

st_chunk = streams[0]
chunk_start, chunk_end = (min(tr.stats.starttime for tr in st_chunk),
                          max(tr.stats.endtime for tr in st_chunk))
st_chunk.trim(starttime=chunk_start, endtime=chunk_end)
for tr in st_chunk:
    if len(tr) > len(st_chunk[0]):
        tr.data = tr.data[0:len(st_chunk[0])]
i=0
end_group = len(templates)
start_group = 0
template_group = [t for t in templates[start_group: end_group]]
### MATCH FILTER HERE
# detections += match_filter(
#     template_names=[t.name for t in template_group],
#     template_list=[t.st for t in template_group], st=st_chunk,
#     xcorr_func=xcorr_func, concurrency=concurrency,
#     threshold=threshold, threshold_type=threshold_type,
#     trig_int=trig_int, plot=plot, plotdir=plotdir, cores=cores,
#     full_peaks=full_peaks, peak_cores=process_cores,
#     )
template_names=[t.name for t in template_group]
template_list=[t.st for t in template_group]
st=st_chunk
xcorr_func=xcorr_func
concurrency=concurrency,
threshold=threshold
threshold_type=threshold_type
trig_int=trig_int
plot=plot
plotdir=plotdir
cores=cores
full_peaks=full_peaks
peak_cores=process_cores



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

detections += match_filter(
    template_names=[t.name for t in template_group],
    template_list=[t.st for t in template_group], st=st_chunk,
    xcorr_func=xcorr_func, concurrency=concurrency,
    threshold=threshold, threshold_type=threshold_type,
    trig_int=trig_int, plot=plot, plotdir=plotdir, cores=cores,
    full_peaks=full_peaks, peak_cores=process_cores,
    **kwargs)

template_names=[t.name for t in template_group]
template_list=[t.st for t in template_group]
st=st_chunk
xcorr_func=xcorr_func
concurrency=concurrency
threshold=threshold
threshold_type=threshold_type
trig_int=trig_int
plot=plot
plotdir=plotdir
cores=cores
full_peaks=full_peaks
peak_cores=process_cores

plot_format='png'
output_cat=False
output_event=True
extract_detections=False
arg_check=True
spike_test=True
copy_data=True,
export_cccsums=False

if arg_check:
    # Check the arguments to be nice - if arguments wrong type the parallel
    # output for the error won't be useful
    if not isinstance(template_names, list):
        raise MatchFilterError('template_names must be of type: list')
    if not isinstance(template_list, list):
        raise MatchFilterError('templates must be of type: list')
    if not len(template_list) == len(template_names):
        raise MatchFilterError('Not the same number of templates as names')
    for template in template_list:
        if not isinstance(template, Stream):
            msg = 'template in template_list must be of type: ' + \
                  'obspy.core.stream.Stream'
            raise MatchFilterError(msg)
    if not isinstance(st, Stream):
        msg = 'st must be of type: obspy.core.stream.Stream'
        raise MatchFilterError(msg)
    if str(threshold_type) not in [str('MAD'), str('absolute'),
                                   str('av_chan_corr')]:
        msg = 'threshold_type must be one of: MAD, absolute, av_chan_corr'
        raise MatchFilterError(msg)
    for tr in st:
        if not tr.stats.sampling_rate == st[0].stats.sampling_rate:
            raise MatchFilterError('Sampling rates are not equal %f: %f' %
                                   (tr.stats.sampling_rate,
                                    st[0].stats.sampling_rate))
    for template in template_list:
        for tr in template:
            if not tr.stats.sampling_rate == st[0].stats.sampling_rate:
                raise MatchFilterError(
                    'Template sampling rate does not '
                    'match continuous data')
    for template in template_list:
        for tr in template:
            if isinstance(tr.data, np.ma.core.MaskedArray):
                raise MatchFilterError(
                    'Template contains masked array, split first')
if spike_test:
    _spike_test(st)
if copy_data:
    # Copy the stream here because we will muck about with it
    stream = st.copy()
    templates = [t.copy() for t in template_list]
    _template_names = template_names.copy()  # This can be a shallow copy

stream, templates, _template_names = _prep_data_for_correlation(
    stream=stream, templates=templates, template_names=_template_names)

#multichannel_normxcorr = get_stream_xcorr(xcorr_func, concurrency)
name_or_func=xcorr_func
func = XCOR_FUNCS[name_or_func or 'default']
concur = concurrency or 'stream_xcorr'
multichannel_normxcorr = getattr(func, concur) # _fftw_stream_xcorr

outtic = default_timer()

#[cccsums, no_chans, chans] = multichannel_normxcorr(
    templates=templates, stream=stream, cores=cores)
templates = templates
stream = stream
cores = cores
stack = True
num_cores_inner = None
from multiprocessing import Pool as ProcessPool, cpu_count
if num_cores_inner is None:
    num_cores_inner = int(os.getenv("OMP_NUM_THREADS", cpu_count()))
chans = [[] for _i in range(len(templates))]
array_dict_tuple = _get_array_dicts(templates, stream, stack=stack)
stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple

# cccsums, tr_chans = fftw_multi_normxcorr(
# template_array = template_dict, stream_array = stream_dict,
# pad_array = pad_dict, seed_ids = seed_ids, cores_inner = num_cores_inner,
# stack = stack, *args, ** kwargs)
template_array = template_dict
stream_array = stream_dict
pad_array = pad_dict
seed_ids = seed_ids
cores_inner = num_cores_inner
stack = stack
utilslib = _load_cdll('libutils')
from future.utils import native_str
import ctypes
utilslib.multi_normxcorr_fftw.argtypes = [
np.ctypeslib.ndpointer(dtype=np.float32,
                       flags=native_str('C_CONTIGUOUS')),
ctypes.c_long, ctypes.c_long, ctypes.c_long,
np.ctypeslib.ndpointer(dtype=np.float32,
                       flags=native_str('C_CONTIGUOUS')),
ctypes.c_long,
np.ctypeslib.ndpointer(dtype=np.float32,
                       flags=native_str('C_CONTIGUOUS')),
ctypes.c_long,
np.ctypeslib.ndpointer(dtype=np.intc,
                       flags=native_str('C_CONTIGUOUS')),
np.ctypeslib.ndpointer(dtype=np.intc,
                       flags=native_str('C_CONTIGUOUS')),
ctypes.c_int,
np.ctypeslib.ndpointer(dtype=np.intc,
                       flags=native_str('C_CONTIGUOUS')),
np.ctypeslib.ndpointer(dtype=np.intc,
                       flags=native_str('C_CONTIGUOUS')),
ctypes.c_int]
utilslib.multi_normxcorr_fftw.restype = ctypes.c_int
###############
# pre processing
used_chans = []
template_len = template_array[seed_ids[0]].shape[1]
for seed_id in seed_ids:
    used_chans.append(~np.isnan(template_array[seed_id]).any(axis=1))
template_array[seed_id] = (
(template_array[seed_id] -
 template_array[seed_id].mean(axis=-1, keepdims=True)) / (
        template_array[seed_id].std(axis=-1, keepdims=True) *
        template_len))
template_array[seed_id] = np.nan_to_num(template_array[seed_id])
n_channels = len(seed_ids)
n_templates = template_array[seed_ids[0]].shape[0]
image_len = stream_array[seed_ids[0]].shape[0]
# In testing, 2**13 consistently comes out fastest - setting to
# default. https://github.com/eqcorrscan/EQcorrscan/pull/285
fft_len = min(2 ** 13, next_fast_len(template_len + image_len - 1))
template_array = np.ascontiguousarray(
[template_array[x] for x in seed_ids], dtype = np.float32)
multipliers = {}
for x in seed_ids:
# Check that stream is non-zero and above variance threshold
    if not np.all(stream_array[x] == 0) and np.var(stream_array[x]) < 1e-8:
        # Apply gain
        stream_array[x] *= MULTIPLIER
        Logger.warning(f"Low variance found for {x}, applying gain "
        "to stabilise correlations")
        multipliers.update({x: MULTIPLIER})
    else:
        multipliers.update({x: 1})
stream_array = np.ascontiguousarray([stream_array[x] for x in seed_ids], dtype = np.float32)

ccc_length = image_len - template_len + 1

if stack:
    cccs = np.zeros((n_templates, ccc_length), np.float32)
else:
    cccs = np.zeros((n_templates, n_channels, ccc_length), dtype = np.float32)

used_chans_np = np.ascontiguousarray(used_chans, dtype=np.intc)
pad_array_np = np.ascontiguousarray([pad_array[seed_id] for seed_id in seed_ids], dtype = np.intc)
variance_warnings = np.ascontiguousarray(np.zeros(n_channels), dtype = np.intc)
missed_correlations = np.ascontiguousarray(np.zeros(n_channels), dtype = np.intc)

# call C function
ret = utilslib.multi_normxcorr_fftw(template_array, n_templates, template_len, n_channels, stream_array,
                                    image_len, cccs, fft_len, used_chans_np, pad_array_np,
                                    cores_inner, variance_warnings, missed_correlations, int(stack))
if ret < 0:
    raise MemoryError("Memory allocation failed in correlation C-code")
elif ret > 0:
Logger.critical(
    'Out-of-range correlation in C-code, see WARNING from C-code.'
    'You are STRONGLY RECOMMENDED to check your data for spikes, '
    'clipping or non-physical artifacts')
# raise CorrelationError("Internal correlation error")
for i, missed_corr in enumerate(missed_correlations):
    if missed_corr:
        Logger.debug(
            f"{missed_corr} correlations not computed on {seed_ids[i]}, "
            f"are there gaps in the data? If not, consider "
            "increasing gain")
for i, variance_warning in enumerate(variance_warnings):
    if variance_warning and variance_warning > template_len:
        Logger.warning(
            f"Low variance found in {variance_warning} places for "
            f"{seed_ids[i]}, check result.")
# Remove gain
for i, x in enumerate(seed_ids):
    stream_array[i] *= multipliers[x]
cccsums, tr_chans = cccs, used_chans


no_chans = np.sum(np.array(tr_chans).astype(np.int), axis=0)
for seed_id, tr_chan in zip(seed_ids, tr_chans):
    for
chan, state in zip(chans, tr_chan):
if state:
    chan.append((seed_id.split('.')[1],
                 seed_id.split('.')[-1].split('_')[0]))
return cccsums, no_chans, chans
########################

outtoc = default_timer()

detections = []
if output_cat:
    det_cat = Catalog()
if str(threshold_type) == str("absolute"):
    thresholds = [threshold for _ in range(len(cccsums))]
elif str(threshold_type) == str('MAD'):
    thresholds = [threshold * np.median(np.abs(cccsum))
                  for cccsum in cccsums]
else:
    thresholds = [threshold * no_chans[i] for i in range(len(cccsums))]
if peak_cores is None:
    peak_cores = cores
outtic = default_timer()
all_peaks = multi_find_peaks(
    arr=cccsums, thresh=thresholds, parallel=parallel,
    trig_int=int(trig_int * stream[0].stats.sampling_rate),
    full_peaks=full_peaks, cores=peak_cores)
outtoc = default_timer()
Logger.info("Finding peaks took {0:.4f}s".format(outtoc - outtic))
for i, cccsum in enumerate(cccsums):
    if export_cccsums:
        fname = (f"{_template_names[i]}-{stream[0].stats.starttime}-"
                 f"{stream[0].stats.endtime}_cccsum.npy")
        np.save(file=fname, arr=cccsum)
        Logger.info(f"Saved correlation statistic to {fname}")
    if np.abs(np.mean(cccsum)) > 0.05:
        Logger.warning('Mean is not zero!  Check this!')
    # Set up a trace object for the cccsum as this is easier to plot and
    # maintains timing
    if plot:
        _match_filter_plot(
            stream=stream, cccsum=cccsum, template_names=_template_names,
            rawthresh=thresholds[i], plotdir=plotdir,
            plot_format=plot_format, i=i)
    if all_peaks[i]:
        Logger.debug("Found {0} peaks for template {1}".format(
            len(all_peaks[i]), _template_names[i]))
        for peak in all_peaks[i]:
            detecttime = (
                stream[0].stats.starttime +
                peak[1] / stream[0].stats.sampling_rate)
            detection = Detection(
                template_name=_template_names[i], detect_time=detecttime,
                no_chans=no_chans[i], detect_val=peak[0],
                threshold=thresholds[i], typeofdet='corr', chans=chans[i],
                threshold_type=threshold_type, threshold_input=threshold)
            if output_cat or output_event:
                detection._calculate_event(template_st=templates[i])
            detections.append(detection)
            if output_cat:
                det_cat.append(detection.event)
    else:
        Logger.debug("Found 0 peaks for template {0}".format(
            _template_names[i]))
Logger.info("Made {0} detections from {1} templates".format(
    len(detections), len(templates)))
if extract_detections:
    detection_streams = extract_from_stream(stream, detections)
del stream, templates

if output_cat and not extract_detections:
    return detections, det_cat
elif not extract_detections:
    return detections
elif extract_detections and not output_cat:
    return detections, detection_streams
else:
    return detections, det_cat, detection_streams