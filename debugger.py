# [calculate_relative_magnitudes] function to calculate magnitudes using CC based on Schaff & Richards (2014)
# Uses template events in PEC as reference magnitudes for relocatable catalog
# This function rethresholds the detected catalog by min_cc while processing it
def calculate_relative_magnitudes(catalog, tribe, data_dir, noise_window, signal_window, min_cc, min_snr,
                                  shift_len, tolerance, resampling_frequency, lowcut, highcut, filt_order,
                                  use_s_picks=False, verbose=False):

    # Import dependencies
    import numpy as np
    from itertools import compress
    from eqcorrscan.utils.mag_calc import relative_magnitude
    from toolbox import prepare_catalog_stream
    from obspy import Catalog, UTCDateTime, Stream
    from obspy.core.event import Magnitude

    # Derive template names from tribe
    template_names = [template.name for template in tribe]

    # Determine all detection-template date combinations
    catalog_times = []
    template_times = []
    for event in catalog:
        # Obtain template event
        template_index = template_names.index(event.comments[0].text.split(' ')[1])
        # Append catalog and template times
        catalog_times.append(event.origins[0].time)
        template_time = np.min([p.time for p in tribe[template_index].event.picks]) - signal_window[0]
        template_times.append(template_time)
    # Get unique list of date pairs
    date_pairs = list(zip([t.date for t in catalog_times], [t.date for t in template_times]))
    date_pairs_unique = list(dict.fromkeys(date_pairs))

    # Create a boolean list to keep track of which events have been attempted
    tracker = np.array([False for i in range(len(catalog))])

    # Print progress count for the number of successfully calculated magnitudes
    if verbose:
        print('Count of successful magnitude computations:')
        count = 0

    # Loop through events in catalog
    for i, event in enumerate(catalog):

        # Check if event has had its magnitude calculation attempted. If True, skip to next
        if tracker[i]:
            continue
        # Otherwise check for other events occuring on the same day so that we load the local data at one go
        else:
            date_pair = date_pairs[i]
            sub_catalog_bool = [dp == date_pair for dp in date_pairs]
            sub_catalog_index = np.flatnonzero(sub_catalog_bool)
            tracker[sub_catalog_index] = True

        # Obtain day-long streams for the sub-catalog
        sub_catalog = Catalog(list(compress(catalog, sub_catalog_bool)))
        master_catalog_stream = prepare_catalog_stream(data_dir, sub_catalog, resampling_frequency, tolerance)
        master_catalog_stream = master_catalog_stream.trim(starttime=UTCDateTime(date_pair[0]),
                                                           endtime=UTCDateTime(date_pair[0]) + 86400, pad=True)

        # Craft catalog from their templates and obtain their corresponding stream
        sub_templates_list = []
        for event in sub_catalog:
            template_index = template_names.index(event.comments[0].text.split(' ')[1])
            sub_templates_list.append(tribe[template_index].event)
        sub_templates = Catalog(sub_templates_list)
        master_template_stream = prepare_catalog_stream(data_dir, sub_templates, resampling_frequency, tolerance)
        master_template_stream = master_template_stream.trim(starttime=UTCDateTime(date_pair[1]),
                                                             endtime=UTCDateTime(date_pair[1]) + 86400, pad=True)

        # Now loop through sub_catalog to compute magnitudes
        for j, event in enumerate(sub_catalog):

            # Do a check on detection threshold. Skip events that record an av_chan_corr < min_cc
            # This makes sure that there are enough picks to compute relative magnitude
            av_chan_corr = float(event.comments[2].text.split('=')[1]) / event.comments[3].text.count(')')
            if abs(av_chan_corr) < min_cc:
                if verbose:
                    print('Skipping event %s. abs(av_chan_corr) is %.2f, lower than min_cc %.2f.' % (
                    event.resource_id.id, abs(av_chan_corr), min_cc))
                continue

            # Obtain template event
            template = sub_templates[j]

            # If the detection value is very high, and the two events are close in time, it is a self-detection.
            # Give target event the same mag as template
            template_time = np.min([p.time for p in template.picks]) - signal_window[0]
            if abs(av_chan_corr) > 0.9 and abs(event.origins[0].time - template_time) < 0.1:
                try:
                    original_mag = template.preferred_magnitude().mag
                except:
                    original_mag = template.magnitudes[0].mag
                catalog[sub_catalog_index[j]].magnitudes.append(Magnitude(mag=original_mag))
                if verbose:
                    count += 1
                    print('%d (Index [i,j] = [%d,%d])' % (count, i, j))
                continue

            # Derive a subset of the master stream that correspond to picked stachans
            event_netstachans = [
                (pick.waveform_id.network_code, pick.waveform_id.station_code, pick.waveform_id.channel_code) for pick
                in event.picks]
            stream_boolean = [(trace.stats.network, trace.stats.station, trace.stats.channel) in event_netstachans for
                              trace in master_catalog_stream]
            sub_catalog_stream = Stream(compress(master_catalog_stream, stream_boolean)).copy()

            # Now trim the sub catalog stream based on pick times and execute fft
            for trace in sub_catalog_stream:
                trace_netstachan = (trace.stats.network, trace.stats.station, trace.stats.channel)
                pick_index = event_netstachans.index(trace_netstachan)
                reference_pick = event.picks[pick_index]
                length = signal_window[1] - signal_window[0]
                trace.trim(reference_pick.time + noise_window[0] - shift_len - 0.05 * length,
                           reference_pick.time + signal_window[1] + shift_len + 0.05 * length)
                # Do some stream checks before detrending and filtering
                if (type(trace.data) == np.ndarray or sum(trace.data.mask) == 0) and len(trace.data) > (
                        length * trace.stats.sampling_rate):
                    trace.detrend('simple')
                    trace.taper(max_percentage=None, max_length=0.05 * length)
                    trace.filter(type='bandpass', freqmin=lowcut, freqmax=highcut, corners=filt_order, zerophase=True)
                    trace.trim(reference_pick.time + noise_window[0], reference_pick.time + signal_window[1])
                else:
                    sub_catalog_stream.remove(trace)

            # Derive a subset of the template stream that corresponds to template picked stachans
            template_netstachans = [
                (pick.waveform_id.network_code, pick.waveform_id.station_code, pick.waveform_id.channel_code) for pick
                in template.picks]
            stream_boolean = [(trace.stats.network, trace.stats.station, trace.stats.channel) in template_netstachans for
                              trace in master_template_stream]
            sub_template_stream = Stream(compress(master_template_stream, stream_boolean)).copy()

            # Now trim the sub template stream based on pick times and execute fft
            for trace in sub_template_stream:
                trace_netstachan = (trace.stats.network, trace.stats.station, trace.stats.channel)
                pick_index = template_netstachans.index(trace_netstachan)
                reference_pick = template.picks[pick_index]
                trace.trim(reference_pick.time + noise_window[0] - shift_len - 0.05 * length,
                           reference_pick.time + signal_window[1] + shift_len + 0.05 * length)
                # Do some stream checks before detrending and filtering
                if (type(trace.data) == np.ndarray or sum(trace.data.mask) == 0) and len(trace.data) > (
                        length * trace.stats.sampling_rate):
                    trace.detrend('simple')
                    trace.taper(max_percentage=None, max_length=0.05 * length)
                    trace.filter(type='bandpass', freqmin=lowcut, freqmax=highcut, corners=filt_order, zerophase=True)
                    trace.trim(reference_pick.time + noise_window[0], reference_pick.time + signal_window[1])
                else:
                    sub_template_stream.remove(trace)

            # Calculate magnitude differences determined by each channel
            delta_mags, ccs = relative_magnitude(sub_catalog_stream, sub_template_stream, event, template,
                                                 noise_window=noise_window, signal_window=signal_window,
                                                 min_snr=0, min_cc=min_cc, use_s_picks=False, correlations=None,
                                                 shift=shift_len, return_correlations=True, correct_mag_bias=True)

            # If the SNR window fails to record above min_snr, we skip the event
            # It is likely that there is a larger event preceding it in the noise window
            # The magnitude estimation will therefore be inaccurate
            if len(delta_mags) == 0:
                if verbose:
                    print('Skipping event %s. No channels record above min_snr.' % (event.resource_id.id))
                continue

            # Now calculate event's magnitude
            # Weight each channel's delta_mag by cc so that well correlated channels hold higher weight in mag change
            delta_mag_values = [delta_mag[1] for delta_mag in delta_mags.items() if not np.isnan(delta_mag[1])]
            if len(delta_mag_values) == 0:
                if verbose:
                    print('Skipping event %s. All calculated delta magnitudes are nan.' % (event.resource_id.id))
                continue
            cc_values = [cc[1] for cc in ccs.items() if (cc[1] > min_cc) and not np.isnan(cc[1])]
            try:
                original_mag = template.preferred_magnitude().mag
            except:
                original_mag = template.magnitudes[0].mag
            target_mag = original_mag + np.dot(delta_mag_values, cc_values) / np.sum(cc_values)

            # Add magnitude to target event
            catalog[sub_catalog_index[j]].magnitudes.append(Magnitude(mag=target_mag))
            if verbose:
                count += 1
                print('%d (Index [i,j] = [%d,%d])' % (count, i, j))

    # Return catalog
    return catalog


##########################

# from obspy import UTCDateTime
# from toolbox import read_hypoi, writer
#
# print('REDOUBT:')
#
# # def read_hypoi(hypoi_file,...):
# hypoi_file = '/home/ptan/enhance_catalog/data/avo/avo_1989-2018_hypoi.txt'
# VOLC_LAT = 60.490  # dome
# VOLC_LON = -152.7626  # dome
# time_lim = (UTCDateTime(2008,4,1),UTCDateTime(2009,9,1))
# radial_lim = (VOLC_LAT,VOLC_LON,25)
# rmse_max = 0.30
# gap_max = 210
# num_P_min = 3
# num_S_min = 2
# channel_convention = True
# summary=True
#
# redoubt = read_hypoi(hypoi_file,
#                      time_lim=time_lim,
#                      radial_lim=radial_lim,
#                      rmse_max=rmse_max,
#                      gap_max=gap_max,
#                      num_P_min=num_P_min,
#                      num_S_min=num_S_min,
#                      channel_convention=True,
#                      summary=True)
#
# writer(outpath+'redoubt_20080401_20090901.xml',redoubt)
#
# print('AUGUSTINE:')
#
# VOLC_LAT = 59.363
# VOLC_LON = -153.43
# time_lim = (UTCDateTime(2005,4,1),UTCDateTime(2006,5,1))
# radial_lim = (VOLC_LAT,VOLC_LON,25)
# rmse_max = 0.30
# gap_max = 190
# num_P_min = 3
# num_S_min = 2
# channel_convention = True
# summary=True
#
# augustine = read_hypoi(hypoi_file,
#                      time_lim=time_lim,
#                      radial_lim=radial_lim,
#                      rmse_max=rmse_max,
#                      gap_max=gap_max,
#                      num_P_min=num_P_min,
#                      num_S_min=num_S_min,
#                      channel_convention=True,
#                      summary=True)
#
# writer(outpath+'augustine_20050401_20060501.xml',augustine)
#
# # #############
# #'
# main_dir = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/redpy_results/'
# max_1 = 40
# max_2 = 81
# max_3 = 18
# max_4 = 115
# col_a = []
# col_b = []
# file_type = 'catalog.txt'
# subfolders = ['augustine_1','augustine_2','augustine_3','augustine_4']
# for subfolder in subfolders:
#     redpy_results_dir = main_dir + 'augustine5/' + subfolder + '/'
#     with open(redpy_results_dir + file_type) as f:
#         if file_type == 'cores.txt' or file_type == 'catalog.txt':
#             lines = f.readlines()
#             subcol_a = [int(line.split(' ')[0]) for line in lines]
#             if subfolder == 'augustine_2':
#                 subcol_a = [v+1+max_1 for v in subcol_a]
#             elif subfolder == 'augustine_3':
#                 subcol_a = [v+2+max_1+max_2 for v in subcol_a]
#             elif subfolder == 'augustine_4':
#                 subcol_a = [v+3+max_1+max_2+max_3 for v in subcol_a]
#             col_a = col_a + subcol_a
#             subcol_b = [line.split(' ')[1][:-1] for line in lines]
#             col_b = col_b + subcol_b
#         elif file_type == 'orphancatalog.txt':
#             lines = f.readlines()
#             subcol_a = [line[:-1] for line in lines]
#             col_a = col_a + subcol_a
#         f.close()
# with open(r'/Users/darrentpk/Desktop/GitHub/enhance_catalog/redpy_results/augustine5/catalog.txt', 'w') as fp:
#     for i in range(len(col_a)):
#         # write each item on a new line
#         if file_type == 'cores.txt' or file_type == 'catalog.txt':
#             fp.write("%d %s\n" % (col_a[i],col_b[i]))
#         elif file_type == 'orphancatalog.txt':
#             fp.write("%s\n" % (col_a[i]))
#     fp.close()
#     print('Done')
# #
# # #########

