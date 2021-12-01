import glob
import obspy
from obspy.io.mseed import InternalMSEEDError

data_dir = '/home/data/augustine/'
file_list = glob.glob(data_dir + '*')

for i, f in enumerate(file_list):
    st_header = obspy.read(f, headonly=True)
    if st_header[0].stats.network == 'AK':
        try:
            st = obspy.read(f)
            for tr in st:
                tr.stats.network = 'AV'
            st.write(f, format="MSEED")
        except InternalMSEEDError:
            print(str(i) + ' failed. InternalMSEEDError')
    if (i / 100).is_integer():
        print(i)

[3047,7521,7793,11465]
'/home/data/augustine/AUI.EHZ.2006:135:00:00:00'
'/home/data/augustine/AUE.BDF.2006:244:00:00:00'
'/home/data/augustine/AUS.EHZ.2006:244:00:00:00'
'/home/data/augustine/AUE.BDL.2006:244:00:00:00'




# [prepare_catalog_stream] function to read streams pertaining to an input catalog
def prepare_catalog_stream(data_dir,catalog,resampling_frequency,tolerance,max_zeros=100,npts_threshold=100):

    # Import dependencies
    import glob
    from obspy import Stream, read
    from toolbox import remove_boxcars, remove_bad_traces

    # Initialize list for all station channel data needed for the day's events
    data_filenames = []

    # Loop through catalog's events
    for event in catalog:

        # Loop through event's picks
        for pick in event.picks:

            # Extract stachan, year and julday of pick
            sta = pick.waveform_id.station_code
            chan = pick.waveform_id.channel_code
            pick_year = pick.time.year
            pick_julday = pick.time.julday

            # Craft file string and append
            data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{pick_julday:03}' + ':*'
            data_filenames.append(data_filename)

            # Add next day if pick occurs in the first 15 minutes of the day
            if pick.time.hour == 0 and pick.time.minute < 15:
                if pick.time.julday == 1:  # special case for first day of the year
                    data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year - 1) + ':365:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_dir + sta + '.' + chan + '.' + str(
                        pick_year) + ':' + f'{(pick_julday - 1):03}' + ':*'
                    data_filenames.append(data_filename)

            # Add previous day if pick occurs in the last 15 minutes of the day
            if pick.time.hour == 23 and pick.time.minute > 45:
                if pick.time.julday == 365:  # special case for last day of the year
                    data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year + 1) + ':001:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_dir + sta + '.' + chan + '.' + str(
                        pick_year) + ':' + f'{(pick_julday + 1):03}' + ':*'
                    data_filenames.append(data_filename)

    # Make list unique, then sort
    data_filenames = list(set(data_filenames))
    data_filenames.sort()

    # Read in all required traces, merging into a stream object
    stream = Stream()
    for data_filename in data_filenames:
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                stream_contribution = read(matching_filename)
                stream = stream + stream_contribution
            except:
                continue

    # Remove bad traces:
    stream = remove_bad_traces(stream,max_zeros=max_zeros,npts_threshold=npts_threshold)

    # Remove boxcar and spikes, resample and merge
    stream = remove_boxcars(stream, tolerance)
    if resampling_frequency:
        stream = stream.resample(resampling_frequency)
    else:
        for tr in stream:
            sr = tr.stats.sampling_rate
            if abs(100-sr) < abs(50-sr):
                tr.resample(100)
            else:
                tr.resample(50)
    stream = stream.merge()

    # Return stream object
    return stream



