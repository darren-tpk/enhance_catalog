#%% Toolbox of all functions

#%% [remove_boxcars] function to remove artificial boxcar-like signals from stream
def remove_boxcars(st,tolerance):

    # Import dependencies
    import numpy as np
    from obspy import Stream

    # Initialize output stream
    st_out = Stream()

    # Loop through traces in stream
    for tr in st:

        # Detrend data and find spikes above tolerance (median(abs(data)) multiplier)
        data = tr.data
        tr_detrended = tr.copy().detrend()
        data_detrended = tr_detrended.data
        spike_indices = np.where(data_detrended > tolerance * (np.median(abs(data_detrended))+1))[0]

        # Initialize isolated spikes list
        isolated_spike_index = []

        # Loop through all spike indices to find isolated spikes
        for i, spike_index in enumerate(spike_indices):

            # If spike indices is 1x1
            if len(spike_indices) == 1:
                isolated_spike_index.append(i)

            # If isolated spike is the first element
            elif i == 0 and (spike_index != spike_indices[i+1]-1):
                isolated_spike_index.append(i)

            # If isolated spike is the last element
            elif i == len(spike_indices)-1 and (spike_index != spike_indices[i-1]+1):
                isolated_spike_index.append(i)

            # If the isolated spike is in the middle of the list
            elif (spike_index != spike_indices[i-1]+1) and (spike_index != spike_indices[i+1]-1):
                isolated_spike_index.append(i)

        # Determine where the isolated spikes and boxcars are by index
        isolated_spikes = [spike_indices[i] for i in isolated_spike_index]
        boxcar_indices = [spike_index for spike_index in spike_indices if spike_index not in spike_indices[isolated_spike_index]]

        # Execute interpolation between isolated spikes
        for isolated_spike in isolated_spikes:
            tr.data[isolated_spike] = (tr.data[isolated_spike-1] + tr.data[isolated_spike+1])/2

        # Initialize boxcar limits list
        boxcar_limits = []

        # Loop through all boxcar indices
        for i in range(len(boxcar_indices)):
            index = boxcar_indices[i]

            # Record index if it is the first boxcar index
            if index == boxcar_indices[0]:
                boxcar_limits.append(index)

            # Record index if it is the last boxcar index
            elif index == boxcar_indices[-1]:
                boxcar_limits.append(index)

            # Record index if it is a discontinuous member
            elif index != boxcar_indices[i - 1] + 1 or index != boxcar_indices[i + 1] - 1:
                boxcar_limits.append(index)

        # if len(boxcar_limits) % 2 != 0:
        #     print('WARNING: Data has a non-even number of spikes, skipping trace.')
        #     continue

        # Use boxcar limits to obtain data limits for extraction
        boxcar_limits = np.array(boxcar_limits)
        extract_limits = [*sum(zip(boxcar_limits[0::2] - 1, boxcar_limits[1::2] + 1), ())]

        # Add our data's head and tail index to extraction limits
        extract_limits = np.append(extract_limits, len(data) - 1)
        extract_limits = np.insert(extract_limits, 0, 0)

        # Use the sampling rate to determine start time and end time of extraction
        tr_starttime = tr.stats.starttime
        dt = 1 / tr.stats.sampling_rate

        # Loop through limit pairs for data extraction
        for j in range(0, len(extract_limits), 2):

            # Determine start and end time of valid data section
            extract_starttime = tr_starttime + (extract_limits[j] * dt)
            extract_endtime = tr_starttime + (extract_limits[j + 1] * dt)

            # Extract trace and store in output stream
            tr_out = tr.slice(extract_starttime, extract_endtime)
            st_out = st_out + tr_out

    # Return de-spiked and de-boxcared trace
    return st_out

# [remove_bad_traces] function to remove traces with too many zeros
def remove_bad_traces(st,max_zeros=100):

    # Recommended value: max_zeros=100 for 50Hz daylong traces

    # Import dependencies
    import numpy as np

    # Check traces
    for tr in st:
        num_zeros = np.sum(tr.data.data==0)
        if num_zeros > max_zeros:
            st.remove(tr)

    # Return filtered stream object
    return st

# [compare_tr] function to compare the same trace on local data vs its counterpart on IRIS
def compare_tr(tr):

    # Import dependencies
    from obspy import Stream
    from waveform_collection import gather_waveforms

    # Define net-sta-chan, starttime, endtime
    network = tr.stats.network
    station = tr.stats.station
    channel = tr.stats.channel
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime

    # Gather waveforms
    tr_client = gather_waveforms(source='IRIS', network=network, station=station,
                          location='*', channel=channel, starttime=starttime,
                          endtime=endtime)

    # Combine both traces and plot on the same figure
    combine = Stream()
    combine = combine + tr
    combine = combine + tr_client
    combine.plot()

# [remove_repeats] function to clean repeating events in a party
# Note that this has been replaced by party.decluster()
def remove_repeats(party,time_interval):

    # Import dependencies
    import numpy as np
    from eqcorrscan import Party, Family

    # Initialize an inclusion list, list of all i&j, and arrays of detection times, values, and thresholds
    include_list = []
    all_i = []
    all_j = []
    all_detection_time = []
    all_detection_val = []
    all_detection_threshold = []

    # Loop through each family in the party
    for i in range(len(party.families)):

        # Populate include list with ones (we start by including all detections in each family)
        include_list.append(np.ones(len(party[i])))
        family = party[i]

        # Loop through each detection in the family, storing detection info (time, val, threshold)
        for j in range(len(family)):
            all_i.append(i)
            all_j.append(j)
            all_detection_time.append(family[j].detect_time)
            all_detection_val.append(family[j].detect_val)
            all_detection_threshold.append(family[j].threshold)

    # Convert all lists to arrays
    all_detection_time = np.array(all_detection_time)
    all_detection_val = np.array(all_detection_val)
    all_detection_threshold = np.array(all_detection_threshold)

    # Loop through all families' detections again, making comparisons
    for i in range(len(party.families)):

        # Extract family
        family = party[i]

        # Loop through family's detections
        for j in range(len(family)):

            # If the event is already excluded, continue
            if include_list[i][j] == 0:
                continue

            # Calculate time difference with all other detections
            time_difference = abs(all_detection_time-family[j].detect_time)
            matching_index = np.where(time_difference<time_interval)[0]

            # If the only matching event is itself, continue
            if len(matching_index) == 1:
                continue

            # Otherwise, find highest detection value, and exclude other events
            else:
                matching_detection_diff = all_detection_val[matching_index] - all_detection_threshold[matching_index]
                max_detection_index = matching_index[np.argmax(abs(matching_detection_diff))]
                exclude_index = matching_index[matching_index!=max_detection_index]
                for k in exclude_index:
                    include_list[all_i[k]][all_j[k]] = 0

    # Populate new party with only unrepeated detections
    party_clean = Party()
    for i in range(len(party.families)):
        family_clean = Family(template=party[i].template)
        for j in range(len(party[i])):
            if include_list[i][j] == 1:
                family_clean = family_clean + party[i][j]
        party_clean = party_clean + family_clean

    # Return cleaned party object
    return party_clean

# [get_local_stations] function to get all stations within a radius of a volcano
def get_local_stations(list_dir,volcano_name,radius):

    # Import dependencies
    import numpy as np
    import pandas as pd
    import geopy.distance

    # Read in volcano lat and lon
    volcano_list = pd.read_csv(list_dir + 'volcano_list.csv')
    station_list = pd.read_csv(list_dir + 'station_list.csv')
    volcano_names = [name.lower() for name in volcano_list.volcano]
    volcano_index = volcano_names.index(volcano_name.lower())
    volcano_lat = volcano_list.latitude[volcano_index]
    volcano_lon = volcano_list.longitude[volcano_index]
    volcano_coord = (volcano_lat,volcano_lon)

    # Get distances from all stations
    station_distances = []
    for i in range(len(station_list)):
        station_lat = station_list.latitude[i]
        station_lon = station_list.longitude[i]
        station_coord = (station_lat,station_lon)
        station_distance = geopy.distance.GeodesicDistance(volcano_coord,station_coord).km
        station_distances.append(station_distance)

    # Use radius condition to derive list
    local_stations = list(station_list.station[(np.array(station_distances) < radius)])

    # Return list of local stations
    return local_stations

# [gen_detect_hist] function to create detection-threshold histogram
def gen_detect_hist(party):

    # Import dependencies
    import numpy as np
    import matplotlib.pyplot as plt

    # Initialize lists for detection and threshold values
    detection_list = []
    threshold_list = []

    # Loop through each family
    for i in range(len(party.families)):
        family = party[i]

        # Loop through detections in each family
        for j in range(len(family)):

            # Append detection & threshold values for every detection
            detection_value = abs(family[j].detect_val)
            detection_list.append(detection_value)
            threshold_value = family[j].threshold
            threshold_list.append(threshold_value)

    # Plot distribution of detections as histogram and save
    detection_array = np.array(detection_list) - np.array(threshold_list)
    detection_floor = np.floor(min(detection_array))
    detection_ceil = np.ceil(max(detection_array))
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.hist(detection_array, bins=np.arange(detection_floor, detection_ceil, 0.1), color='teal', edgecolor='black')
    ax.set_xlim([detection_floor, max(detection_array)])
    ax.set_xlabel('Detection-Threshold Gap')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of Detection Values')
    fig.show()

# [writer] function to write and save a catalog/tribe/party
def writer(outpath,object):

    # Run quick check for file extension
    if (outpath[-4:] != '.tgz') and (outpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')

    # Import dependencies
    import os
    from obspy import Catalog
    from eqcorrscan import Party
    from eqcorrscan.core.match_filter.tribe import Tribe

    # If the file already exists, remove it
    if os.path.exists(outpath):
        os.remove(outpath)

    # Initialize class samplers
    catalog_type = Catalog()
    party_type = Party()
    tribe_type = Tribe()

    # Use if condition to choose saving style
    if type(object) == type(catalog_type):
        object.write(outpath, format='QUAKEML')
    elif type(object) == type(party_type):
        object.write(outpath, format='tar')
    elif type(object) == type(tribe_type):
        object.write(outpath)
    else:
        raise ValueError('Review writer(): input has invalid object type!')

# [reader] function to read a catalog/tribe/party
def reader(inpath):

    # Run quick check for file extension
    if (inpath[-4:] != '.tgz') and (inpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')

    # Import dependencies
    from obspy import read_events
    from eqcorrscan.core.match_filter import read_party
    from eqcorrscan.core.match_filter.tribe import read_tribe

    # Use if condition to choose reading style
    if ('catalog' or 'events' or 'templates' in inpath) and (inpath[-4:] == '.xml'):
        object = read_events(inpath, format='QUAKEML')
    elif ('party' in inpath) and (inpath[-4:] == '.tgz'):
        object = read_party(inpath, read_detection_catalog=False)
    elif ('tribe' in inpath) and (inpath[-4:] == '.tgz'):
        object = read_tribe(inpath)
    else:
        raise ValueError('Review reader(): input path does not follow naming convention.')

    # Return object that was read
    return object

# [read_trace] read trace from local data files (limited to less than 24h in duration)
def read_trace(data_dir,station,channel,starttime,endtime,tolerance=4e4):

    # Import dependencies
    import glob
    from obspy import Stream, read
    from toolbox import remove_boxcars

    # Initialize list of data filenames for reading
    data_filenames = []

    # Extract year and julday from starttime
    start_year = starttime.year
    start_julday = starttime.julday
    end_year = endtime.year
    end_julday = endtime.julday

    # Craft file string and append
    data_filename = data_dir + station + '.' + channel + '.' + str(start_year) + ':' + f'{start_julday:03}' + ':*'
    data_filenames.append(data_filename)

    # Add another day if starttime and endtime are on different data files
    if start_year != end_year or start_julday != end_julday:
        data_filename = data_dir + station + '.' + channel + '.' + str(end_year) + ':' + f'{end_julday:03}' + ':*'
        data_filenames.append(data_filename)

    # Read in all traces we need
    tr_unmerged = Stream()
    for data_filename in data_filenames:
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                tr_contribution = read(matching_filename)
                tr_unmerged = tr_unmerged + tr_contribution
            except:
                continue

    # Remove boxcars and spikes, detrend, merge, and trim to desired starttime and endtime
    tr_merged = tr_unmerged.copy()
    tr_merged.trim(starttime=starttime, endtime=endtime)
    tr_merged = remove_boxcars(tr_merged, tolerance)
    tr_merged.detrend("simple").merge()

    # Since a single station and channel is provided, the traces should merge into 1
    if len(tr_merged) == 0:
        output = Stream()
    elif len(tr_merged) == 1:
        output = tr_merged[0]
    else:
        raise ValueError('Error in read_trace(), function returns stream rather than trace.')

    # Return trace
    return output

# [prepare_catalog_stream] function to read streams pertaining to an input catalog
def prepare_catalog_stream(data_dir,catalog,resampling_frequency,tolerance):

    # Import dependencies
    import glob
    from obspy import Stream, read
    from toolbox import remove_boxcars

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

    # Remove boxcar and spikes, resample and merge
    stream = remove_boxcars(stream, tolerance)
    stream = stream.resample(resampling_frequency)
    stream = stream.merge()

    # Return stream object
    return stream

# [client_download] function to use IRIS to download the same input streams
def client_download(stream_in,source='IRIS'):

    # Import dependencies
    from obspy import Stream
    from waveform_collection import gather_waveforms

    # Initialize output stream
    stream_out = Stream()

    # Loop through input stream's traces
    for trace_in in stream_in:

        # Gather net-sta-chan and time limits to use gather_waveforms
        network = trace_in.stats.network
        station = trace_in.stats.station
        location = trace_in.stats.location
        channel = trace_in.stats.channel
        starttime = trace_in.stats.starttime
        endtime = trace_in.stats.endtime

        # Use gather_waveforms to pull data from IRIS
        trace_out = gather_waveforms(source=source, network=network, station=station,
                      location=location, channel=channel, starttime=starttime,
                      endtime=endtime)

        # Add trace contribution to output stream
        stream_out = stream_out + trace_out

    # Return stream object
    return(stream_out)

# [get_detection] get stream related to detection by downloading necessary streams
def get_detection(detection,data_dir='/home/data/redoubt/',client_name='IRIS',length=10,resampling_frequency=50,tolerance=4e4,lowcut=1,highcut=10,plot=False):

    # Import dependencies
    import glob
    from obspy import Stream, read
    from obspy.clients.fdsn import Client
    from toolbox import remove_boxcars

    # Extract detection time and sta-chan combinations
    detect_time = detection.detect_time
    year = detect_time.year
    julday = detect_time.julday
    stachans = detection.chans

    # If a data directory is provided,
    if data_dir is not None:

        # Initialize data filename list
        data_filenames = []

        # Loop through sta-chan combinations
        for sta, chan in stachans:

            # Stitch data filename and append
            data_filename = data_dir + sta + '.' + chan + '.' + str(year) + ':' + f'{julday:03}' + ':*'
            data_filenames.append(data_filename)

        # Read in all required streams
        stream = Stream()
        for data_filename in data_filenames:
            matching_filenames = (glob.glob(data_filename))
            for matching_filename in matching_filenames:
                try:
                    stream_contribution = read(matching_filename)
                    stream = stream + stream_contribution
                except:
                    continue

    # Otherwise, download data from IRIS
    else:

        # Define client and target network
        client = Client(client_name)
        net = detection.event.picks[0].waveform_id.network_code

        # Get waveforms
        stream = Stream()

        # Download all waveforms needed
        for sta, chan in stachans:
            stream_contribution = client.get_waveforms(net, sta, "*", chan, detect_time, detect_time+length)
            stream = stream + stream_contribution


    # Enforce desired length
    stream.trim(starttime=detect_time,endtime=detect_time+length)

    # Resample all traces to prevent inconsistencies
    stream.resample(resampling_frequency)

    # Remove boxcar and spikes, detrend and merge
    if data_dir is None:
        stream = remove_boxcars(stream, tolerance)
        stream = stream.detrend("simple")

    # Filter and taper
    stream = stream.filter('bandpass',freqmax=highcut,freqmin=lowcut)
    stream = stream.taper(0.05, type='hann', max_length=(0.75 * 1024 / 50))
    stream = stream.merge()

    # Plot if desired
    if plot:
        stream.plot(color='b',equal_scale=False, size=(800, 600))

    # Return stream object
    return stream

# [prepare_stream_dict] prepare stream dictionary pertaining to input catalog
def prepare_stream_dict(catalog,pre_pick,length,local=False,client_name="IRIS",data_dir=None,resampling_frequency=50):

    # Import dependencies
    from toolbox import read_trace
    from obspy import Stream
    from obspy.clients.fdsn import Client

    # Define client if not using local
    if not local:
        client = Client(client_name)

    # Initialize list to contain tuples
    stream_tuples = []

    # Loop through events in catalog
    for event in catalog:

        # Extract resource id and initialize stream object
        event_id = event.resource_id.id
        event_st = Stream()

        # Loop through picks in event
        for pick in event.picks:

            # Extract waveform details from pick object
            network = pick.waveform_id.network_code
            station = pick.waveform_id.station_code
            channel = pick.waveform_id.channel_code
            starttime = pick.time - pre_pick
            endtime = starttime + length

            # Read trace from local data directory if local, use IRIS if not local
            if local:
                tr = read_trace(data_dir, station, channel, starttime, endtime, tolerance=4e4)
            else:
                tr = client.get_waveforms(network, station, '*', channel, starttime, endtime)
                tr = tr.resample(resampling_frequency)  # for my EQcorrscan attempt

            # Add trace to stream object
            event_st += tr

        # Split stream and append by alongside event id
        event_st = event_st.split()
        stream_tuples.append((event_id,event_st))

    # Convert tuple object to dict object
    stream_dict = dict(stream_tuples)

    # Return stream dictionary
    return stream_dict

# raster2array from https://www.neonscience.org/resources/learning-hub/tutorials/merge-lidar-geotiff-py
def raster2array(geotif_file):
    metadata = {}
    dataset = gdal.Open(geotif_file)
    metadata['array_rows'] = dataset.RasterYSize
    metadata['array_cols'] = dataset.RasterXSize
    metadata['bands'] = dataset.RasterCount
    metadata['driver'] = dataset.GetDriver().LongName
    metadata['projection'] = dataset.GetProjection()
    metadata['geotransform'] = dataset.GetGeoTransform()

    mapinfo = dataset.GetGeoTransform()
    metadata['pixelWidth'] = mapinfo[1]
    metadata['pixelHeight'] = mapinfo[5]

    xMin = mapinfo[0]
    xMax = mapinfo[0] + dataset.RasterXSize/mapinfo[1]
    yMin = mapinfo[3] + dataset.RasterYSize/mapinfo[5]
    yMax = mapinfo[3]

    metadata['extent'] = (xMin,xMax,yMin,yMax)

    raster = dataset.GetRasterBand(1)
    array_shape = raster.ReadAsArray(0,0,metadata['array_cols'],metadata['array_rows']).astype(np.float).shape
    metadata['noDataValue'] = raster.GetNoDataValue()
    metadata['scaleFactor'] = raster.GetScale() or 1

    array = np.zeros((array_shape[0],array_shape[1],dataset.RasterCount),'uint8') #pre-allocate stackedArray matrix

    if metadata['bands'] == 1:
        raster = dataset.GetRasterBand(1)
        metadata['noDataValue'] = raster.GetNoDataValue()
        metadata['scaleFactor'] = raster.GetScale() or 1

        array = dataset.GetRasterBand(1).ReadAsArray(0,0,metadata['array_cols'],metadata['array_rows']).astype(np.float)
        #array[np.where(array==metadata['noDataValue'])]=np.nan
        array = array/metadata['scaleFactor']

    elif metadata['bands'] > 1:
        for i in range(1, dataset.RasterCount+1):
            band = dataset.GetRasterBand(i).ReadAsArray(0,0,metadata['array_cols'],metadata['array_rows']).astype(np.float)
            #band[np.where(band==metadata['noDataValue'])]=np.nan
            band = band/metadata['scaleFactor']
            array[...,i-1] = band

    return array, metadata