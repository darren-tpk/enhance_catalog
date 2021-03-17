# function to remove artificial boxcar-like signals from stream
def remove_boxcars(st,tolerance):
    # import packages
    import numpy as np
    from obspy import Stream
    # initialize output stream
    st_out = Stream()
    # loop through traces in stream
    for tr in st:
        # detrend data and find spikes above tolerance (median(abs(data)) multiplier)
        data = tr.data
        tr_detrended = tr.copy().detrend()
        data_detrended = tr_detrended.data
        spike_index = np.where(data_detrended > tolerance * np.median(abs(data_detrended)))[0]
        # extract spike limits
        spike_limits = []
        for i in range(len(spike_index)):
            index = spike_index[i]
            if index == spike_index[0]:
                spike_limits.append(index)
            elif index == spike_index[-1]:
                spike_limits.append(index)
            elif index != spike_index[i - 1] + 1 or index != spike_index[i + 1] - 1:
                spike_limits.append(index)
        if len(spike_limits) % 2 != 0:
            raise ValueError('Review remove_boxcars: The number of spike limits is not even!')
        # knowing our spike limits, we get data limits for extraction
        spike_limits = np.array(spike_limits)
        extract_limits = [*sum(zip(spike_limits[0::2] - 1, spike_limits[1::2] + 1), ())]
        # add data's head and tail index to extraction limits
        extract_limits = np.append(extract_limits, len(data) - 1)
        extract_limits = np.insert(extract_limits, 0, 0)
        # using the sampling rate, we determine start time and end time of extraction
        tr_starttime = tr.stats.starttime
        dt = 1 / tr.stats.sampling_rate
        # loop through limit pairs
        for j in range(0, len(extract_limits), 2):
            extract_starttime = tr_starttime + (extract_limits[j] * dt)
            extract_endtime = tr_starttime + (extract_limits[j + 1] * dt)
            # extract trace and store in output stream
            tr_out = tr.slice(extract_starttime, extract_endtime)
            st_out = st_out + tr_out
    return st_out

# function to compare trace with data on IRIS
def compare_tr(tr):
    # import packages
    from obspy import Stream
    from waveform_collection import gather_waveforms
    # define net-sta-chan, starttime, endtime
    network = tr.stats.network
    station = tr.stats.station
    channel = tr.stats.channel
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    # gather waveforms
    tr_client = gather_waveforms(source='IRIS', network=network, station=station,
                          location='*', channel=channel, starttime=starttime,
                          endtime=endtime)
    # combine and plot on same figure
    combine = Stream()
    combine = combine + tr
    combine = combine + tr_client
    combine.plot()

# function to clean repeating events in a party
def remove_repeats(party,time_interval):
    # import packages
    import numpy as np
    from eqcorrscan import Party, Family
    # initialize an inclusion list, list of all i&j, and arrays of detection times, values, and thresholds
    include_list = []
    all_i = []
    all_j = []
    all_detection_time = []
    all_detection_val = []
    all_detection_threshold = []
    # loop through each family in the party
    for i in range(len(party.families)):
        # populate include list with ones (we start by including all detections in each family)
        include_list.append(np.ones(len(party[i])))
        family = party[i]
        # loop through each detection in the family, storing key info
        for j in range(len(family)):
            all_i.append(i)
            all_j.append(j)
            all_detection_time.append(family[j].detect_time)
            all_detection_val.append(family[j].detect_val)
            all_detection_threshold.append(family[j].threshold)
    # convert to arrays
    all_detection_time = np.array(all_detection_time)
    all_detection_val = np.array(all_detection_val)
    all_detection_threshold = np.array(all_detection_threshold)
    # now loop through all families' detections again, making comparisons
    for i in range(len(party.families)):
        # extract family
        family = party[i]
        # loop through family's detections
        for j in range(len(family)):
            # if the event is already excluded, continue
            if include_list[i][j] == 0:
                continue
            # calculate time difference for all other detections
            time_difference = abs(all_detection_time-family[j].detect_time)
            matching_index = np.where(time_difference<time_interval)[0]
            # if the only matching event is itself, continue
            if len(matching_index) == 1:
                continue
            # otherwise, find highest detection value, and exclude other events
            else:
                matching_detection_diff = all_detection_val[matching_index] - all_detection_threshold[matching_index]
                max_detection_index = matching_index[np.argmax(abs(matching_detection_diff))]
                exclude_index = matching_index[matching_index!=max_detection_index]
                for k in exclude_index:
                    include_list[all_i[k]][all_j[k]] = 0
    # lastly, populate new party with only unrepeated detections
    party_clean = Party()
    for i in range(len(party.families)):
        family_clean = Family(template=party[i].template)
        for j in range(len(party[i])):
            if include_list[i][j] == 1:
                family_clean = family_clean + party[i][j]
        party_clean = party_clean + family_clean
    return party_clean

# function to get all stations within a radius of a volcano
def get_local_stations(volcano_name,radius):
    # import packages
    import numpy as np
    import pandas as pd
    import geopy.distance
    # get volcano lat and lon
    volcano_list = pd.read_csv('/home/ptan/project/avo_data/volcano_list.csv')
    station_list = pd.read_csv('/home/ptan/project/avo_data/station_list.csv')
    volcano_names = [name.lower() for name in volcano_list.volcano]
    volcano_index = volcano_names.index(volcano_name.lower())
    volcano_lat = volcano_list.latitude[volcano_index]
    volcano_lon = volcano_list.longitude[volcano_index]
    volcano_coord = (volcano_lat,volcano_lon)
    # now get distances with all stations
    station_distances = [];
    for i in range(len(station_list)):
        station_lat = station_list.latitude[i]
        station_lon = station_list.longitude[i]
        station_coord = (station_lat,station_lon)
        station_distance = geopy.distance.GeodesicDistance(volcano_coord,station_coord).km
        station_distances.append(station_distance)
    # use radius condition to derive list
    local_stations = list(station_list.station[(np.array(station_distances) < radius)])
    return local_stations

# function to create detection-threshold histogram
def gen_detect_hist(party):
    # import pacakges
    import numpy as np
    import matplotlib.pyplot as plt
    # Initialize lists for detection and threshold values
    detection_list = []
    threshold_list = []
    # loop through each family
    for i in range(len(party.families)):
        family = party[i]
        # loop through detections in each family
        for j in range(len(family)):
            # append detection & threshold values for every detection
            detection_value = abs(family[j].detect_val)
            detection_list.append(detection_value)
            threshold_value = family[j].threshold
            threshold_list.append(threshold_value)
    # plot distribution of detections as histogram and save
    detection_array = np.array(detection_list) - np.array(threshold_list)
    detection_floor = np.floor(min(detection_array))
    detection_ceil = np.ceil(max(detection_array))
    threshold_array = np.unique(threshold_list)
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.hist(detection_array, bins=np.arange(detection_floor, detection_ceil, 0.1), color='teal', edgecolor='black')
    ax.set_xlim([detection_floor, max(detection_array)])
    ax.set_xlabel('Detection-Threshold Gap')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of Detection Values')
    fig.show()

def writer(outpath,object):
    # check for file extension
    if (outpath[-4:] != '.tgz') and (outpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')
    # import pacakges and object classes
    import os
    from obspy import Catalog
    from eqcorrscan import Party
    from eqcorrscan.core.match_filter.tribe import Tribe
    # if the file already exists, remove it
    if os.path.exists(outpath):
        os.remove(outpath)
    # initialize class samplers
    catalog_type = Catalog()
    party_type = Party()
    tribe_type = Tribe()
    # use if condition to choose saving style
    if type(object) == type(catalog_type):
        object.write(outpath, format='QUAKEML')
    elif type(object) == type(party_type):
        object.write(outpath, format='tar')
    elif type(object) == type(tribe_type):
        object.write(outpath)
    else:
        raise ValueError('Review writer(): input has invalid object type!')

def reader(inpath):
    # check for file extension
    if (inpath[-4:] != '.tgz') and (inpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')
    # import packages
    from obspy import read_events
    from eqcorrscan import Party
    from eqcorrscan.core.match_filter.tribe import read_tribe
    # use if condition to choose reading style
    if ('catalog' in inpath) and (inpath[-4:] == '.xml'):
        object = read_events(inpath, format='QUAKEML')
    elif ('party' in inpath) and (inpath[-4:] == '.tgz'):
        object = Party().read(inpath)
    elif ('tribe' in inpath) and (inpath[-4:] == '.tgz'):
        object = read_tribe(inpath)
    else:
        raise ValueError('Review reader(): input path does not follow naming convention.')
    return object

