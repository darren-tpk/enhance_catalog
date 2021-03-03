
def remove_boxcars(st,tolerance):
    import numpy as np
    from obspy import Stream
    st_out = Stream()
    for tr in st:
        data = tr.data
        tr_detrended = tr.copy().detrend()
        data_detrended = tr_detrended.data
        spike_index = np.where(data_detrended > tolerance * np.median(abs(data_detrended)))[0]
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
        spike_limits = np.array(spike_limits)
        extract_limits = [*sum(zip(spike_limits[0::2] - 1, spike_limits[1::2] + 1), ())]
        extract_limits = np.append(extract_limits, len(data) - 1)
        extract_limits = np.insert(extract_limits, 0, 0)
        tr_starttime = tr.stats.starttime
        dt = 1 / tr.stats.sampling_rate
        for j in range(0, len(extract_limits), 2):
            extract_starttime = tr_starttime + (extract_limits[j] * dt)
            extract_endtime = tr_starttime + (extract_limits[j + 1] * dt)
            tr_out = tr.slice(extract_starttime, extract_endtime)
            st_out = st_out + tr_out
    return st_out

def compare_tr(tr):
    from obspy import Stream
    from waveform_collection import gather_waveforms
    network = tr.stats.network
    station = tr.stats.station
    channel = tr.stats.channel
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    tr_client = gather_waveforms(source='IRIS', network=network, station=station,
                          location='*', channel=channel, starttime=starttime,
                          endtime=endtime)
    combine = Stream()
    combine = combine + tr
    combine = combine + tr_client
    combine.plot()

### clean repeating events
def remove_repeats(party,time_interval):
    import numpy as np
    from eqcorrscan import Party, Family
    # craft an inclusion list, list of indices i&j
    # and arrays of detection times & values
    include_list = []
    all_i = []
    all_j = []
    all_detection_time = []
    all_detection_val = []
    all_detection_threshold = []
    for i in range(len(party.families)):
        include_list.append(np.ones(len(party[i])))
        family = party[i]
        for j in range(len(family)):
            all_i.append(i)
            all_j.append(j)
            all_detection_time.append(family[j].detect_time)
            all_detection_val.append(family[j].detect_val)
            all_detection_threshold.append(family[j].threshold)
    all_detection_time = np.array(all_detection_time)
    all_detection_val = np.array(all_detection_val)
    all_detection_threshold = np.array(all_detection_threshold)
    # now loop through all detections again, making comparisons
    for i in range(len(party.families)):
        family = party[i]
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
            # otherwise, find highest detection value, and set exclude other events
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

def get_local_stations(volcano_name,radius):
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
    # now get all distances with all stations
    station_distances = [];
    for i in range(len(station_list)):
        station_lat = station_list.latitude[i]
        station_lon = station_list.longitude[i]
        station_coord = (station_lat,station_lon)
        station_distance = geopy.distance.GeodesicDistance(volcano_coord,station_coord).km
        station_distances.append(station_distance)
    local_stations = list(station_list.station[(np.array(station_distances) < radius)])
    return local_stations

# for i in range(len(tribe)):
#     template = tribe[i]
#     for trace in template.st:
#         station = trace.stats.station
#         if station not in local_stations:
#             print(i, station)



