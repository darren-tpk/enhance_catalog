#%% Toolbox of all supporting functions

#%% [pull_cores] Pulls out catalog of cores from a mixed catalog
def pull_cores(full_catalog):

    # Import all dependencices
    import time
    from obspy import Catalog

    print('\nExtracting core events...')
    time_start = time.time()

    # Initialize core catalog object
    all_cores = Catalog()

    # Loop through full catalog to check for cores
    for detection in full_catalog:
        if detection.origins[0].comments[0].text.split(' ')[2] == 'core':
            all_cores.append(detection)

    # Find repeats and save their indices
    remove_index = []
    for ii in range(1,len(all_cores)):
        if all_cores[ii].origins[0].time == all_cores[ii - 1].origins[0].time:
            remove_index.append(ii)

    # Reappend detections to core catalog, avoiding repeat indices
    core_catalog = Catalog()
    for ii, core in enumerate(all_cores):
        if ii not in remove_index:
            core_catalog.append(core)

    # Conclude process
    time_stop = time.time()
    print('Core extraction complete, processing time: %.2f s' % (time_stop - time_start))
    return core_catalog

#%% [read_hypoi] Reads in a hypoinverse.pha file as an obspy catalog with numerous catalog filtering options.
def read_hypoi(hypoi_file,
               time_lim=None,
               radial_lim=None,
               bbox_lim=None,
               nsta_min=None,
               nobs_min=None,
               rmse_max=None,
               gap_max=None,
               dist2sta_max=None,
               eh_max=None,
               ev_max=None,
               depth_req=True,
               depth_lim=None,
               mag_req=True,
               mag_lim=None,
               num_P_min=None,
               num_S_min=None,
               channel_convention=True,
               summary=True):
    """
    Reads in a hypoinverse.pha file as an obspy catalog with numerous catalog filtering options.
    :param hypoi_file (str): file path for hypoinverse.pha file
    :param time_lim (tuple): time bounds for events in output catalog (length 2 tuple with UTCDateTime entries)
    :param radial_lim (tuple): radial distance bound for events in output catalog (length 3 tuple storing latitude, longitude, in decimal degrees and radius in km)
    :param bbox_lim (tuple): bounding box for events in output catalog (length 4 tuple storing latitude and longitude of the top left corner followed by the bottom right corner)
    :param nsta_min (int): minimum number of picked stations for an event to be kept
    :param nobs_min (int): minimum number of valid picks (weight >= 0.1) for an event to be kept
    :param rmse_max (float): maximum travel time root-mean-squared-error for an event to be kept
    :param gap_max (degrees): maximum azimuthal gap between picked stations for an event to be kept
    :param dist2sta_max (float): maximum event-to-nearest-station distance allowed for an event to be kept
    :param eh_max (float): maximum horizontal error allowed to be kept (km)
    :param ev_max (float): maximum vertical error allowed for an event to be kept (km)
    :param depth_req (bool): if `True`, only events with a valid depth entry are kept
    :param depth_lim (tuple): depth bounds for an event to be kept (length 2 tuple storing minimum and maximum depth in m)
    :param mag_req (bool): if `True`, only events with a valid magnitude entry are kept
    :param mag_lim (tuple): magnitude bounds for an event to be kept (length 2 tuple storing minimum and maximum magnitude)
    :param num_P_min (int): minimum number of valid (weight >= 0.1) P phases needed for an event to be kept
    :param num_S_min (int): minimum number of valid (weight >= 0.1) S phases needed for an event to be kept
    :param channel_convention (bool): if `True`, enforces channel convention for P and S picks (P on vertical, S on horizontal)
    :param summary (bool): if `True`, prints out a summary of all conditions imposed
    :return: obspy Catalog object
    """

    # Import dependencies
    import pandas as pd
    import numpy as np
    from obspy import UTCDateTime, Catalog
    from obspy.core.event import Event, Origin, Magnitude, Comment, Pick, WaveformStreamID, Arrival, ResourceIdentifier, \
        OriginQuality
    import geopy.distance
    import time

    # Start timer for code
    time_start = time.time()

    # Read input file
    data_frame = pd.read_csv(hypoi_file, header=None, skipfooter=1, engine='python')
    input = data_frame[0]

    # Initialize catalog object and header count
    catalog = Catalog()
    total_headers = 0

    # Initialize valid event
    valid_event = False

    # Loop through every line in input file
    for i in range(len(input)):

        # Extract line
        line = input[i]

        # Check what type of line it is

        # If it is a new earthquake location line
        if '19' in line[0:2] or '20' in line[0:2]:

            # Add header count
            total_headers += 1

            # Initialize valid event bool
            valid_event = True

            # Construct latitude (convert to decimal degrees)
            if 'S' in line[18]:
                latitude = -1 * (float(line[16:18]) + (0.01 * float(line[19:23]) / 60))
            elif ' ' in line[18]:
                latitude = (float(line[16:18]) + (0.01 * float(line[19:23]) / 60))
            else:
                raise ValueError("Hypoinverse file has an invalid latitude entry.")

            # Construct longitude (convert to decimal degrees)
            if 'E' in line[26]:
                longitude = (float(line[23:26]) + (0.01 * float(line[27:31]) / 60))
            elif ' ' in line[26]:
                longitude = -1 * (float(line[23:26]) + (0.01 * float(line[27:31]) / 60))
            else:
                raise ValueError("Hypoinverse file has an invalid longitude entry.")

            # Check if event falls within radial limits
            if radial_lim is not None:
                if len(radial_lim) != 3:
                    raise ValueError(
                        'Input radial_lim provided is not of length 3. Please follow the (latitude,longitude,radius) format.')
                circle_center = radial_lim[0:2]
                radial_distance = geopy.distance.GeodesicDistance((latitude, longitude), circle_center).km
                if radial_distance > radial_lim[-1]:
                    valid_event = False
                    continue

            # Check if event falls within bounding box limits
            if bbox_lim is not None:
                if type(bbox_lim) != tuple and len(bbox_lims) != 4:
                    raise ValueError(
                        'Input bbox_lim provided is not of length 4. Please follow the (latitude1,longitude1,latitude2,longitude2) format.')
                if latitude < bbox_lim[0] or latitude > bbox_lim[2] or longitude > bbox_lim[1] or longitude < bbox_lim[
                    3]:
                    valid_event = False
                    continue

            # Construct event UTCDateTime
            year = int(line[0:4])
            month = int(line[4:6])
            day = int(line[6:8])
            hour = int(line[8:10])
            minute = int(line[10:12])
            second = 0.01 * float(line[12:16])
            datetime = UTCDateTime(year, month, day, hour, minute, second)

            # Check if event falls within temporal limits
            if time_lim is not None:
                if datetime < time_lim[0] or datetime > time_lim[1]:
                    valid_event = False
                    continue

            # Extract most filter params (except for mag and depth requirement)
            nobs = int(line[39:42])
            rmse = 0.01 * float(line[48:52])  # travel time residual
            gap = int(line[42:45])
            dist2sta = int(line[45:48])
            eh = 0.01 * float(line[85:89])
            ev = 0.01 * float(line[89:93])
            num_S = int(line[82:85])
            num_P = nobs - num_S

            # Now do filtering
            if nobs_min is not None and nobs < nobs_min:
                valid_event = False
                continue
            if rmse_max is not None and rmse > rmse_max:
                valid_event = False
                continue
            if gap_max is not None and gap > gap_max:
                valid_event = False
                continue
            if dist2sta_max is not None and dist2sta > dist2sta_max:
                valid_event = False
                continue
            if eh_max is not None and eh > eh_max:
                valid_event = False
                continue
            if ev_max is not None and ev > ev_max:
                valid_event = False
                continue
            if num_P_min is not None and num_P < num_P_min:
                valid_event = False
                continue
            if num_S_min is not None and num_S < num_S_min:
                valid_event = False
                continue

            # Now do check for depth
            if (line[31:36]).isspace():
                if depth_req:
                    valid_event = False
                    continue
                else:
                    depth = -999
                    print('WARNING: The catalog has an event that does not have a depth entry. Setting depth=-999.')
            else:
                depth = float(line[31:36]) / 100
                if depth_lim is not None and (depth < depth_lim[0] or depth > depth_lim[1]):
                    valid_event = False
                    continue

            # Now do check for magnitude
            if (line[147:150]).isspace():
                if mag_req:
                    valid_event = False
                    continue
                else:
                    mag = -999
                    print('WARNING: The catalog has an event that does not have a depth entry. Setting magnitude=-999.')
            else:
                mag = float(line[147:150]) / 100
                if mag_lim is not None and (mag < mag_lim[0] or mag > mag_lim[1]):
                    valid_event = False
                    continue

            # With all checks done, we craft the event
            event = Event(resource_id=ResourceIdentifier(id="smi:local/event/" + line[136:146].strip()),
                          preferred_origin_id=ResourceIdentifier(id="smi:local/origin/" + line[136:146].strip()),
                          preferred_magnitude_id=ResourceIdentifier(id="smi:local/magnitude/" + line[136:146].strip()),
                          origins=[
                              Origin(resource_id=ResourceIdentifier(id="smi:local/origin/" + line[136:146].strip()),
                                     time=datetime, latitude=latitude, longitude=longitude, depth=depth,
                                     quality=OriginQuality(standard_error=rmse))],
                          magnitudes=[Magnitude(
                              resource_id=ResourceIdentifier(id="smi:local/magnitude/" + line[136:146].strip()),
                              mag=mag)])

        # If it is a new phase line with potential S picks:
        elif line[0:2].isalpha():

            # Check if the event header is valid in the first place. If event was kicked out, skip all phase lines
            if not valid_event:
                continue

            # Check for phase hint
            if channel_convention == True:

                # Check if vertical component has P arrival
                if (line[11] == 'Z') and (line[14] == 'P'):

                    # Convert P weight
                    p_weight_code = int(line[16])
                    if p_weight_code == 0:
                        p_weight = 1.0
                    elif p_weight_code == 1:
                        p_weight = 0.5
                    elif p_weight_code == 2:
                        p_weight = 0.2
                    elif p_weight_code == 3:
                        p_weight = 0.1
                    else:
                        # p weight is zero and we skip the phase
                        continue

                    # Construct pick UTCDateTime
                    phase_year = int(line[17:21])
                    phase_month = int(line[21:23])
                    phase_day = int(line[23:25])
                    phase_hour = int(line[25:27])
                    phase_minute = int(line[27:29])
                    phase_second = 0.01 * float(line[29:34])
                    # Fix if phase seconds spill over to the next minute (>= 60s)
                    if phase_second >= 60:
                        if phase_minute < 59:
                            phase_minute = phase_minute + 1
                            phase_second = phase_second - 60
                        else:
                            phase_minute = 0
                            phase_second = phase_second - 60
                    phase_datetime = UTCDateTime(phase_year, phase_month, phase_day, phase_hour, phase_minute,
                                                 phase_second)

                    # Construct waveform id
                    network = line[5:7]
                    station = (line[0:5]).strip()  # input is left justified
                    channel = line[9:12]
                    location = line[111:113]
                    if location == '--':
                        location = ''

                    # Construct and add pick
                    pick = Pick(
                        waveform_id=WaveformStreamID(network_code=network, station_code=station, channel_code=channel,
                                                     location_code=location),
                        time=phase_datetime, phase_hint='P')
                    event.picks.append(pick)

                    # Read in other arrival information
                    p_travel_time_residual = 0.01 * float(line[34:38])
                    # p_weight_used = 0.01 * float(line[38:41])
                    # p_delay_time = 0.01 * float(line[66:70])
                    phase_epicentral_distance = 0.1 * float(line[74:78])
                    phase_emergence_angle = int(line[78:81])
                    phase_azimuth_to_station = int(line[91:94])
                    p_importance = 0.001 * float(line[100:104])
                    # If phase importance is larger than 0.5, hypoinverse weight is flagged as negative for auto-inclusion
                    # See HypoDD manual for Hat matrix explanation
                    if p_importance > 0.5:
                        p_weight = -abs(p_weight)

                    # Construct and add arrival
                    arrival = Arrival(phase='P',
                                      azimuth=phase_azimuth_to_station,
                                      distance=phase_epicentral_distance,
                                      takeoff_angle=phase_emergence_angle,
                                      time_residual=p_travel_time_residual,
                                      time_weight=p_weight)
                    event.origins[0].arrivals.append(arrival)

                # Check if horizontal component has S arrival
                elif (line[11] in ['E', 'N', '1', '2']) and len(line) > 46 and (line[47] == 'S'):

                    # Convert S weight
                    s_weight_code = int(line[49])
                    if s_weight_code == 0:
                        s_weight = 1.0
                    elif s_weight_code == 1:
                        s_weight = 0.5
                    elif s_weight_code == 2:
                        s_weight = 0.2
                    elif s_weight_code == 3:
                        s_weight = 0.1
                    else:
                        # s weight is zero and we skip the phase
                        continue

                    # Construct pick UTCDateTime
                    phase_year = int(line[17:21])
                    phase_month = int(line[21:23])
                    phase_day = int(line[23:25])
                    phase_hour = int(line[25:27])
                    phase_minute = int(line[27:29])
                    phase_second = 0.01 * float(line[41:46])
                    # Fix if phase seconds spill over to the next minute (>= 60s)
                    if phase_second >= 60:
                        if phase_minute < 59:
                            phase_minute = phase_minute + 1
                            phase_second = phase_second - 60
                        else:
                            phase_minute = 0
                            phase_second = phase_second - 60
                    phase_datetime = UTCDateTime(phase_year, phase_month, phase_day, phase_hour, phase_minute,
                                                 phase_second)

                    # Construct waveform id
                    network = line[5:7]
                    station = (line[0:5]).strip()  # input is left justified
                    channel = line[9:12]
                    location = line[111:113]
                    if location == '--':
                        location = ''

                    # Construct and add pick
                    pick = Pick(
                        waveform_id=WaveformStreamID(network_code=network, station_code=station, channel_code=channel,
                                                     location_code=location),
                        time=phase_datetime, phase_hint='S')
                    event.picks.append(pick)

                    # Read in other arrival information
                    s_travel_time_residual = 0.01 * float(line[50:54])
                    # s_weight_used = 0.01 * float(line[63:66])
                    # s_delay_time = 0.01 * float(line[70:74])
                    phase_epicentral_distance = 0.1 * float(line[74:78])
                    phase_emergence_angle = int(line[78:81])
                    phase_azimuth_to_station = int(line[91:94])
                    s_importance = 0.001 * float(line[104:108])
                    # If phase importance is larger than 0.5, hypoinverse weight is flagged as negative for auto-inclusion
                    # See HypoDD manual for Hat matrix explanation
                    if s_importance > 0.5:
                        s_weight = -abs(s_weight)

                    # Construct and add arrival
                    arrival = Arrival(phase='S',
                                      azimuth=phase_azimuth_to_station,
                                      distance=phase_epicentral_distance,
                                      takeoff_angle=phase_emergence_angle,
                                      time_residual=s_travel_time_residual,
                                      time_weight=s_weight)
                    event.origins[0].arrivals.append(arrival)

            # If we are not using channel conventions, then we can have P and S picks on any component
            else:

                # Check if this line has a P pick
                if line[14] == 'P':

                    # Convert P weight
                    p_weight_code = int(line[16])
                    if p_weight_code == 0:
                        p_weight = 1.0
                    elif p_weight_code == 1:
                        p_weight = 0.5
                    elif p_weight_code == 2:
                        p_weight = 0.2
                    elif p_weight_code == 3:
                        p_weight = 0.1
                    else:
                        # p weight is zero and we skip the phase
                        continue

                    # Construct pick UTCDateTime
                    phase_year = int(line[17:21])
                    phase_month = int(line[21:23])
                    phase_day = int(line[23:25])
                    phase_hour = int(line[25:27])
                    phase_minute = int(line[27:29])
                    phase_second = 0.01 * float(line[29:34])
                    # Fix if phase seconds spill over to the next minute (>= 60s)
                    if phase_second >= 60:
                        if phase_minute < 59:
                            phase_minute = phase_minute + 1
                            phase_second = phase_second - 60
                        else:
                            phase_minute = 0
                            phase_second = phase_second - 60
                    phase_datetime = UTCDateTime(phase_year, phase_month, phase_day, phase_hour, phase_minute,
                                                 phase_second)

                    # Construct waveform id
                    network = line[5:7]
                    station = (line[0:5]).strip()  # input is left justified
                    channel = line[9:12]
                    location = line[111:113]
                    if location == '--':
                        location = ''

                    # Construct and add pick
                    pick = Pick(
                        waveform_id=WaveformStreamID(network_code=network, station_code=station, channel_code=channel,
                                                     location_code=location),
                        time=phase_datetime, phase_hint='P')
                    event.picks.append(pick)

                    # Read in other arrival information
                    p_travel_time_residual = 0.01 * float(line[34:38])
                    # p_weight_used = 0.01 * float(line[38:41])
                    # p_delay_time = 0.01 * float(line[66:70])
                    phase_epicentral_distance = 0.1 * float(line[74:78])
                    phase_emergence_angle = int(line[78:81])
                    phase_azimuth_to_station = int(line[91:94])
                    p_importance = 0.001 * float(line[100:104])
                    # If phase importance is larger than 0.5, hypoinverse weight is flagged as negative for auto-inclusion
                    # See HypoDD manual for Hat matrix explanation
                    if p_importance > 0.5:
                        p_weight = -abs(p_weight)

                    # Construct and add arrival
                    arrival = Arrival(phase='P',
                                      azimuth=phase_azimuth_to_station,
                                      distance=phase_epicentral_distance,
                                      takeoff_angle=phase_emergence_angle,
                                      time_residual=p_travel_time_residual,
                                      time_weight=p_weight)
                    event.origins[0].arrivals.append(arrival)

                # Check if this line has an S pick:
                if len(line) > 46 and line[47] == 'S':

                    # Convert S weight
                    s_weight_code = int(line[49])
                    if s_weight_code == 0:
                        s_weight = 1.0
                    elif s_weight_code == 1:
                        s_weight = 0.5
                    elif s_weight_code == 2:
                        s_weight = 0.2
                    elif s_weight_code == 3:
                        s_weight = 0.1
                    else:
                        # s weight is zero and we skip the phase
                        continue

                    # Construct pick UTCDateTime
                    phase_year = int(line[17:21])
                    phase_month = int(line[21:23])
                    phase_day = int(line[23:25])
                    phase_hour = int(line[25:27])
                    phase_minute = int(line[27:29])
                    phase_second = 0.01 * float(line[41:46])
                    # Fix if phase seconds spill over to the next minute (>= 60s)
                    if phase_second >= 60:
                        if phase_minute < 59:
                            phase_minute = phase_minute + 1
                            phase_second = phase_second - 60
                        else:
                            phase_minute = 0
                            phase_second = phase_second - 60
                    phase_datetime = UTCDateTime(phase_year, phase_month, phase_day, phase_hour, phase_minute,
                                                 phase_second)

                    # Construct waveform id
                    network = line[5:7]
                    station = (line[0:5]).strip()  # input is left justified
                    channel = line[9:12]
                    location = line[111:113]
                    if location == '--':
                        location = ''

                    # Construct and add pick
                    pick = Pick(
                        waveform_id=WaveformStreamID(network_code=network, station_code=station, channel_code=channel,
                                                     location_code=location),
                        time=phase_datetime, phase_hint='S')
                    event.picks.append(pick)

                    # Read in other arrival information
                    s_travel_time_residual = 0.01 * float(line[50:54])
                    # s_weight_used = 0.01 * float(line[63:66])
                    # s_delay_time = 0.01 * float(line[70:74])
                    phase_epicentral_distance = 0.1 * float(line[74:78])
                    phase_emergence_angle = int(line[78:81])
                    phase_azimuth_to_station = int(line[91:94])
                    s_importance = 0.001 * float(line[104:108])
                    # If phase importance is larger than 0.5, hypoinverse weight is flagged as negative for auto-inclusion
                    # See HypoDD manual for Hat matrix explanation
                    if s_importance > 0.5:
                        s_weight = -abs(s_weight)

                    # Construct and add arrival
                    arrival = Arrival(phase='S',
                                      azimuth=phase_azimuth_to_station,
                                      distance=phase_epicentral_distance,
                                      takeoff_angle=phase_emergence_angle,
                                      time_residual=s_travel_time_residual,
                                      time_weight=s_weight)
                    event.origins[0].arrivals.append(arrival)

        # If it is shadow line or blank line, skip:
        elif '$' in line[0] or '  ' in line[0:2]:

            # If it has passed all the checks, valid_event = True
            if valid_event:

                # Execute final check of nsta
                nsta = len(np.unique([pick.waveform_id.station_code for pick in event.picks]))
                if nsta_min is not None and nsta < nsta_min:
                    continue
                # If it passes, append the event
                else:
                    catalog.events.append(event)

            # Reset valid event boolean
            valid_event = False

    # If a summary is desired, print a summary of the output event count and the conditions imposed
    if summary:

        print('The output catalog retained %d out of %d events from the hypoi.pha input file.' % (
        len(catalog), total_headers))
        print('The following are the conditions used:')
        if time_lim is not None:
            print('- Temporal limits (in UTC): %s to %s' % (
            time_lim[0].strftime('%y/%m/%dT%H:%M:%S'), time_lim[1].strftime('%y/%m/%dT%H:%M:%S')))
        if radial_lim is not None:
            print('- Radial limit of %.2f km from coordinates (%f,%f)' % (radial_lim[2], radial_lim[0], radial_lim[1]))
        if bbox_lim is not None:
            print('- Bounding box from (%f,%f) to (%f,%f)' % (bbox_lim[0], bbox_lim[1], bbox_lim[2], bbox_lim[3]))
        if nsta_min is not None:
            print('- Minimum of %d stations with valid picks' % nsta_min)
        if nobs_min is not None:
            print('- Minimum number of %d valid picks' % nobs_min)
        if rmse_max is not None:
            print('- Maximum allowed travel time RMSE of %.2f s' % rmse_max)
        if gap_max is not None:
            print('- Maximum azimuthal gap of %.2f degrees' % gap_max)
        if dist2sta_max is not None:
            print('- Maximum event-to-nearest-picked-station distance of %.2f km')
        if eh_max is not None:
            print('- Maximum horizontal error of %.2f km')
        if ev_max is not None:
            print('- Maximum vertical error of %.2f km')
        if mag_req is not None:
            print('- Events must have a valid magnitude assigned')
        if depth_req is not None:
            print('- Events must have a valid depth assigned')
        if num_P_min is not None:
            print('- Minimum of %d P picks' % num_P_min)
        if num_S_min is not None:
            print('- Minimum of %d S picks' % num_S_min)
        if channel_convention == True:
            print('- Retained picks must conform to their channel convention (P on vertical, S on horizontal)')

    # End timer for code
    print('Total time taken: %.2f s' % (time.time() - time_start))

    return catalog


#%% [download_catalog] function to download catalog along with its phase data using libcomcat
def download_catalog(start_time,end_time,contributor,latitude,longitude,max_radius,max_depth,review_status='reviewed',verbose=False):

    from libcomcat.search import search
    import requests
    from libcomcat.utils import HEADERS, TIMEOUT
    from obspy.io.quakeml.core import Unpickler
    from obspy import Catalog

    # Pull earthquakes from libcomcat
    earthquakes = search(starttime=start_time.datetime,
                         endtime=end_time.datetime,
                         contributor=contributor.lower(),
                         latitude=latitude,
                         longitude=longitude,
                         maxradiuskm=max_radius,
                         maxdepth=max_depth,
                         reviewstatus=review_status)

    print('Found %d earthquakes. Processing phases...' % (len(earthquakes)))

    # Initialize catalog object to populate
    catalog = Catalog()

    # Loop over earthquakes to get phase products
    for i,eq in enumerate(earthquakes):
        if verbose:
            print('Event #%d: %s' % (i+1,eq))
        detail = eq.getDetailEvent()
        phasedata = detail.getProducts('phase-data', source=str(eq)[0:2])[0]
        quakeurl = phasedata.getContentURL('quakeml.xml')
        response = requests.get(quakeurl, timeout=TIMEOUT, headers=HEADERS)
        data = response.text.encode('utf-8')
        unpickler = Unpickler()
        catalog += unpickler.loads(data)

    return catalog


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
def remove_bad_traces(st,max_zeros=100,npts_threshold=100):

    # Recommended value: max_zeros=100 for 50Hz daylong traces

    # Import dependencies
    import numpy as np

    # Check traces
    for tr in st:
        num_zeros = np.sum(tr.data.data==0)
        if num_zeros > max_zeros or tr.stats.npts < npts_threshold:
            st.remove(tr)

    # Return filtered stream object
    return st

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
            nsta = len(family[j].chans)
            detection_value = abs(family[j].detect_val)
            detection_list.append(abs(detection_value)/nsta)
            threshold_value = family[j].threshold
            threshold_list.append(abs(threshold_value)/nsta)

    # Plot distribution of detections as histogram and save
    detection_array = np.array(detection_list) - np.array(threshold_list)
    detection_floor = np.floor(min(detection_array))
    detection_ceil = np.ceil(max(detection_array))
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.hist(detection_array, bins=np.arange(detection_floor, detection_ceil, 0.01), color='teal', edgecolor='black')
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
def read_trace(data_path,station,channel,starttime,endtime,tolerance=4e4):

    # Import dependencies
    import glob
    import numpy as np
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
    data_filename = data_path + station + '.' + channel + '.' + str(start_year) + ':' + f'{start_julday:03}' + ':*'
    data_filenames.append(data_filename)

    # Add another day if starttime and endtime are on different data files
    if start_year != end_year or start_julday != end_julday:
        data_filename = data_path + station + '.' + channel + '.' + str(end_year) + ':' + f'{end_julday:03}' + ':*'
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
    tr_split = tr_unmerged.copy()
    for tr in tr_split:
        tr.stats.sampling_rate = np.round(tr.stats.sampling_rate)
    tr_split.trim(starttime=starttime, endtime=endtime)
    tr_split = remove_boxcars(tr_split, tolerance)
    tr_split = tr_split.detrend("simple")
    tr_merged = tr_split.merge()

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
def prepare_catalog_stream(data_path,catalog,resampling_frequency,tolerance,max_zeros=100,npts_threshold=100):

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
            data_filename = data_path + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{pick_julday:03}' + ':*'
            data_filenames.append(data_filename)

            # Add next day if pick occurs in the first minute of the day
            if pick.time.hour == 0 and pick.time.minute <= 1:
                if pick.time.julday == 1:  # special case for first day of the year
                    data_filename = data_path + sta + '.' + chan + '.' + str(pick_year - 1) + ':365:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_path + sta + '.' + chan + '.' + str(
                        pick_year) + ':' + f'{(pick_julday - 1):03}' + ':*'
                    data_filenames.append(data_filename)

            # Add previous day if pick occurs in the last minute of the day
            if pick.time.hour == 23 and pick.time.minute >= 59:
                if pick.time.julday == 365:  # special case for last day of the year
                    data_filename = data_path + sta + '.' + chan + '.' + str(pick_year + 1) + ':001:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_path + sta + '.' + chan + '.' + str(
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
    stream = stream.merge()

    # Return stream object
    return stream

# [prepare_stream_dict] prepare stream dictionary pertaining to input catalog
def prepare_stream_dict(catalog,pre_pick,length,local=False,client_name="IRIS",data_path=None):

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
                tr = read_trace(data_path, station, channel, starttime, endtime, tolerance=4e4)
            else:
                tr = client.get_waveforms(network, station, '*', channel, starttime, endtime)

            # Add trace to stream object
            event_st += tr

        # Split stream and append by alongside event id
        event_st = event_st.split()
        stream_tuples.append((event_id,event_st))

    # Convert tuple object to dict object
    stream_dict = dict(stream_tuples)

    # Return stream dictionary
    return stream_dict

# [adjust_weights] adjust weights within a dt.cc file and write/append to another file
def adjust_weights(dtcc_filepath,target_filepath,dt_cap=None,min_link=0,append=False,weight_func=None):

    # Import dependencies
    import os
    import numpy as np

    # Open dt.cc file to read line by line
    dtcc_textfile = open(dtcc_filepath, 'r')
    dtcc_lines = dtcc_textfile.readlines()

    # Initialize list to store target file lines
    target_lines = []

    # Loop over dt.cc file lines
    for dtcc_line in dtcc_lines:

        # If the dt.cc line is an event pair header
        if dtcc_line[0] == '#' or len(dtcc_line) != 25:

            # Append dt.cc line directly
            target_lines.append(dtcc_line)

        # If the dt.cc line is not an event pair header, it stores a channel's dt and weight
        else:

            # Check for dt_cap if defined:
            if dt_cap:
                dt_str = dtcc_line.split()[1]
                dt = float(dt_str)
                if abs(dt) > dt_cap:
                    continue

            # Extract old weight and apply sqrt
            old_weight_str = dtcc_line.split()[2]
            old_weight = float(old_weight_str)
            if weight_func and callable(weight_func):
                new_weight = weight_func(old_weight)
            elif weight_func:
                raise ValueError('weight_func input is not a callable function!')
            else:
                new_weight = old_weight
            new_weight_str = '%0.4f' % new_weight

            # Construct a new line and append
            target_line = dtcc_line.replace(old_weight_str,new_weight_str)
            target_lines.append(target_line)

    # Remove event pairs that have less than min_link observations
    header_indices = np.flatnonzero([target_line[0] == '#' for target_line in target_lines])
    num_link = np.diff(np.append(header_indices, len(target_lines)-1))
    failed_pair_indices = header_indices[np.flatnonzero(num_link < min_link)]
    for failed_pair_index in failed_pair_indices:
        target_lines[failed_pair_index] = ''
        if failed_pair_index != (len(target_lines) - 1):
            i = 1
            while target_lines[failed_pair_index+i][0] != '#' and (failed_pair_index+i) != (len(target_lines)-1):
                target_lines[failed_pair_index+i] = ''
                i += 1
    target_lines = [target_line for target_line in target_lines if target_line != '']

    # If the target filepath does not exist
    if not os.path.exists(target_filepath):

        # Write target file
        target_textfile = open(target_filepath, 'w')
        target_textfile.writelines(target_lines)

    # If the target filepath exists
    else:

        # And if we want to append the information to the target file
        if append:

            # Open file in append mode and write
            target_textfile = open(target_filepath, 'a')
            #target_lines.insert(0,'\n')
            target_textfile.writelines(target_lines)

        # If we choose not to append
        else:

            # Delete the old filepath and write target file
            os.remove(target_filepath)
            target_textfile = open(target_filepath, 'w')
            target_textfile.writelines(target_lines)

# [calculate_catalog_FI] Calculate the FIs of events in a catalog using a reference sta-chan.
def calculate_catalog_FI(catalog, data_path, reference_station, reference_channel, min_match, resampling_frequency,
                         prepick, length, lowcut, highcut, filomin=1, filomax=2.5, fiupmin=5, fiupmax=10, tolerance=4e4,
                         verbose=False, histogram=False):

    # Import all dependencies
    import glob
    import numpy as np
    from itertools import compress
    from obspy import Stream, Catalog, read, UTCDateTime
    from obspy.core.event import Comment
    from scipy.fftpack import fft, fftfreq
    from toolbox import prepare_catalog_stream
    import matplotlib.pyplot as plt

    # If reference_station and reference_channel are comma separated strings, convert to lists
    if type(reference_station) == type(reference_channel) == str:
        reference_station = reference_station.split(',')
        reference_channel = reference_channel.split(',')

    # Zip reference station and channels
    reference_stachans = tuple(zip(reference_station, reference_channel))

    # Create a boolean list to keep track of which events have been attempted
    tracker = np.array([False for i in range(len(catalog))])

    # Print progress count for the number of successfully calculated magnitudes
    if verbose:
        print('Count of successful FI computations:')
        count = 0

    # Initialize array of all event times
    detection_times = np.array([UTCDateTime(event.resource_id.id.split('_')[-1]) for event in catalog])

    # Loop over events in catalog
    for i, event in enumerate(catalog):

        # Check if event has had its FI calculation attempted. If True, skip to next
        if tracker[i]:
            continue
        # Otherwise check for other events occuring on the same day so that we load the local data at one go
        else:
            event_daystart = UTCDateTime(UTCDateTime(event.resource_id.id.split('_')[-1]).date)
            event_dayend = event_daystart + 86400
            sub_catalog_bool = (np.array(detection_times) > event_daystart) & (np.array(detection_times) < event_dayend)
            sub_catalog_index = np.flatnonzero(sub_catalog_bool)
            sub_catalog = Catalog(list(compress(catalog, sub_catalog_bool)))
            tracker[sub_catalog_index] = True

        # Remove picks in sub_catalog events that are not on reference stachans
        for event in sub_catalog:
            for pick in event.picks:
                if (pick.waveform_id.station_code, pick.waveform_id.channel_code) not in reference_stachans:
                    event.picks.remove(pick)

        # Load the appropriate data for sub_catalog
        master_stream = prepare_catalog_stream(data_path, sub_catalog, resampling_frequency, tolerance)
        master_stream = master_stream.trim(starttime=event_daystart, endtime=event_dayend, pad=True)

        # Now loop through sub_catalog to compute FI
        for j, event in enumerate(sub_catalog):

            # Derive a subset of the master stream that correspond to picked stachans
            event_netstachans = [
                (pick.waveform_id.network_code, pick.waveform_id.station_code, pick.waveform_id.channel_code) for pick
                in event.picks]
            stream_boolean = [(trace.stats.network, trace.stats.station, trace.stats.channel) in event_netstachans for
                              trace in master_stream]
            sub_stream = Stream(compress(master_stream, stream_boolean)).copy()

            # Now trim the sub_stream based on pick times and execute fft
            for trace in sub_stream:
                trace_netstachan = (trace.stats.network, trace.stats.station, trace.stats.channel)
                pick_index = event_netstachans.index(trace_netstachan)
                reference_pick = event.picks[pick_index]
                trace.trim(reference_pick.time - prepick - 0.05 * length,
                           reference_pick.time - prepick + length + 0.05 * length)
                # Do some stream checks before detrending and filtering
                if (type(trace.data) == np.ndarray or sum(trace.data.mask) == 0) and len(trace.data) > (
                        length * trace.stats.sampling_rate):
                    trace.detrend('simple')
                    trace.taper(max_percentage=None, max_length=0.05 * length)
                    trace.filter(type='bandpass', freqmin=lowcut, freqmax=highcut, corners=4, zerophase=True)
                    trace.trim(reference_pick.time - prepick, reference_pick.time - prepick + length)
                else:
                    sub_stream.remove(trace)

            # If sub_stream is empty or if the number of matching stachans is below the user-specified threshold
            if len(sub_stream) == 0 or len(sub_stream) < min_match:
                comment_text = 'FI=None'
                catalog[sub_catalog_index[j]].comments.append(Comment(text=comment_text))
                if verbose:
                    print('Inserting FI=None for %s. Sub stream object is empty.' % (event.resource_id.id))
                continue

            # Execute fft on data
            fft_amplitude_stack = np.zeros((sub_stream[0].data.size,))
            sample_frequency = fftfreq(sub_stream[0].data.size, d=sub_stream[0].stats.delta)
            for trace in sub_stream:
                fft_amplitude = np.abs(fft(trace.data))
                fft_amplitude_stack += fft_amplitude[0:sub_stream[0].data.size]
            fft_amplitude_average = fft_amplitude_stack * (1 / len(sub_stream))

            # Calculate FI
            FI_highband_indices = np.where((sample_frequency > fiupmin) & (sample_frequency < fiupmax))
            FI_lowband_indices = np.where((sample_frequency > filomin) & (sample_frequency < filomax))
            FI = np.log10(np.mean(fft_amplitude_average[FI_highband_indices]) / np.mean(
                fft_amplitude_average[FI_lowband_indices]))

            # Add FI to event as comment
            comment_text = 'FI=%.5f' % FI
            catalog[sub_catalog_index[j]].comments.append(Comment(text=comment_text))
            if verbose:
                count += 1
                print('%d (Index [i,j] = [%d,%d])' % (count,i,j))

    # If a histogram is desired, plot a histogram of FI
    if histogram:

        # Extract FI values from catalog
        all_FIs = []
        for event in catalog:
            FI_str = event.comments[-1].text.split('=')[1]
            if FI_str != 'None':
                all_FIs.append(float(FI_str))

        # Plot histogram
        fig, ax = plt.subplots()
        ax.grid(True)
        ax.hist(all_FIs, bins=np.arange(-1.25, 1.25, 0.05), color='teal', edgecolor='black')
        ax.set_xlim([-1.25, 1.25])
        ax.set_xlabel('Frequency Index (FI)')
        ax.set_ylabel('Number of Events')
        ax.set_title('Frequency Index for %d out of %d Events' % (len(all_FIs), len(catalog)))
        fig.show()

    # Return catalog
    return catalog

# [calculate_relative_magnitudes] function to calculate magnitudes using CC based on Schaff & Richards (2014)
# Uses template events in PEC as reference magnitudes for relocatable catalog
# This function rethresholds the detected catalog by min_cc while processing it
def calculate_relative_magnitudes(catalog, tribe, data_path, noise_window, signal_window, min_cc, min_snr,
                                  shift_len, resampling_frequency, lowcut, highcut, tolerance=4e4, use_s_picks=False,
                                  verbose=False):

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
        master_catalog_stream = prepare_catalog_stream(data_path, sub_catalog, resampling_frequency, tolerance)
        master_catalog_stream = master_catalog_stream.trim(starttime=UTCDateTime(date_pair[0]),
                                                           endtime=UTCDateTime(date_pair[0]) + 86400, pad=True)

        # Craft catalog from their templates and obtain their corresponding stream
        sub_templates_list = []
        for event in sub_catalog:
            template_index = template_names.index(event.comments[0].text.split(' ')[1])
            sub_templates_list.append(tribe[template_index].event)
        sub_templates = Catalog(sub_templates_list)
        master_template_stream = prepare_catalog_stream(data_path, sub_templates, resampling_frequency, tolerance)
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
                    trace.filter(type='bandpass', freqmin=lowcut, freqmax=highcut, corners=4, zerophase=True)
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
                    trace.filter(type='bandpass', freqmin=lowcut, freqmax=highcut, corners=4, zerophase=True)
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

# [estimate_s_picks] function to go through a catalog to add estimated s picks
def estimate_s_picks(catalog, vpvs_ratio=1.73, s_pick_component='N'):
    # Loop over events in catalog
    for event in catalog:
        # Get list of picked stations
        picked_stations = list(np.unique([p.waveform_id.station_code for p in event.picks]))
        # Loop over picked stations
        for picked_station in picked_stations:
            # Get a sub-list of this station's picks
            station_picks = [p for p in event.picks if p.waveform_id.station_code == picked_station]
            # If there is only one P-pick, make an S-pick
            if len(station_picks)==1 and station_picks[0].phase_hint == 'P':
                # Copy P pick to work with
                pick_copy = station_picks[0].copy()
                # Calculate S pick traveltime
                p_traveltime =  pick_copy.time - event.origins[0].time
                s_traveltime = vpvs_ratio * p_traveltime
                # Craft S pick and append
                pick_copy.phase_hint = 'S'
                pick_copy.waveform_id.channel_code = pick_copy.waveform_id.channel_code[:-1] + s_pick_component
                pick_copy.time = event.origins[0].time + s_traveltime
                event.picks.append(pick_copy)
    return catalog

# [clean_cc_file] removes repeated event pairs and observations in dt.cc file
# WARNING: this function works in place.
def clean_cc_file(dtcc_filepath):
    # Open dt.cc file
    with open(dtcc_filepath, 'r') as dtcc_file:
        # Read all lines and split by line break
        all_lines = dtcc_file.read()
        lines = all_lines.split('\n')
        # Initialize list storing new lines and event pairs
        new_lines = []
        event_pairs = []
        # Loop through lines and add lines to list if they are new
        for line in lines:
            # Check if line is a header line
            if line != '' and line[0] == '#':
                # Switch off phase inclusion and define both permutations of the current event pair
                include_phases = False
                event_pair = (int(line.split()[1]),int(line.split()[2]))
                event_pair_flipped = (int(line.split()[2]), int(line.split()[1]))
                # Check against current list
                if event_pair not in event_pairs and event_pair_flipped not in event_pairs:
                    # Append if unique, and switch on phase inclusion
                    new_lines.append(line)
                    event_pairs.append(event_pair)
                    include_phases = True
            # Append phase line if phase inclusion is switched on
            elif include_phases:
                new_lines.append(line)
            # Otherwise, skip line
            else:
                continue
        # Add final empty line and join by line break
        new_lines.append('')
        new_all_lines = '\n'.join(new_lines)
    # Overwrite input dt.cc file
    with open(dtcc_filepath, 'w') as dtcc_file:
        dtcc_file.write(new_all_lines)


# [loc2cat] converts .loc and .reloc into Catalog objects
def loc2cat(loc_filepath, input_catalog=None, type='loc', depth_correction=0):

    # Import all dependencies
    from obspy import Catalog, UTCDateTime
    from obspy.core.event import Event, Origin, Magnitude, Comment
    from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper

    # Get list of input event resource ids
    if input_catalog:
        all_ids = [event.resource_id.id for event in input_catalog]
        event_id_mapper = _generate_event_id_mapper(input_catalog, event_id_mapper=None)

    # Initialize output catalog
    outcat = Catalog()

    # Open loc file and read lines
    with open(loc_filepath, "r") as loc:
        lines = loc.read()
        lines = lines.split('\n')

        # Loop over lines
        for line in lines[:-1]:

            # Extract information and craft event
            if type == 'loc':
                (evid, lat, lon, dep, _, _, _, _, _, _, yr, mo, dy, hr, mm, ss, mag, cid) = line.split()
            elif type == 'reloc':
                (evid, lat, lon, dep, _, _, _, _, _, _, yr, mo, dy, hr, mm, ss, mag, _, _, _, _, _, _,
                 cid) = line.split()
                if lat == 'NaN':
                    continue
            else:
                raise ValueError('type argument is not loc/reloc!')
            yr = int(yr)
            mo = int(mo)
            dy = int(dy)
            hr = int(hr)
            mm = int(hr)
            ss = float(ss)
            if ss < 60:
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            else:
                ss -= 60
                mm += 1
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            lon = float(lon)
            lat = float(lat)
            dep = (float(dep) - depth_correction)
            mag = float(mag)
            ev = Event(origins=[Origin(time=time, longitude=lon, latitude=lat, depth=dep)],
                       magnitudes=[Magnitude(mag=mag)])

            # Find base event and copy over comments
            if input_catalog:
                old_event_id = [id for id, num in event_id_mapper.items() if num == int(evid)]
                if len(old_event_id) != 1:
                    raise ValueError('Multiple matched on event id mapper. Check event ids!')
                else:
                    old_event_id = old_event_id[0]
                    old_event = input_catalog[all_ids.index(old_event_id)]
                    ev.comments = old_event.comments
                    evid_origin_comment = 'Event ID: %s' % evid
                    ev.origins[0].comments.append(Comment(text=evid_origin_comment))

            # Append event to catalog
            outcat.append(ev)

    # Return catalog
    return outcat


# [dat2cat] converts .dat and .sel to Catalog objects
def dat2cat(dat_filepath, input_catalog=None, depth_correction=0):

    # Import all dependencies
    from obspy import Catalog, UTCDateTime
    from obspy.core.event import Event, Origin, Magnitude, Comment
    from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper

    # Get list of input event resource ids
    if input_catalog:
        all_ids = [event.resource_id.id for event in input_catalog]
        event_id_mapper = _generate_event_id_mapper(input_catalog, event_id_mapper=None)

    # Initialize output catalog
    outcat = Catalog()

    # Open loc file and read lines
    with open(dat_filepath, "r") as dat:
        lines = dat.read()
        lines = lines.split('\n')

        # Loop over lines
        for line in lines[:-1]:

            # Extract information and craft event
            (yrmody, hrmmsscs, lat, lon, dep, mag, _, _, _, evid) = line.split()
            yr = int(yrmody[0:4])
            mo = int(yrmody[4:6])
            dy = int(yrmody[6:8])
            hrmmsscs = '%08d' % int(hrmmsscs)
            hr = int(hrmmsscs[:-6])
            mm = int(hrmmsscs[-6:-4])
            ss = float(hrmmsscs[-4:-2]) + float(hrmmsscs[-2:]) / 100
            if ss < 60:
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            else:
                ss -= 60
                mm += 1
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            ev = Event(origins=[Origin(time=time, longitude=float(lon), latitude=float(lat),
                                       depth=(float(dep) - depth_correction))],
                       magnitudes=[Magnitude(mag=float(mag))])

            # Find base event and copy over resource id and comments
            if input_catalog:
                old_event_id = [id for id, num in event_id_mapper.items() if num == int(evid)]
                if len(old_event_id) != 1:
                    raise ValueError('Multiple matched on event id mapper. Check event ids!')
                else:
                    old_event_id = old_event_id[0]
                    old_event = input_catalog[all_ids.index(old_event_id)]
                    ev.comments = old_event.comments
                    evid_origin_comment = 'Event ID: %s' % evid
                    ev.origins[0].comments.append(Comment(text=evid_origin_comment))

            # Append event to catalog
            outcat.append(ev)

    # Return catalog
    return outcat


# [remove_catalog_repeats] returns a copied catalogA, but repeat events in catalogB are removed
def remove_catalog_repeats(catalogA, catalogB):

    # Parse event ids
    catalogA_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogA]
    catalogB_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogB]

    # Copy catalog
    catalogC = catalogA.copy()

    # Loop over event ids, and if there is a repeat, remove events from the copied catalog
    for i, catalogA_evid in reversed(list(enumerate(catalogA_evids))):
        if catalogA_evid in catalogB_evids:
            catalogC.events.pop(i)

    return catalogC