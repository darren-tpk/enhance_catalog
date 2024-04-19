def initialize_run(name):

    """
    Initialize run by creating nested output and data repositories.
    :param name: desired name of workflow run
    :return: N/A
    """

    # Import all dependencies
    import os

    # Create all data and output directories and subdirectories
    data_main = './data/'
    data_dir = data_main + name + '/'
    output_main = './output/'
    output_dir = output_main + name + '/'
    print('\nCreating subdirectories for workflow outputs...')
    output_subdirs = [data_main, data_dir,
                      output_main, output_dir,
                      output_dir + 'run_redpy',
                      output_dir + 'convert_redpy/',
                      output_dir + 'create_tribe/',
                      output_dir + 'scan_data/',
                      output_dir + 'relocate_catalog/']
    for output_subdir in output_subdirs:
        try:
            os.mkdir(output_subdir)
        except FileExistsError:
            print('%s already exists.' % output_subdir)
    print('All data and output subdirectories created.')

def download_data(data_destination,starttime,endtime,client,network,station,channel,location):

    """
    Downloads data across a user-defined duration from a Client, using ObsPy's st.read() and st.write().
    :param data_destination (str): miniseed output path
    :param starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): start time for desired data pull
    :param endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): end time for desired data pull
    :param client (str): Name of desired FDSN client (e.g., IRIS)
    :param network (str or list): SEED network code(s)
    :param station (str or list): SEED station code(s)
    :param channel (str or list): SEED channel code(s)
    :param location (str or list): SEED location code(s)
    :return: N/A
    """

    # Import all dependencies
    import glob
    import time
    import numpy as np
    from obspy import UTCDateTime, Stream
    from obspy.clients.fdsn import Client

    # Point to client
    client = Client(client)

    # If input for net/sta/chan/loc are str, reformat them into lists
    if type(network) == type(station) == type(channel) == type(location) == str:
        network = network.split(',')
        station = station.split(',')
        channel = channel.split(',')
        location = location.split(',')

    # Round down starttime to the nearest day and round up endtime to the nearest day
    starttime = UTCDateTime(starttime.date)
    if endtime != UTCDateTime(endtime.date):
        endtime = UTCDateTime(endtime.date) + 86400

    # Dissect duration into days
    num_days = int(np.floor((endtime - starttime) / 86400))
    print("\nCommencing data fetch...")
    time_start = time.time()

    # Commence loop over days
    for i in range(num_days):

        # Define temporal boundaries for data fetch
        t1 = starttime + (i * 86400)
        t2 = starttime + ((i + 1) * 86400)
        print('\nNow at %s...' % str(t1.date))

        # Loop over station list and extract relevant seed information
        for j in range(len(station)):

            net = network[j]
            sta = station[j]
            chan = channel[j]
            loc = location[j]

            # Check if file already exists
            julday_str = '%03d' % t1.julday
            seed_sample_filepath = data_destination + sta + '.' + chan + '.' + str(t1.year) + ':' + julday_str + ':*'
            matching_data_files = glob.glob(seed_sample_filepath)
            if len(matching_data_files) != 0:
                print('%s already has data files for %s, skipping.' % (sta, t1.date))
                continue

            # Get waveforms by querying client
            try:
                st = client.get_waveforms(net, sta, loc, chan, t1, t2)

                # Save every trace separately
                for tr in st:

                    # Craft the seed file name
                    trace_year = str(tr.stats.starttime.year)
                    trace_julday = str(tr.stats.starttime.julday)
                    trace_time = str(tr.stats.starttime.time)[0:8]
                    trace_datetime_str = ":".join([trace_year, str(trace_julday).zfill(3), trace_time])
                    seed_filename = tr.stats.station + '.' + tr.stats.channel + '.' + trace_datetime_str

                    # Write to seed file
                    seed_filepath = data_destination + seed_filename
                    tr.write(seed_filepath, format="MSEED")

                    # Print success
                    print('%s successfully saved.' % seed_filename)

            # continue if waveform retrieval fails
            except:
                print('%s.%s failed.' % (sta,chan))
                continue

        # Print progression
        time_current = time.time()
        print('Data collection for the day complete. Elapsed time: %.2f hours\n' % ((time_current - time_start) / 3600))

    # Conclude process
    time_end = time.time()
    print('\nData collection complete. Time taken: %.2f hours' % ((time_end - time_start) / 3600))

def run_redpy(run_title,
              output_destination,
              data_path,
              starttime,
              endtime,
              network,
              station,
              channel,
              location,
              stalats,
              stalons,
              samprate,
              fmin,
              fmax,
              nstaC,
              lwin,
              swin,
              trigon,
              trigoff,
              winlen,
              cmin,
              ncor,
              minorph,
              maxorph,
              nsec=3600,
              filomin=1,
              filomax=2.5,
              fiupmin=5,
              fiupmax=10,
              fispanlow=-0.5,
              fispanhigh=0.5,
              trigalg='recstalta',
              offset=0,
              plotformat='eqrate,fi,occurrence+occurrencefi,longevity',
              printsta=0,
              minplot=5,
              dybin=1,
              hrbin=1.,
              occurbin=24,
              recbin=24,
              fixedheight=False,
              recplot=14.,
              mminplot=0,
              mhrbin=1.,
              mrecbin=1.,
              mrecplot=30.,
              verbosecatalog=False,
              anotfile='annotation.csv',
              amplims='global',
              checkComCat=False,
              serr=5.,
              locdeg=0.5,
              regdeg=2.,
              regmag=2.5,
              telemag=4.5,
              matchMax=0,
              kurtwin=5.,
              kurtmax=80.,
              kurtfmax=150.,
              oratiomax=0.15,
              telefi=-1.0,
              teleok=2):

    """
    Runs REDPy on pre-downloaded data across a user-defined date range
    :param run_title: title of the REDPy (use unique names for unique runs to prevent overwriting)
    :param output_destination (str): output directory for REDPy
    :param data_path (str): path to where relevant miniseed files can be found
    :param starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): start time for REDPy scan
    :param endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): end time for REDPy scan
    :param network (str or list): SEED network code(s)
    :param station (str or list): SEED station code(s)
    :param channel (str or list): SEED channel code(s)
    :param location (str or list): SEED location code(s)
    :param stalats (str or list or typle): list (or tuple or comma separated string) of station latitudes
    :param stalons (str or list or tuple): list (or tuple or comma separated string) of station latitudes
    :param samprate (int): sampling rate to store waveforms in (data will be resampled if there is a mismatch)
    :param fmin (float): bandpass filter lower bound (Hz)
    :param fmax (float): bandpass filter upper bound (Hz)
    :param nstaC (int): number of stations where a coincident trigger must be observed on to make a detection
    :param lwin (float): long window for STA/LTA (s)
    :param swin (float): short window for STA/LTA (s)
    :param trigon (float): STA/LTA ratio for trigger to turn on
    :param trigoff (float): STA/LTA ratio for trigger to turn off
    :param winlen (int): cross-correlation window length in samples (2^n is preferred)
    :param cmin (float): minimum cross-correlation coefficient to be considered a repeater
    :param ncor (int): number of stations where cmin must be observed to be counted as a repeater
    :param minorph (float): amount of days to keep orphans in the queue when it just triggers above threshold
    :param maxorph (float): amount of days to keep orphans in the queue when it triggers way above threshold
    :param nsec (float): maximum data length (in seconds) to process at once (defaults to 3600 s)
    :param filomin (float): frequency index lower bound minimum (defaults to 1.0 Hz)
    :param filomax (float): frequency index lower bound maximum (defaults to 2.5 Hz)
    :param fiupmin (float): frequency index upper bound minimum (defaults to 5.0 Hz)
    :param fiupmax (float): frequency index upper bound minimum (defaults to 10.0 Hz)
    :param fispanlow (float): FI span minimum for plot colorbar (dimensionless, defaults to -0.5)
    :param fispanhigh (float): FI span maximum for plot colorbar (dimensionless, defaults to +0.5)
    :param trigalg (str): chosen STA/LTA algorithm supported by ObsPy's network coincidence trigger (defaults to 'recstalta')
    :param offset (float): optional time offset (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg) (defaults to 0)
    :param plotformat (str): list and order of plots to be included (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg)
    :param printsta (int): index of the station that is shown on the main plot, and used for plotting cores and amplitudes (defaults to 0)
    :param minplot (int): minimum number of members within a REDPy family that is needed before it appears on the timeline plot (defaults to 5)
    :param dybin (float): width (in days) of bins for the histogram subplot (defaults to 1)
    :param hrbin (float): width (in hours) of bins for histogram subplot on 'recent' timeline (defaults to 1)
    :param occurbin (float): width (in hours) in occurrence plot (defaults to 24)
    :param recbin (float): width (in hours) in occurrence plot (defaults to 24)
    :param fixedheight (bool): option for the occurrence plot to have the same height as the other subplots (defaults to False)
    :param recplot (float): number of days prior to last repeater to show on 'recent' timeline (defaults to 14)
    :param mminplot (float): settings for separate meta view customizations (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg)
    :param mhrbin (float): settings for separate meta view customizations (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg)
    :param mrecbin (float): settings for separate meta view customizations (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg)
    :param mrecplot (float): settings for separate meta view customizations (see https://github.com/ahotovec/REDPy/blob/master/settings.cfg)
    :param verbosecatalog (bool): option to print a more verbose version of the catalog (defaults to False)
    :param amplims (str): use 'global' or 'family' to define amplitude plot limits (defaults to 'global')
    :param checkComCat (bool): if `True`, search ANSS Comprehensive Catalog for matches (defaults to `False`)
    :param serr (float): seconds of error allowable to be considered a match (defaults to 5 s)
    :param locdeg (float): distance in degrees to search for local matches (defaults to 0.5 degrees)
    :param regdeg (float): distance in degrees to search for regional matches (defaults to 5 degrees)
    :param regmag (float): minimum magnitude to allow for regional matches (defaults to M2.5)
    :param telemag (float): minimum magnitude to allow for teleseismic matches (defaults to M4.5)
    :param matchMax (int): maximum number of events to try to match (defaults to 0, which matches all)
    :param kurtwin (float): minimum kurtosis of waveform in a window for spike removal (defaults to 5)
    :param kurtmax (float): maximum kurtosis of waveform in a window for spike removal (defaults to 80)
    :param oratiomax (float): maximum ratio of outliers to total number of datapoints in trace (defaults to 0.15)
    :param telefi (float): frequency index minimum threshold for teleseisms (defaults to -1.0)
    :param teleok (float): number of stations allowed to pass that exceed telefi threshold (defaults to 2)
    :return: N/A
    """

    # Import all dependencies
    import subprocess
    import sys
    from obspy import UTCDateTime

    # If input net/sta/chan/loc/stalats/stalons are lists, merge into comma separated string
    if type(network) == type(station) == type(channel) == type(location) == type(stalats) == type(stalons) == list:
        network = ','.join(network)
        station = ','.join(station)
        channel = ','.join(channel)
        location = ','.join(location)
    # Reformat location if not process (convert '--' to blank)
    location = location.replace('--','')

    # If stalats/stalons are lists/tuples, merge into comma separated string
    if type(stalats) == type(stalons) == list:
        stalats = ','.join([str(f) for f in stalats])
        stalons = ','.join([str(f) for f in stalons])
    elif type(stalats) == type(stalons) == tuple:
        stalats = ','.join(map(str, stalats))
        stalons = ','.join(map(str, stalons))

    # Round down starttime to the nearest day and round up endtime to the nearest day
    starttime = UTCDateTime(starttime.date)
    if endtime != UTCDateTime(endtime.date):
        endtime = UTCDateTime(endtime.date) + 86400

    # Craft groupname and filename
    group_name = "".join(run_title.split())
    cfg_name = group_name+'.cfg'
    h5_name = group_name+'.h5'

    # Craft config file
    config_text = ['[Settings]',
                   '\n## RUN PARAMETERS ##',
                   'title='+run_title,
                   'outputPath='+output_destination,
                   'groupName='+group_name,
                   'filename='+output_destination+h5_name,
                   'minorph='+str(minorph),
                   'maxorph='+str(maxorph),
                   'nsec='+str(nsec),
                   '\n## STATION PARAMETERS ##',
                   'nsta='+str(len(station.split(','))),
                   'station='+station,
                   'channel='+channel,
                   'network='+network,
                   'location='+location,
                   'samprate='+str(samprate),
                   'fmin='+str(fmin),
                   'fmax='+str(fmax),
                   'filomin='+str(filomin),
                   'filomax='+str(filomax),
                   'fiupmin='+str(fiupmin),
                   'fiupmax='+str(fiupmax),
                   'fispanlow='+str(fispanlow),
                   'fispanhigh='+str(fispanhigh),
                   '\n## DATA SOURCE ##',
                   'server=file',
                   'searchdir='+data_path,
                   'filepattern=*.*.*:*:*:*:*',
                   '\n## TRIGGERING SETTINGS ##',
                   'trigalg='+trigalg,
                   'nstaC='+str(nstaC),
                   'lwin='+str(lwin),
                   'swin='+str(swin),
                   'trigon='+str(trigon),
                   'trigoff='+str(trigoff),
                   'offset='+str(offset),
                   '\n## CROSS-CORRELATION PARAMETERS ##',
                   'winlen='+str(winlen),
                   'cmin='+str(cmin),
                   'ncor='+str(ncor),
                   '\n## PLOTTING PARAMETERS ##',
                   'plotformat='+plotformat,
                   'printsta='+str(printsta),
                   'minplot='+str(minplot),
                   'dybin='+str(dybin),
                   'hrbin='+str(hrbin),
                   'occurbin='+str(occurbin),
                   'recbin='+str(recbin),
                   'fixedheight='+str(fixedheight),
                   'recplot='+str(recplot),
                   'mminplot='+str(mminplot),
                   'mhrbin='+str(mhrbin),
                   'mrecbin='+str(mrecbin),
                   'mrecplot='+str(mrecplot),
                   'verbosecatalog='+str(verbosecatalog),
                   'amplims='+amplims,
                   '\n## CHECK COMCAT FOR MATCHES ##',
                   'checkComCat='+str(checkComCat),
                   'stalats='+stalats,
                   'stalons='+stalons,
                   'serr='+str(serr),
                   'locdeg='+str(locdeg),
                   'regdeg='+str(regdeg),
                   'regmag='+str(regmag),
                   'telemag='+str(telemag),
                   'matchMax='+str(matchMax),
                   '\n## AUTOMATED SPIKE AND TELESEISM REMOVAL ##',
                   'kurtwin='+str(kurtwin),
                   'kurtmax='+str(kurtmax),
                   'kurtfmax='+str(kurtfmax),
                   'oratiomax='+str(oratiomax),
                   'telefi='+str(telefi),
                   'teleok='+str(teleok)]
    with open(output_destination + cfg_name, 'w') as f:
        f.write('\n'.join(config_text))

    # Call initialize command
    initialize_command = '%s ./redpy/initialize.py -v -c %s' % (sys.executable, output_destination + cfg_name)
    subprocess.call(initialize_command, shell=True)

    # Call backfill command
    initialize_command = '%s ./redpy/backfill.py -v -c %s -s %s -e %s' % (sys.executable,
                                                                          output_destination + cfg_name,
                                                                          starttime.strftime('%Y-%m-%d'),
                                                                          endtime.strftime('%Y-%m-%d'))
    subprocess.call(initialize_command, shell=True)

def convert_redpy(analyst_catalog,
                  redpy_output_path,
                  convert_redpy_output_dir,
                  data_path,
                  redpy_station,
                  redpy_channel,
                  fmin,
                  fmax,
                  nstaC,
                  lwin,
                  swin,
                  trigon,
                  trigoff,
                  max_dt=4,
                  tolerance=4e4,
                  add_redpy_pick_to_associated=True,
                  add_campaign_pick_to_associated=False,
                  campaign_starttime=None,
                  campaign_endtime=None,
                  campaign_station=None,
                  campaign_channel=None):

    """
    :param analyst_catalog (:class:`~obspy.core.event.Catalog`): analyst-derived catalog to reference
    :param redpy_output_path (str): output directory from run_redpy step
    :param convert_redpy_output_dir (str): output directory for the convert_redpy step
    :param data_path (str): path to where relevant miniseed files can be found
    :param redpy_station (str or list): SEED station code(s) used in REDPy
    :param redpy_channel (str or list): SEED channel code(s) used in REDPy
    :param fmin (float): bandpass filter lower bound (Hz)
    :param fmax (float): bandpass filter upper bound (Hz)
    :param nstaC (int): number of stations where a coincident trigger must be observed on to make a detection
    :param lwin (float): long window for STA/LTA (s)
    :param swin (float): short window for STA/LTA (s)
    :param trigon (float): STA/LTA ratio for trigger to turn on
    :param trigoff (float): STA/LTA ratio for trigger to turn off
    :param max_dt (float): maximum time difference between REDPy detections and AVO events allowed (s)
    :param tolerance (float): factor to median tolerance for boxcar removal from data (as a factor to median)
    :param add_redpy_pick_to_associated (bool): if `True`, translate analyst picks on REDPy stations to associated REDPy detections
    :param add_campaign_pick_to_associated (bool): if `True`, use STA/LTA to make coincident pick on campaign stations
    :param campaign_starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): start time of available campaign/backfilled data
    :param campaign_endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): end time of available campaign/backfilled data
    :param campaign_station (str or list): SEED station code(s) available from campaign/backfilled data
    :param campaign_channel (str or list): SEED channel code(s) available from campaign/backfilled data
    :return: N/A
    """

    # Import all dependencies
    import time
    import pandas as pd
    import numpy as np
    from obspy import UTCDateTime, Catalog, Stream
    from obspy.core.event import Event, Origin, Comment, WaveformStreamID, Pick
    from obspy.signal.trigger import coincidence_trigger, classic_sta_lta, trigger_onset
    from toolbox import writer, pull_cores, read_trace

    # If input for net/sta/chan/loc are str, reformat them into lists
    if type(redpy_station) == type(redpy_channel) == str:
        redpy_station = redpy_station.split(',')
        redpy_channel = redpy_channel.split(',')
    if add_campaign_pick_to_associated and (type(campaign_station) == type(campaign_channel) == list):
        campaign_station = campaign_station.split(',')
        campaign_channel = campaign_channel.split(',')

    # Read redpy text files using pandas
    redpy_detections = pd.read_csv(redpy_output_path + 'catalog.txt', sep=' ', names=['Cluster', 'DateTime'])
    redpy_cores = pd.read_csv(redpy_output_path + 'cores.txt', sep=' ', names=['Cluster', 'DateTime'])

    # Prepare core and analyst-derived catalog datetimes as arrays
    core_event_times = np.array([UTCDateTime(t) for t in redpy_cores.DateTime])
    analyst_event_times = np.array([analyst_event.origins[0].time for analyst_event in analyst_catalog])

    # Initialize catalog object and lists before commencing loop
    redpy_catalog = Catalog()  # empty catalog to populate
    associated_cluster_list = []  # store unique cluster numbers that have associated catalog events
    unmatched_indices_redpy = list(range(len(analyst_catalog)))  # list of avo event indices. matched indices are removed
    unmatched_indices_core = list(range(len(analyst_catalog)))  # list of avo event indices. matched indices are removed

    # Loop through redpy detection list
    print('\nCreating ObsPy catalog for REDPy results...')
    time_start = time.time()
    for i in range(len(redpy_detections)):

        # Extract values
        cluster = redpy_detections.Cluster[i]
        detection_time = UTCDateTime(redpy_detections.DateTime[i])

        # Check if event is a core
        if min(abs(core_event_times - detection_time)) == 0:
            event_tag = 'Cluster ' + str(cluster) + ' core event;'
            redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text=event_tag)])])

        # Otherwise it is a cluster event
        else:
            event_tag = 'Cluster ' + str(cluster) + ' event;'
            redpy_event = Event(origins=[Origin(time=detection_time, comments=[Comment(text=event_tag)])])

        # Check if event is part of pre-existing catalog, using an inter-event tolerance defined by user
        if min(abs(analyst_event_times - detection_time)) < max_dt:

            # Find the closest event in time from the pre-existing catalog
            analyst_index = np.argmin(abs(analyst_event_times - detection_time))
            analyst_event = analyst_catalog[analyst_index]

            # Use the closest event's information to fill up event object
            redpy_event = analyst_event.copy()
            full_event_tag = event_tag + ' in analyst catalog, REDPy detection time = %s' % str(detection_time)
            redpy_event.origins[0].comments.append(Comment(text=full_event_tag))

            # Add cluster number to valid clusters list
            associated_cluster_list.append(cluster)

            # Remove avo index from unmatched list if it still exists
            if analyst_index in unmatched_indices_redpy:
                unmatched_indices_redpy.remove(analyst_index)
            if analyst_index in unmatched_indices_core and event_tag.split(' ')[2] == 'core':
                unmatched_indices_core.remove(analyst_index)

        # If event is not part of the analyst, we tag it as so
        else:
            redpy_event.origins[0].comments[0].text += ' not in analyst catalog'

        # Finally we append the event
        redpy_catalog.append(redpy_event)

    # Write the redpy catalog to an xml file
    writer(convert_redpy_output_dir + 'redpy_catalog.xml', redpy_catalog)

    # Get unique list of AVO-associated clusters and non-associated clusters
    associated_clusters = list(np.unique(np.array(associated_cluster_list)))
    unassociated_clusters = [cluster for cluster in list(np.unique(np.array(redpy_detections.Cluster)))
                             if cluster not in associated_clusters]

    # Write them as npy files
    np.save(convert_redpy_output_dir + 'unassociated_clusters.npy', unassociated_clusters)
    np.save(convert_redpy_output_dir + 'associated_clusters.npy', associated_clusters)
    np.save(convert_redpy_output_dir + 'unmatched_indices_redpy.npy', unmatched_indices_redpy)
    np.save(convert_redpy_output_dir + 'unmatched_indices_core.npy', unmatched_indices_core)

    # Also generate unmatched analyst catalogs (that can be used as templates later)
    unmatched_analyst_events_redpy = Catalog() # analyst events that don't match redpy catalog
    unmatched_analyst_events_core = Catalog() # analyst events that don't match redpy cores
    for j, analyst_event in enumerate(analyst_catalog):
        if j in unmatched_indices_redpy:
            unmatched_analyst_events_redpy += analyst_event
        if j in unmatched_indices_core:
            unmatched_analyst_events_core += analyst_event

    # Write out unmatched analyst catalog to .xml file
    writer(convert_redpy_output_dir + 'unmatched_analyst_events_redpy.xml', unmatched_analyst_events_redpy)
    writer(convert_redpy_output_dir + 'unmatched_analyst_events_core.xml', unmatched_analyst_events_core)

    # Conclude process
    time_stop = time.time()
    print('Catalog object created, processing time: %.2f s' % (time_stop-time_start))

    # Call pull_cores function
    core_catalog = pull_cores(redpy_catalog)

    # Write core catalog to .xml file
    writer(convert_redpy_output_dir + 'core_catalog.xml', core_catalog)

    # Commence process
    print('\nMaking picks for non-associated cluster events...')
    time_start = time.time()

    # Revisit redpy catalog to add picks using single channel STA/LTA
    # Loop through unassociated clusters
    for unassociated_cluster in unassociated_clusters:

        # Find all detections within unassociated cluster
        contribution_list = list(np.where(np.array(redpy_detections.Cluster) == unassociated_cluster)[0])

        # Loop through these unassociated detections to add picks
        for contribution_index in contribution_list:

            # Determine time difference to offset pick times
            contribution_time = redpy_catalog[contribution_index].origins[0].time

            # Retrieve time limits for data fetch (+/- 12s window)
            starttime = contribution_time - 12
            endtime = contribution_time + 12

            # Gather data from local machine in a +/- 12s window, filter, and taper
            stream = Stream()
            for k in range(len(redpy_station)):
                station_tr = read_trace(data_path=data_path, station=redpy_station[k], channel=redpy_channel[k],
                                        starttime=starttime, endtime=endtime, tolerance=tolerance)
                stream = stream + station_tr
            stream = stream.split()
            stream = stream.filter('bandpass',freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
            stream = stream.taper(0.05, type='hann', max_length=(0.75*1024/100))  # [HARD CODED]
            stream = stream.merge(method=1, fill_value=0)
            stream = stream.trim(starttime=starttime, endtime=endtime)

            # Use coincidence trigger to get a pick time estimate
            try:
                coin_trigger = []
                coin_trigger = coincidence_trigger('recstalta', trigon, trigoff, stream, nstaC, sta=swin, lta=lwin, details=True)
            # If it fails (usually due to data being too gappy), move to next event
            except:
                continue

            # If there are no coincidence triggers, move to next event
            if not coin_trigger:
                continue

            # Otherwise, continue to single channel STA/LTA pick making
            else:

                # Extract coincidence trigger time
                coin_trigger_time = coin_trigger[0]["time"]

                # For each channel,
                for tr in stream:

                    # Calculate the value of the characteristic function
                    sampling_rate = tr.stats.sampling_rate
                    cft = classic_sta_lta(tr.data, int(swin * sampling_rate), int(lwin * sampling_rate))

                    # Obtain trigger limits
                    trigger_limits = np.array(trigger_onset(cft, trigon, trigoff))

                    # If there exists some trigger limits
                    if trigger_limits.size != 0:

                        # Convert to UTCDateTime and find the trigger-on time closest to the coincidence trigger
                        trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
                        pick_time = trigger_on[np.argmin(abs(trigger_on-coin_trigger_time))]

                        # Craft waveform stream ID
                        tr_id = tr.id.split('.')
                        waveform_id = WaveformStreamID(tr_id[0],tr_id[1],tr_id[2],tr_id[3])

                        # Create ObsPy pick and ObsPy arrival objects and add to redpy_catalog
                        add_pick = Pick(time=pick_time,waveform_id=waveform_id,phase_hint='P',comments=[Comment(text='derived from REDPy STA/LTA')])
                        redpy_catalog[contribution_index].picks.append(add_pick)

    # For associated clusters, give their contribution list redpy station picks based on their closest event in time
    if add_redpy_pick_to_associated:

        # Commence process
        print('Translating redpy picks for associated cluster events...')
        time_start = time.time()

        # Loop through associated clusters
        for associated_cluster in associated_clusters:

            # Find all detections within associated cluster
            contribution_list = list(np.where(np.array(redpy_detections.Cluster) == associated_cluster)[0])

            # Find indices of events with picks and without picks within cluster
            contribution_indices_with_picks = []
            contribution_indices_without_picks = []
            for contribution_index in contribution_list:
                if redpy_catalog[contribution_index].picks == []:
                    contribution_indices_without_picks.append(contribution_index)
                else:
                    contribution_indices_with_picks.append(contribution_index)

            # Loop over event indices that need picks
            for contribution_index in contribution_indices_without_picks:

                # Find the pick sharer event and use its REDPy detection time to determine pick offset
                contribution_time = redpy_catalog[contribution_index].origins[0].time
                contribution_time_diffs = [(UTCDateTime(redpy_catalog[ind].origins[0].comments[0].text.split()[-1])
                                            - contribution_time) for ind in contribution_indices_with_picks]
                contribution_time_diffs_abs = [abs(diff) for diff in contribution_time_diffs]
                pick_sharer_event = redpy_catalog[contribution_indices_with_picks[np.argmin(contribution_time_diffs_abs)]]
                contribution_time_diff = contribution_time_diffs[np.argmin(contribution_time_diffs_abs)]

                # Loop through picks and give picks if it is from a redpy station (translating it in time)
                for pick in pick_sharer_event.picks:
                    if pick.waveform_id.station_code in redpy_station:
                        add_pick = pick.copy()
                        add_pick.time -= contribution_time_diff
                        redpy_catalog[contribution_index].picks.append(add_pick)

    if add_campaign_pick_to_associated:

        # Commence process
        print('Adding campaign picks for associated cluster events...')
        time_start = time.time()

        # Define some coincidence_trigger arguments

        # Loop through associated clusters
        for associated_cluster in associated_clusters:

            # Find all detections within unassociated cluster
            contribution_list = list(np.where(np.array(redpy_detections.Cluster) == associated_cluster)[0])

            # Loop through these unassociated detections to add picks
            for contribution_index in contribution_list:

                # Extract contribution time
                contribution_time = redpy_catalog[contribution_index].origins[0].time

                # Retrieve time limits for data fetch (+/- 12s window)
                starttime = contribution_time - 12
                endtime = contribution_time + 12

                # if the contribution is from before the campaign was installed, we skip the pick-adding for the event
                if contribution_time < campaign_starttime or contribution_time > campaign_endtime:
                    continue

                # Gather data from local machine in a +/- 12s window, filter, and taper
                stream = Stream()
                for k in range(len(campaign_station)):
                    station_tr = read_trace(data_path=data_path, station=campaign_station[k], channel=campaign_channel[k],
                                            starttime=starttime, endtime=endtime, tolerance=tolerance)
                    stream = stream + station_tr
                stream = stream.split()
                stream = stream.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
                stream = stream.taper(0.05, type='hann', max_length=(0.75 * 1024 / 100))
                stream = stream.merge(method=1, fill_value=0)
                stream = stream.trim(starttime=starttime, endtime=endtime)

                # Use coincidence trigger to get a pick time estimate
                try:
                    coin_trigger = []
                    coin_trigger = coincidence_trigger('recstalta', trigon, trigoff, stream, nstaC, sta=swin, lta=lwin, details=True)
                # If it fails (usually due to data being too gappy), move to next event
                except:
                    continue

                # If there are no coincidence triggers, move to next event
                if not coin_trigger:
                    continue

                # Otherwise, continue to single channel STA/LTA pick making
                else:

                    # Extract coincidence trigger time
                    coin_trigger_time = coin_trigger[0]["time"]

                    # Remove traces that belong to channels that already have been picked
                    existing_station_list = [pick.waveform_id.station_code for pick in redpy_catalog[contribution_index].picks]
                    for tr in stream:
                        if tr.stats.station in existing_station_list:
                            stream.remove(tr)

                    # For each channel,
                    for tr in stream:

                        # Calculate the value of the characteristic function
                        sampling_rate = tr.stats.sampling_rate
                        cft = classic_sta_lta(tr.data, int(swin * sampling_rate), int(lwin * sampling_rate))
                        # To plot function, use: plot_trigger(tr, cft, 3, 2)

                        # Obtain trigger limits
                        trigger_limits = np.array(trigger_onset(cft, trigon, trigoff))

                        # If there exists some trigger limits
                        if trigger_limits.size != 0:

                            # Convert to UTCDateTime and find the trigger-on time closest to the coincidence trigger
                            trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
                            pick_time = trigger_on[np.argmin(abs(trigger_on - coin_trigger_time))]

                            # Craft waveform stream ID
                            tr_id = tr.id.split('.')
                            waveform_id = WaveformStreamID(tr_id[0], tr_id[1], tr_id[2], tr_id[3])

                            # Create ObsPy pick and ObsPy arrival objects and add to redpy_catalog
                            add_pick = Pick(time=pick_time, waveform_id=waveform_id, phase_hint='P')
                            redpy_catalog[contribution_index].picks.append(add_pick)

    # Write the fully picked redpy catalog to an xml file
    writer(convert_redpy_output_dir + 'redpy_catalog_picked.xml', redpy_catalog)

    # Conclude process
    time_stop = time.time()
    print('Pick making complete, processing time: %.2f s' % (time_stop-time_start))

    #%% Pull the fully picked core catalog from picked redpy catalog

    # Call pull_cores function
    core_catalog = pull_cores(redpy_catalog)

    # Write picked core catalog to .xml file
    writer(convert_redpy_output_dir + 'core_catalog_picked.xml', core_catalog)

def create_tribe(convert_redpy_output_dir,
                 create_tribe_output_dir,
                 data_path,
                 template_stations,
                 resampling_frequency,
                 lowcut,
                 highcut,
                 prepick,
                 length,
                 min_snr,
                 process_len=86400,
                 tolerance=4e4,
                 channel_convention=True,
                 use_all_analyst_events=False):

    """
    Creates a tribe of templates for EQcorrscan
    :param convert_redpy_output_dir (str): output directory for the convert_redpy step
    :param create_tribe_output_dir (str): output directory for the create_tribe step
    :param data_path (str): path to where relevant miniseed files can be found
    :param template_stations (str or list): SEED station code(s) accepted for template creation
    :param resampling_frequency (float): desired sampling rate for templates
    :param lowcut (float): bandpass filter lower bound (Hz)
    :param highcut (float): bandpass filter upper bound (Hz)
    :param prepick (float): time before pick time to start template waveform trim (s)
    :param length (float): time from pre-pick to stop template waveform trim (s)
    :param min_snr (float): minimum signal-to-noise to accept waveform into template
    :param process_len (float): maximum data length (in seconds) to process at once (defaults to 86400 s)
    :param tolerance (float): factor to median tolerance for boxcar removal from data (as a factor to median)
    :param channel_convention (bool): if `True`, enforce strict compliance for P/S picks to be on vertical/horizontal components
    :param use_all_analyst_events (bool): if `True`, disregard REDPy clustering and create templates from all analyst-derived events
    :return: tribe (:class:`~core.match_filter.tribe`): tribe of successfully crafted templates
    """

    # Import all dependencies
    import time
    import numpy as np
    from itertools import compress
    from obspy import UTCDateTime, Catalog
    from eqcorrscan.core.match_filter.tribe import Tribe
    from toolbox import prepare_catalog_stream, reader, writer

    # Commence tribe creation
    print('Commencing tribe creation...')

    # If template_stations is a comma separated string, convert to list
    if type(template_stations) == list:
        template_stations = template_stations.split(',')

    # Read in, and combine, the picked core catalog and unmatched analyst catalog
    core_catalog_picked = reader(convert_redpy_output_dir + 'core_catalog_picked.xml')
    if use_all_analyst_events:
        unmatched_analyst_events = reader(convert_redpy_output_dir + 'unmatched_analyst_events_core.xml')
    else:
        unmatched_analyst_events = reader(convert_redpy_output_dir + 'unmatched_analyst_events_redpy.xml')
    template_catalog = core_catalog_picked + unmatched_analyst_events

    # Clean catalog to only include picks from our station list
    print('\nRemoving picks on stations not in input station list...')
    for i, event in enumerate(template_catalog):
        for pick in reversed(event.picks):
            if pick.waveform_id.station_code not in template_stations:
                print('Removed ' + pick.waveform_id.station_code + ' pick from template catalog index ' + str(i))
                event.picks.remove(pick)
    print('Done')

    # If channel convention needs to be enforced, remove off-direction picks
    if channel_convention:
        print('\nRemoving picks that do not conform to channel convention...')
        for i, event in enumerate(template_catalog):
            for pick in reversed(event.picks):
                if pick.phase_hint == 'P' and pick.waveform_id.channel_code[-1] != 'Z':
                    print('Removed ' + pick.waveform_id.station_code + 'P pick from template catalog index ' + str(i))
                    event.picks.remove(pick)
                elif pick.phase_hint == 'S' and pick.waveform_id.channel_code[-1] == 'Z':
                    print('Removed ' + pick.waveform_id.station_code + 'S pick from template catalog index ' + str(i))
                    event.picks.remove(pick)
        print('Done')

    # Initialize tribe and tracker for valid events
    tribe = Tribe()
    time_start = time.time()

    # Create a boolean list to keep track of which events have been attempted
    tracker = np.array([False for i in range(len(template_catalog))])

    # Initialize array of all event times
    catalog_times = np.array([event.origins[0].time for event in template_catalog])

    # Loop through events in catalog
    print('\nCommencing tribe creation...')
    for k, event in enumerate(template_catalog):

        # Check if event has had its template creation attempted. If True, skip to next
        if tracker[k]:
            continue
        # Otherwise check for other events occuring on the same day so that data can be loaded at one go
        else:
            event_daystart = UTCDateTime(event.origins[0].time.date)
            event_dayend = event_daystart + 86400
            sub_catalog_index = (np.array(catalog_times) > event_daystart) & (np.array(catalog_times) < event_dayend)
            sub_catalog = Catalog(list(compress(template_catalog, sub_catalog_index))).copy()
            tracker[np.where(sub_catalog_index==True)] = True

        # Prepare catalog stream, trim to the start and end of the day to enforce daylong processing
        stream = prepare_catalog_stream(data_path, sub_catalog, resampling_frequency, tolerance)
        stream = stream.trim(starttime=event_daystart, endtime=event_dayend, pad=True)

        # Do stream checks
        for trace in stream:
            # Check if trace is masked with too many zeros (more than half of samples are masked)
            if hasattr(trace.data, 'mask') and (np.sum(trace.data.mask) / len(trace.data.mask)) > 0.5:
                print('%s.%s got removed due to overly large data gaps.' % (trace.stats.station, trace.stats.channel))
                stream.remove(trace)

        # Construct mini tribe out of sub catalog
        try:
            sub_tribe = Tribe().construct(
                method="from_meta_file", lowcut=lowcut, highcut=highcut, samp_rate=resampling_frequency, length=length,
                filt_order=4, prepick=prepick, meta_file=sub_catalog, st=stream, process=True,
                process_len=process_len, min_snr=min_snr, parallel=True)
            if sub_tribe is not None or len(sub_tribe) > 0:
                tribe += sub_tribe

        # If tribe fails then we add the current sub catalog to the client catalog to try the client method later
        except:
            print('Tribe creation attempt on %s to %s failed.' % (event_daystart, event_dayend))

    # Conclude process
    time_end = time.time()
    print('\nTemplate creation attempts complete. Time elapsed: %.2f s' % (time_end-time_start))

    print('%d out of %d events in the catalog were converted to templates.' % (len(tribe),len(template_catalog)))
    writer(create_tribe_output_dir + 'tribe.tgz', tribe)

    return tribe

def scan_data(tribe,
              scan_data_output_dir,
              data_path,
              min_stations,
              min_picks,
              starttime,
              endtime,
              resampling_frequency,
              threshold_type,
              threshold,
              trig_int,
              decluster=True,
              decluster_metric='avg_cor',
              max_zeros=100,
              npts_threshold=100,
              tolerance=4e4,
              parallel_process=True):

    """
    Conduct a matched-filter scan on a user-specified duration of seismic data using the input Tribe of templates.
    :param tribe (:class:`~core.match_filter.tribe`): tribe of templates for matched-filter scan
    :param convert_redpy_output_dir (str): output directory for the convert_redpy step
    :param scan_data_output_dir (str): output directory for the scan_data step
    :param data_path (str): path to where relevant miniseed files can be found
    :param min_stations (int): minimum number of stations where each template must be observed on before being used for matched-filter
    :param min_picks (int): minimum number of picks that each template must possess before being used for matched-filter
    :param starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): start time for EQcorrscan matched-filter scan
    :param endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): end time for EQcorrscan matched-filter scan
    :param resampling_frequency (float): standardized sampling rate used for seismic data matched-filter scan (make sure this is equal to template resampling_frequency)
    :param threshold_type (str): EQcorrscan threshold type -- choose between 'MAD', 'absolute', 'av_chan_corr'
    :param threshold (float): threshold value used for matched-filter detections
    :param trig_int (float): minimum trigger interval for individual template (s)
    :param decluster (bool): if `True`, decluster matched-filter output such that different templates are not allowed to trigger within trig_int
    :param decluster_metric (str): metric to sort peaks by. 'avg_cor' takes the single station average correlation, 'cor_sum' takes the total correlation sum across all channels, 'thresh_exc' takes the factor by how much the detection exceeded the input threshold
    :param max_zeros (int): maximum number of zeros allowed in trace segment prior to merging (defaults to 100)
    :param npts_threshold (int): minimum number of samples needed for each trace segment prior to merging (defaults to 100)
    :param tolerance (float): factor to median tolerance for boxcar removal from data (as a factor to median)
    :param parallel_process (bool): if `True`, use parallel processing for matched-filter scan
    :return: party (:class:`~core.match_filter.party`): Party of all matched-filter detections grouped by parent template
    :return: detected_catalog (:class:`~obspy.core.event.Catalog`): Catalog of all matched-filter detections a.k.a. the temporally enhanced event list
    :return: relocatable_catalog (:class:`~obspy.core.event.Catalog`): Catalog of all matched-filter detections with adopted locations
    """

    # Import all dependencies
    import time
    import glob
    import numpy as np
    from obspy import read, Stream, Catalog, UTCDateTime
    from obspy.core.event import Origin
    from eqcorrscan import Party
    from toolbox import remove_boxcars, remove_bad_traces, writer

    # Clean tribe off templates using min_stations and min_picks
    print('\nBefore cleaning:', tribe)

    # First remove template traces that are all nans
    for t in tribe:
        for tr in t.st:
            if np.isnan(tr.data).all():
                t.st.traces.remove(tr)

    # Remove based on min_stations
    tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= min_stations]
    print('After removing templates with < %d stations:' % min_stations, tribe)

    # Remove based on min_picks
    tribe.templates = [t for t in tribe if len(t.event.picks) >= min_picks]
    print('After removing templates with < %d valid picks:' % min_picks, tribe)

    # Get a unique list of all template station-channel combinations for data fetching
    data_fileheads = []

    # Loop through trace information in every template's stream
    for template in tribe:
        for trace in template.st:
            # Extract station and channel names from trace
            sta = trace.stats.station
            chan = trace.stats.channel

            # Craft filehead and append
            data_filehead = data_path + sta + '.' + chan + '.'
            data_fileheads.append(data_filehead)

    # Now compile unique and sort
    data_fileheads = list(set(data_fileheads))
    data_fileheads.sort()

    # Prepare to loop over days for scan
    party_all = Party()
    num_days = int(np.floor((endtime - starttime) / 86400))
    time_start = time.time()

    # Commence loop over days
    for i in range(num_days):

        # Define temporal boundaries of our day's scan
        t1 = starttime + (i * 86400)
        t2 = starttime + ((i + 1) * 86400)
        print('\nNow at %s...' % str(t1))

        # Initialize stream object and fetch data from data_path
        stream = Stream()

        # Loop over data fileheads, using glob.glob to match data files we want in data_path
        for data_filehead in data_fileheads:
            data_filename = data_filehead + str(t1.year) + ':' + f'{t1.julday :03}' + ':*'
            matching_filenames = (glob.glob(data_filename))

            # Try to read the miniseed data file and add it to the stream (occasionally fails)
            for matching_filename in matching_filenames:
                try:
                    stream_contribution = read(matching_filename)
                    stream = stream + stream_contribution
                except:
                    continue

        # Remove bad traces (too many zeros or too little samples)
        stream = remove_bad_traces(stream, max_zeros=max_zeros, npts_threshold=npts_threshold)

        # Process stream (remove spikes, downsample to match tribe, detrend, merge)
        stream = remove_boxcars(stream, tolerance)
        stream = stream.resample(sampling_rate=resampling_frequency)
        stream = stream.detrend("simple")
        stream = stream.merge()
        stream = stream.trim(starttime=t1, endtime=t2, pad=True)

        # Check if stream is already empty. If yes, declare failure and skip
        if stream is None:
            print('The stream is already empty.')
            time_current = time.time()
            print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))
            continue

        # Remove traces that are overly masked (i.e. more zeros than actual data points)
        for tr in stream:
            if len(np.nonzero(tr.data)[0]) < 0.5 * len(tr.data):
                stream.remove(tr)

        print('Stream despiked, resampled, merged and trimmed. Getting party of detections...')

        # Attempt to scan the current day
        try:
            party = tribe.detect(
                stream=stream, threshold=threshold, threshold_type=threshold_type, trig_int=trig_int, daylong=True,
                overlap=None, parallel_process=parallel_process)

            # Append party to party_all
            party_all = party_all + party
            time_current = time.time()
            print('Party created, appending. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))

        # If the scan for the day fails, print a notification
        except:
            time_current = time.time()
            print('Party failed, skipping. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))
            continue

    # Conclude process
    time_end = time.time()
    print('\nParty creation complete. Time taken: %.2f hours' % ((time_end - time_start) / 3600))

    # Remove detections that are anchored by too little stations
    print('\nFiltering away detections anchored by too little stations...')
    for family in party_all:
        for detection_index in reversed(range(len(family.detections))):
            detection = family[detection_index]
            if detection.no_chans < min_stations:
                family.detections.pop(detection_index)

    # Clean the party off of repeats (different templates that detect the "same" event)
    if decluster:
        print('Declustering party with leftover detections...')
        party_all = party_all.decluster(trig_int=trig_int, metric=decluster_metric)

    # Convert party into an ObsPy catalog
    detected_catalog = party_all.get_catalog()
    print('The matched-filter scan detected a total of %d events.' % len(detected_catalog))

    # Extract relocated catalog from full detected catalog
    relocatable_catalog = Catalog()

    # Extract all template names from tribe
    template_names = [template.name for template in tribe]

    # Loop over the detected catalog, copying over events if their templates have locations
    # Also give each event pick phase information based on channel convention
    for event in detected_catalog:

        # Get source event
        source_template_name = event.comments[0].text.split(' ')[1]
        source_template_index = template_names.index(source_template_name)
        source_template = tribe[source_template_index]
        source_event = source_template.event

        # Check if the source event has a location. Only continue if it does
        if source_event.origins[0].latitude is not None:

            # We extract the source event's location (lat, lon, dep)
            lat = source_event.origins[0].latitude
            lon = source_event.origins[0].longitude
            dep = source_event.origins[0].depth

            # Reformat the event's detection time to UTC and calculate estimated origin time
            time_text = event.resource_id.id.split('_')[-1]
            detection_utc = UTCDateTime(time_text)
            template_origin_dt = np.min([tr.stats.starttime for tr in source_template.st]) - source_event.origins[0].time
            origin_utc = detection_utc - template_origin_dt

            # Create an ObsPy Origin object and store it in the event
            event.origins = [Origin(time=origin_utc, latitude=lat, longitude=lon, depth=dep)]

            # Loop over picks in event
            for pick in event.picks:

                # Check for non-vertical channel, and give an 'S' hint
                if pick.waveform_id.channel_code[-1] == 'N' or pick.waveform_id.channel_code[-1] == 'E':
                    pick.phase_hint = 'S'

                # Also check for a vertical channel, and give an 'P' hint
                elif pick.waveform_id.channel_code[-1] == 'Z':
                    pick.phase_hint = 'P'

                # Throw out a ValueError if the channel code is not recognized
                else:
                    raise ValueError

            # Append to new catalog
            relocatable_catalog += event

    print('The relocatable catalog consists of %d events.' % len(relocatable_catalog))

    # Write out party and catalogs
    writer(scan_data_output_dir + 'party.tgz', party_all)
    writer(scan_data_output_dir + 'detected_catalog.xml', detected_catalog)
    writer(scan_data_output_dir + 'relocatable_catalog.xml', relocatable_catalog)

    return party, detected_catalog, relocatable_catalog

def rethreshold_results(tribe,
                        party,
                        threshold_type,
                        new_threshold,
                        decluster,
                        trig_int):

    """
    Rethreshold party of detections using inbuilt function and re-extract detected and relocated catalogs
    :param tribe (:class:`~core.match_filter.tribe`): tribe of templates for matched-filter scan:
    :param party (:class:`~core.match_filter.party`): Party of all matched-filter detections grouped by parent template
    :param threshold_type (str): EQcorrscan threshold type -- choose between 'MAD', 'absolute', 'av_chan_corr'
    :param new_threshold (float): new threshold used to discard unwanted detections
    :param decluster (bool): if `True`, decluster matched-filter output such that different templates are not allowed to trigger within trig_int
    :param trig_int (float): minimum trigger interval for individual template (s)
    :return: party (:class:`~core.match_filter.party`): Party of all matched-filter detections after rethresholding, grouped by parent template
    :return: detected_catalog (:class:`~obspy.core.event.Catalog`): Catalog of all matched-filter detections a.k.a. the temporally enhanced event list
    :return: relocatable_catalog (:class:`~obspy.core.event.Catalog`): Catalog of all matched-filter detections with adopted locations
    """

    # Import all dependencies
    from obspy import Catalog, UTCDateTime
    from obspy.core.event import Origin

    # Use EQcorrscan's rethreshold function
    print('\nRethresholding outputs using %s = %.2f...' % (threshold_type, new_threshold))
    party = party.rethreshold(new_threshold=new_threshold,
                              new_threshold_type=threshold_type)

    # Clean the party off of repeats (different templates that detect the "same" event)
    if decluster:
        print('Declustering party with leftover detections...')
        party = party.decluster(trig_int=trig_int)

    # Convert party into an ObsPy catalog
    detected_catalog = party.get_catalog()
    print('The rethresholded catalog has a total of %d events.' % len(detected_catalog))

    # Extract relocated catalog from full detected catalog
    relocatable_catalog = Catalog()

    # Extract all template names from tribe
    template_names = [template.name for template in tribe]

    # Loop over the detected catalog, copying over events if their templates have locations
    # Also give each event pick phase information based on channel convention
    for event in detected_catalog:

        # Get source event
        source_template_name = event.comments[0].text.split(' ')[1]
        source_template_index = template_names.index(source_template_name)
        source_template = tribe[source_template_index]
        source_event = source_template.event

        # Check if the source event has a location. Only continue if it does
        if source_event.origins[0].latitude is not None:

            # We extract the source event's location (lat, lon, dep)
            lat = source_event.origins[0].latitude
            lon = source_event.origins[0].longitude
            dep = source_event.origins[0].depth

            # Reformat the event's detection time to UTC and calculate estimated origin time
            time_text = event.resource_id.id.split('_')[-1]
            detection_utc = UTCDateTime(time_text)
            template_origin_dt = np.min([tr.stats.starttime for tr in source_template.st]) - source_event.origins[0].time
            origin_utc = detection_utc - template_origin_dt

            # Create an ObsPy Origin object and store it in the event
            event.origins = [Origin(time=origin_utc, latitude=lat, longitude=lon, depth=dep)]

            # Loop over picks in event
            for pick in event.picks:

                # Check for non-vertical channel, and give an 'S' hint
                if pick.waveform_id.channel_code[-1] == 'N' or pick.waveform_id.channel_code[-1] == 'E':
                    pick.phase_hint = 'S'

                # Also check for a vertical channel, and give an 'P' hint
                elif pick.waveform_id.channel_code[-1] == 'Z':
                    pick.phase_hint = 'P'

                # Throw out a ValueError if the channel code is not recognized
                else:
                    raise ValueError

            # Append to new catalog
            relocatable_catalog += event

    print('The relocatable catalog consists of %d events.' % len(relocatable_catalog))

    return party, detected_catalog, relocatable_catalog

def generate_dtcc(catalog,
                  relocate_catalog_output_dir,
                  data_path,
                  pre_pick_actual,
                  pre_pick_excess,
                  length_actual,
                  length_excess,
                  shift_len,
                  resampling_frequency,
                  lowcut,
                  highcut,
                  min_cc,
                  min_link,
                  max_sep=100,
                  weight_func=None,
                  make_s_picks=False,
                  pre_pick_actual_S=None,
                  length_actual_S=None,
                  shift_len_S=None,
                  parallel_process=True):

    """
    Generate HypoDD/GrowClust cross-correlation differential times file (dt.cc)
    :param catalog (:class:`~obspy.core.event.Catalog`): catalog of relocatable events
    :param relocate_catalog_output_dir (str): output directory for the relocate_catalog step
    :param data_path (str): path to where relevant miniseed files can be found
    :param pre_pick_actual (float): length of pre-pick to include in phase segment (s)
    :param pre_pick_excess (float): length in excess of pre-pick for waveform retrieval (s)
    :param length_actual (float): length of phase segment from pre-pick (s)
    :param length_excess (float): length in excess of phase segment from pre-pick for waveform retrieval (s)
    :param shift_len (float): length of time shift to find the maximum cross-correlation coefficient between phase segments (s)
    :param resampling_frequency (float): desired sampling rate for stream dictionary (Hz) (upsample to 100Hz for 0.01s precision)
    :param lowcut (float): bandpass filter lower bound (Hz)
    :param highcut (float): bandpass filter upper bound (Hz)
    :param min_cc (float): minimum cross-correlation coefficient for a differential time to be retained for an event pair
    :param min_link (int): minimum number of differential time observations for an event pair to be included in the output
    :param max_sep (float): maximum separation between event pairs allowed (km) (defaults to 100 to include all event pairs)
    :param weight_func (function): function to be applied on raw cross-correlation coefficient squared value computed by write_correlations()
    :param make_s_picks (bool): if `True`, make theoretical S picks from P picks using an vp/vs ratio of 1.73 before constructing dt.cc
    :param pre_pick_actual_S (float): if defined, use a different pre_pick_actual for S picks. pre_pick_actual argument will be used for P picks only.
    :param length_actual_S (float): if defined, use a different length_actual for S picks. length_actual argument will be used for P picks only.
    :param shift_len_S (float): if defined, use a different shift_len for S picks. shift_len argument will be used for P picks only.
    :param parallel_process (bool): if `True`, use parallel processing to calculate cross-correlation differential times
    :return:
    """

    # Import all dependencies
    import os
    import numpy as np
    from eqcorrscan.utils.catalog_to_dd import write_correlations, _generate_event_id_mapper
    from obspy import Catalog
    from obspy.core.event import Event
    from toolbox import prepare_stream_dict, adjust_weights, estimate_s_picks

    # Make estimated S pick if user specifies
    if make_s_picks:
        catalog = estimate_s_picks(catalog)

    # Define S-pick specific params if they were provided, else copy parent param
    pre_pick_actual_S = pre_pick_actual_S or pre_pick_actual
    length_actual_S = length_actual_S or length_actual
    shift_len_S = shift_len_S or shift_len

    # Separate the catalog into two sub-catalogs (one with Ps, one with Ss)
    catalog_P = Catalog()
    catalog_S = Catalog()
    for e in catalog:
        ev1 = Event(origins=e.origins, comments=e.comments, magnitudes=e.magnitudes)
        ev2 = Event(origins=e.origins, comments=e.comments, magnitudes=e.magnitudes)
        for pick in e.picks:
            if pick.phase_hint == 'P':
                ev1.picks.append(pick)
            else:
                ev2.picks.append(pick)
        catalog_P.events.append(ev1)
        catalog_S.events.append(ev2)

    # Process the entire catalog at one go
    print('\nPreparing stream dictionaries for differential times computation...')

    # Generate stream dictionary (refer to toolbox.py)
    # Note that we want pre_pick and length to be in excess, since write_correlations trims the data for us
    stream_dict_P = prepare_stream_dict(catalog_P, pre_pick=pre_pick_excess, length=length_excess,
                                        resampling_frequency=resampling_frequency, local=True, data_path=data_path)
    stream_dict_S = prepare_stream_dict(catalog_S, pre_pick=pre_pick_excess, length=length_excess,
                                        resampling_frequency=resampling_frequency, local=True, data_path=data_path)

    # Execute cross correlations and write out a .cc file using write_correlations (refer to EQcorrscan docs)
    # Note this stores a file called "dt.cc" in your current working directory
    print('Correlating P-arrivals...')
    _ = write_correlations(catalog=catalog_P, stream_dict=stream_dict_P, extract_len=length_actual,
                           pre_pick=pre_pick_actual, shift_len=shift_len, lowcut=lowcut,
                           highcut=highcut, max_sep=max_sep, min_link=0, min_cc=min_cc,
                           interpolate=False, max_workers=None, parallel_process=parallel_process)

    # Define source and target cc filepaths, then copy over to relocate_catalog output directory
    original_dt_path = './dt.cc'
    target_dt_path_P = relocate_catalog_output_dir + 'dt.ccP'
    adjust_weights(original_dt_path, target_dt_path_P, dt_cap=shift_len, min_link=1,
                   append=False, weight_func=weight_func)
    os.remove(original_dt_path)

    # Correlate S-arrivals
    print('Correlating S-arrivals...')
    _ = write_correlations(catalog=catalog_S, stream_dict=stream_dict_S, extract_len=length_actual_S,
                           pre_pick=pre_pick_actual_S, shift_len=shift_len_S, lowcut=lowcut,
                           highcut=highcut, max_sep=max_sep, min_link=0, min_cc=min_cc,
                           interpolate=False, max_workers=None, parallel_process=parallel_process)
    target_dt_path_S = relocate_catalog_output_dir + 'dt.ccS'
    adjust_weights(original_dt_path, target_dt_path_S, dt_cap=shift_len_S, min_link=1,
                   append=False, weight_func=weight_func)
    os.remove(original_dt_path)

    # Combine dt.cc files for P and S
    dtccP_textfile = open(target_dt_path_P, 'r')
    dtccP_lines = dtccP_textfile.readlines()
    dtccS_textfile = open(target_dt_path_S, 'r')
    dtccS_lines = dtccS_textfile.readlines()
    if dtccP_lines != [] and dtccS_lines == []:
        dtcc_lines = dtccP_lines
    elif dtccP_lines == [] and dtccS_lines != []:
        dtcc_lines = dtccS_lines
    else:
        dtccP_dict = {l[0:25]:l[25:] for l in ''.join(dtccP_lines).split('#')}
        dtccS_dict = {l[0:25]:l[25:] for l in ''.join(dtccS_lines).split('#')}
        dtcc_dict = {key: dtccP_dict.get(key,'') + dtccS_dict.get(key,'')
                    for key in set(list(dtccP_dict.keys())+list(dtccS_dict.keys()))}
        dtcc_lines = '#'.join([key+val for key,val in dtcc_dict.items()]).split('\n')
        dtcc_lines = [l+'\n' for l in (dtcc_lines)]

    # Remove event pairs that have less than min_link observations
    header_indices = np.flatnonzero([l[0] == '#' for l in dtcc_lines])
    num_link = np.diff(np.append(header_indices, len(dtcc_lines)-1)) - 1
    failed_pair_indices = header_indices[np.flatnonzero(num_link < min_link)]
    for failed_pair_index in failed_pair_indices:
        dtcc_lines[failed_pair_index] = ''
        if failed_pair_index != (len(dtcc_lines) - 1):
            i = 1
            while dtcc_lines[failed_pair_index+i][0] != '#' and (failed_pair_index+i) != (len(dtcc_lines)-1):
                dtcc_lines[failed_pair_index+i] = ''
                i += 1
    dtcc_lines = [l for l in dtcc_lines if (l != '' and l != '\n')]

    with open(relocate_catalog_output_dir + 'dt.cc', "w") as open_file:
        open_file.write(''.join(dtcc_lines))
    open_file.close()

    print('Cross correlations done!')


def run_hypoDD(catalog,
               relocate_catalog_output_dir,
               stations,
               stalats,
               stalons,
               staelevs,
               vzmodel_path,
               has_ps_ratio,
               correct_depths,
               ph2dt_inc_dict,
               ph2dt_inp_dict,
               hypoDD_inc_dict,
               hypoDD_inp_dict,
               hypoDD_dir='./hypoDD/'):
    """
    Run hypoDD using input .inc and .inp files
    :param catalog (:class:`~obspy.core.event.Catalog`): catalog of all relocatable events
    :param relocate_catalog_output_dir (str): output directory for the relocate_catalog step
    :param stations (str or list): SEED code(s) for stations allowed for differential time calculation
    :param stalats (str or list or typle): list (or tuple or comma separated string) of station latitudes [-90,90]
    :param stalons (str or list or typle): list (or tuple or comma separated string) of station longitudes [-180,180]
    :param staelevs (str or list or typle): list (or tuple or comma separated string) of station elevations (m)
    :param vzmodel_path (str): path to velocity model text file
    :param has_ps_ratio (bool): set to `True` if velocity model text file has a column of Vp/Vs ratios. If `False`, assume 1.75 for all layers.
    :param correct_depths (bool): set to `True` to shift the entire relocation problem downwards using the velocity model ceiling to avoid erroneous airquake elimination
    :param ph2dt_inc_dict (dict): dictionary of all ph2dt.inc arguments
    :param ph2dt_inp_dict (dict): dictionary of all ph2dt.inp arguments
    :param hypoDD_inc_dict (dict): dictionary of all hypoDD.inc arguments
    :param hypoDD_inp_dict (dict): dictionary of all hypoDD.inp arguments
    :param hypoDD_dir (str): path to hypoDD directory
    :return: hypoDD_loc (:class:`~obspy.core.event.Catalog`): catalog of all relocation candidate events
    :return: hypoDD_reloc (:class:`~obspy.core.event.Catalog`): catalog of all successfully relocated events
    """

    # Import all dependencies
    import os
    import tarfile
    import subprocess
    import pandas as pd
    from obspy import UTCDateTime
    from obspy.core.event import Event
    from toolbox import writer, loc2cat
    from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper

    # Initialize hypoDD sub directories
    print('\nCreating sub directories...')
    for dir in [hypoDD_dir + 'IN', hypoDD_dir + 'INP', hypoDD_dir + 'OUT_ph2dt', hypoDD_dir + 'OUT_hypoDD']:
        try:
            os.mkdir(dir)
        except FileExistsError:
            print('Subdirectory %s already exists.' % dir)

    # If depth_correction = True, we want to shift the entire problem downwards
    if correct_depths:
        print('\nWARNING: correcting elevations and depths for input catalog, vzmodel and stations')
        # Read the input vzmodel and use the
        vzmodel_read = open(vzmodel_path, "r")
        vzmodel_lines = vzmodel_read.readlines()
        model_ceiling = float(vzmodel_lines[0].split()[0])
        depth_correction = -1 * model_ceiling * 1000  # m
        vzmodel_read.close()
    else:
        depth_correction = 0

    # Deconstruct .xml file
    print('\nDeconstructing catalog into phase.dat and stations.dat ...')

    # Prepare output filepaths
    STATION_DAT_FILEPATH = hypoDD_dir + 'IN/stations.dat'
    PHASE_DAT_FILEPATH = hypoDD_dir + 'IN/phase.dat'
    DETECTION_DAT_FILEPATH = hypoDD_dir + 'IN/detection.dat'

    # Prepare file formats
    STATION_DAT_FORMAT = '%s %.4f %.4f %.1f\n'
    PHASE_DAT_FORMAT_1 = '#  %4d %2d %2d %2d %2d %2d.%02d %8.4f %9.4f %8.3f %5.2f %.2f %.2f %.2f %10d\n'
    PHASE_DAT_FORMAT_2 = '%s %.3f %.2f  %2s\n'
    DETECTION_DAT_FORMAT = '%4d%02d%02d  %2d%02d%02d%02d %9.4f %10.4f %10.3f %6.2f %7.2f %7.2f %6.2f %10d\n'

    # Write station.dat file
    if type(stations) == type(stalats) == type(stalons) == type(staelevs) == str:
        stations = stations.split(',')
        stalats = [float(n) for n in stalats.split(',')]
        stalons = [float(n) for n in stalons.split(',')]
        staelevs = [float(n) for n in staelevs.split(',')]
    elif type(stations) == type(stalats) == type(stalons) == type(staelevs) == tuple:
        stations = list(stations)
        stalats = [float(n) for n in stalats]
        stalons = [float(n) for n in stalons]
        staelevs = [float(n) for n in staelevs]
    if correct_depths:
        staelevs = [n - depth_correction for n in staelevs]
    STATION_DAT_FILE = open(STATION_DAT_FILEPATH, 'w')
    for i in range(len(stations)):
        STATION_DAT_FILE.write(STATION_DAT_FORMAT % (stations[i],
                                                     stalats[i],
                                                     stalons[i],
                                                     staelevs[i]))
    STATION_DAT_FILE.close()

    # Generate event id mapper
    event_id_mapper = _generate_event_id_mapper(catalog, event_id_mapper=None)
    # Initialize phase.dat and detection.dat and write both file simultaneously
    PHASE_DAT_FILE = open(PHASE_DAT_FILEPATH, 'w')
    DETECTION_DAT_FILE = open(DETECTION_DAT_FILEPATH, 'w')
    # Loop through events in catalog
    for event in catalog:
        # Extract key information
        yr = int(event.origins[0].time.year)
        mo = int(event.origins[0].time.month)
        dy = int(event.origins[0].time.day)
        hh = int(event.origins[0].time.hour)
        mm = int(event.origins[0].time.minute)
        ss = int(event.origins[0].time.second)
        cs = int(event.origins[0].time.microsecond / 1e4)
        lat = event.origins[0].latitude
        lon = event.origins[0].longitude
        if correct_depths:
            dep = (event.origins[0].depth + depth_correction) / 1000
        else:
            dep = (event.origins[0].depth) / 1000
        if event.magnitudes == []:
            print('Event index %d, UTC Time %s has no magnitude, skipping.' % (
            catalog.events.index(event), event.origins[0].time))
            continue
        else:
            mag_object = event.preferred_magnitude() or event.magnitudes[0]
            mag = mag_object.mag
        eh = 0  # Change accordingly if event object contains error info
        ez = 0  # Change accordingly if event object contains error info
        rms = 0  # Change accordingly if event object contains error info
        evid = event_id_mapper[event.resource_id.id]
        # Check event time against template to see if it is a self detection
        event_time = event.origins[0].time
        try:
            template_name = event.comments[0].text.split()[1]
            template_time_segments = [int(numstr) for numstr in template_name.replace('t', '_').split('_')]
            template_time = UTCDateTime(*template_time_segments)
            time_diff = abs(event_time - template_time)
        except:
            time_diff = 0
        # If it is a self detection, write info into phase file for ph2dt
        if time_diff < 1:
            PHASE_DAT_FILE.write(
                PHASE_DAT_FORMAT_1 % (yr, mo, dy, hh, mm, ss, cs, lat, lon, dep, mag, eh, ez, rms, evid))
            for i, pick in enumerate(event.picks):
                sta = pick.waveform_id.station_code
                tt = pick.time - event_time
                if event.origins[0].arrivals == []:
                    wght = 0.1  # Implement weighting scheme
                else:
                    wght = event.origins[0].arrivals[i].time_weight
                pha = pick.phase_hint
                PHASE_DAT_FILE.write(PHASE_DAT_FORMAT_2 % (sta, tt, wght, pha))
        else:
            # Otherwise, it is a detection, and we write it into our detection list file for appending later
            DETECTION_DAT_FILE.write(
                DETECTION_DAT_FORMAT % (yr, mo, dy, hh, mm, ss, cs, lat, lon, dep, mag, eh, ez, rms, evid))
    # Close files
    PHASE_DAT_FILE.close()
    DETECTION_DAT_FILE.close()
    print('Catalog parsed into phase.dat and detection.dat successfully.')

    # Extract contents from hypoDD tarfile
    print('\nExtracting hypoDD tarfile contents ...')
    hypoDD_tar = tarfile.open(hypoDD_dir + 'HYPODD_2.1b.tar.gz', "r:gz")
    hypoDD_tar.extractall(hypoDD_dir)
    print('hypoDD tarfile extracted.')

    # Configure ph2dt.inc
    ph2dt_inc_text = ['c ph2dt.inc: Stores parameters that define array dimensions in ph2dt.',
                      'c            Modify to fit size of problem and available computer memory.',
                      'c Parameter Description:',
                      'c MEV:   Max number of events.',
                      'c MSTA:  Max number of stations.',
                      'c MOBS:  Max number of phases (P&S) per eventer event.',
                      '',
                      '      integer	MEV, MSTA, MOBS',
                      '',
                      '      parameter(MEV=    %d,' % ph2dt_inc_dict['MEV'],
                      '     &		MSTA=   %d,' % ph2dt_inc_dict['MSTA'],
                      '     &		MOBS=   %d)' % ph2dt_inc_dict['MOBS']]
    with open(hypoDD_dir + '/HYPODD/include/ph2dt.inc', "w") as open_file:
        open_file.write('\n'.join(ph2dt_inc_text))
    open_file.close()

    # Configure hypoDD.inc
    hypoDD_inc_text = ['c hypoDD.inc: Stores parameters that define array dimensions in hypoDD.',
                       'c             Modify to fit size of problem and available computer memory.',
                       'c	      If 3D raytracing is used, also set model parameters in vel3d.inc.',
                       'c Parameter Description:',
                       'c MAXEVE:   Max number of events (must be at least the size of the number ',
                       'c           of events listed in the event file)',
                       'c MAXDATA:  Max number of observations (must be at least the size of the ',
                       'c           number of observations).  ',
                       'c MAXEVE0:  Max number of events used for SVD. If only LSQR is used, ',
                       'c           MAXEVE0 can be set to 2 to free up memory. ',
                       'c MAXDATA0: Max number of observations used for SVD. If only LSQR is used, ',
                       'c           MAXDATA0 can be set to 1 to free up memory. ',
                       'c MAXLAY:   Max number of model layers.',
                       'c MAXSTA:   Max number of stations.',
                       'c MAXCL:    Max number of clusters allowed. ',
                       '      integer*4 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0',
                       '      integer*4 MAXCL',
                       '',
                       '      parameter(MAXEVE   = %d,' % hypoDD_inc_dict['MAXEVE'],
                       '     &          MAXDATA  = %d,' % hypoDD_inc_dict['MAXDATA'],
                       '     &          MAXEVE0  = %d,' % hypoDD_inc_dict['MAXEVE0'],
                       '     &          MAXDATA0 = %d,' % hypoDD_inc_dict['MAXDATA0'],
                       '     &          MAXLAY   = %d,' % hypoDD_inc_dict['MAXLAY'],
                       '     &          MAXSTA   = %d,' % hypoDD_inc_dict['MAXSTA'],
                       '     &          MAXCL    = %d)' % hypoDD_inc_dict['MAXCL']]
    with open(hypoDD_dir + '/HYPODD/include/hypoDD.inc', "w") as open_file:
        open_file.write('\n'.join(hypoDD_inc_text))
    open_file.close()

    # Clean and recompile within src dir
    print('Compiling ...')
    clean_command = 'make clean -C %s' % (hypoDD_dir + 'HYPODD/src/')
    make_command = 'make -C %s' % (hypoDD_dir + 'HYPODD/src/')
    subprocess.call(clean_command, shell=True)
    subprocess.call(make_command, shell=True)  # not working for PyCharm IDE?

    # Configure ph2dt.inp
    ph2dt_inp_text = ['* ph2dt.inp - input control file for program ph2dt',
                      '* Input station file:',
                      '%s' % STATION_DAT_FILEPATH,
                      '* Input phase file:',
                      '%s' % PHASE_DAT_FILEPATH,
                      '*MINWGHT: min. pick weight allowed [0]',
                      '*MAXDIST: max. distance in km between event pair and stations [200]',
                      '*MAXSEP: max. hypocentral separation in km [10]',
                      '*MAXNGH: max. number of neighbors per event [10]',
                      '*MINLNK: min. number of links required to define a neighbor [8]',
                      '*MINOBS: min. number of links per pair saved [8]',
                      '*MAXOBS: max. number of links per pair saved [20]',
                      '*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS',
                      '   %d      %d     %d     %d     %d      %d     %d' %
                      (ph2dt_inp_dict['MINWGHT'],
                       ph2dt_inp_dict['MAXDIST'],
                       ph2dt_inp_dict['MAXSEP'],
                       ph2dt_inp_dict['MAXNGH'],
                       ph2dt_inp_dict['MINLNK'],
                       ph2dt_inp_dict['MINOBS'],
                       ph2dt_inp_dict['MAXOBS'])]
    with open(hypoDD_dir + '/INP/ph2dt.inp', "w") as open_file:
        open_file.write('\n'.join(ph2dt_inp_text))
    open_file.close()

    # Run ph2dt
    ph2dt_command = hypoDD_dir + 'HYPODD/src/ph2dt/ph2dt' + ' ' + hypoDD_dir + '/INP/ph2dt.inp'
    subprocess.call(ph2dt_command, shell=True)

    # Move ph2dt outputs to "hypoDDver/OUT_ph2dt/" and move log to relocate_catalog_output_dir
    to_move = ['station.sel', 'event.sel', 'event.dat', 'dt.ct']
    for filename in to_move:
        move_command = 'mv ./%s %sOUT_ph2dt/%s' % (filename, hypoDD_dir, filename)
        subprocess.call(move_command, shell=True)
    move_command = 'mv ./ph2dt.log %sph2dt.log' % relocate_catalog_output_dir
    subprocess.call(move_command, shell=True)

    # Stitch event.dat and detection.dat to create event_hypoDD.dat
    with open(hypoDD_dir + 'OUT_ph2dt/event.dat') as open_file:
        EVENT_DAT = open_file.read()
    with open(DETECTION_DAT_FILEPATH) as open_file2:
        DETECTION_DAT = open_file2.read()
    with open(hypoDD_dir + 'OUT_ph2dt/event_hypoDD.dat', 'w') as open_file3:
        open_file3.write(EVENT_DAT + DETECTION_DAT)

    # Prepare velocity model information for hypoDD.inp
    vzmodel = pd.read_csv(vzmodel_path, header=None)
    TOP = []
    VELP = []
    RAT = []
    for i in range(len(vzmodel.values)):
        line = vzmodel.values[i][0]
        if correct_depths:
            TOP.append(str(float(line.split(' ')[0]) + depth_correction/1000))
        else:
            TOP.append(str(line.split(' ')[0]))
        VELP.append(str(line.split(' ')[1]))
        if has_ps_ratio:
            RAT.append(str(line.split(' ')[2]))
        else:
            RAT.append(str(1.75))
    TOP = ' '.join(TOP)
    VELP = ' '.join(VELP)
    RAT = ' '.join(RAT)

    # Configure weighting
    weight_params = [hypoDD_inp_dict['NITER'],
                     hypoDD_inp_dict['WTCCP'],
                     hypoDD_inp_dict['WTCCS'],
                     hypoDD_inp_dict['WRCC'],
                     hypoDD_inp_dict['WDCC'],
                     hypoDD_inp_dict['WTCTP'],
                     hypoDD_inp_dict['WTCTS'],
                     hypoDD_inp_dict['WRCT'],
                     hypoDD_inp_dict['WDCT'],
                     hypoDD_inp_dict['DAMP']]
    # Check if each parameter has the same number of values
    if not all(len(weight_param) == len(weight_params[0]) for weight_param in weight_params[1:]):
        raise ValueError('Not all weight parameters have the same length!')
    # Now craft rows accordingly
    WEIGHT_TABLE = ''
    for i in range(len(weight_params[0])):
        WEIGHT_TABLE = WEIGHT_TABLE + ' %2d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %5d' % tuple(
            [weight_param[i] for weight_param in weight_params])
        if i != len(weight_params[0]) - 1:
            WEIGHT_TABLE = WEIGHT_TABLE + '\n'

    # Configure hypoDD.inp
    hypoDD_inp_text = ['hypoDD_2',
                       '* Parameter control file hypoDD.inp:',
                       '*--- INPUT FILE SELECTION',
                       '* filename of cross-corr diff. time input(blank if not available):',
                       '%s' % (relocate_catalog_output_dir + 'dt.cc'),
                       '* filename of catalog travel time input(blank if not available):',
                       '%s' % (hypoDD_dir + 'OUT_ph2dt/dt.ct'),
                       '* filename of initial hypocenter input:',
                       '%s' % (hypoDD_dir + 'OUT_ph2dt/event_hypoDD.dat'),
                       '* filename of initial hypocenter input:',
                       '%s' % (hypoDD_dir + 'IN/stations.dat'),
                       '*',
                       '*--- OUTPUT FILE SELECTION',
                       '* filename of initial hypocenter output (if blank: output to hypoDD.loc):',
                       '%s' % (hypoDD_dir + 'OUT_hypoDD/hypoDD.loc'),
                       '* filename of relocated hypocenter output (if blank: output to hypoDD.reloc):',
                       '%s' % (hypoDD_dir + 'OUT_hypoDD/hypoDD.reloc'),
                       '* filename of station residual output (if blank: no output written):',
                       '%s' % (hypoDD_dir + 'OUT_hypoDD/hypoDD.sta'),
                       '* filename of data residual output (if blank: no output written):',
                       '%s' % (hypoDD_dir + 'OUT_hypoDD/hypoDD.res'),
                       '* filename of takeoff angle output (if blank: no output written):',
                       '%s' % (hypoDD_dir + 'OUT_hypoDD/hypoDD.src'),
                       '*',
                       '*--- DATA SELECTION:',
                       '* IDAT IPHA DIST',
                       '    %d     %d     %d' % (hypoDD_inp_dict['IDAT'],
                                                 hypoDD_inp_dict['IPHA'],
                                                 hypoDD_inp_dict['DIST']),
                       '*',
                       '*--- EVENT CLUSTERING:',
                       '* OBSCC OBSCT MINDS MAXDS MAXGAP ',
                       '    %d    %d    %d    %d    %d        ' % (hypoDD_inp_dict['OBSCC'],
                                                                   hypoDD_inp_dict['OBSCT'],
                                                                   hypoDD_inp_dict['MINDS'],
                                                                   hypoDD_inp_dict['MAXDS'],
                                                                   hypoDD_inp_dict['MAXGAP']),
                       '*',
                       '*--- SOLUTION CONTROL:',
                       '* ISTART ISOLV IAQ NSET',
                       '    %d    %d    %d    %d ' % (hypoDD_inp_dict['ISTART'],
                                                      hypoDD_inp_dict['ISOLVE'],
                                                      hypoDD_inp_dict['IAQ'],
                                                      hypoDD_inp_dict['NSET']),
                       '*',
                       '*--- DATA WEIGHTING AND REWEIGHTING:',
                       '* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP',
                       '%s' % WEIGHT_TABLE,
                       '*',
                       '*--- FORWARD MODEL SPECIFICATIONS:',
                       '* IMOD',
                       '1',
                       '* TOP: depths of top of layer (km) ',
                       '%s' % TOP,
                       '* VELP: layer P velocities (km/s)',
                       '%s' % VELP,
                       '* RAT: p/s ratio ',
                       '%s' % RAT,
                       '*',
                       '*--- CLUSTER/EVENT SELECTION:',
                       '* CID',
                       '    0',
                       '* ID',
                       '*',
                       '* end of hypoDD.inp']

    # Write out hypoDD.inp in inputs directory
    with open(hypoDD_dir + 'INP/hypoDD.inp', "w") as open_file:
        open_file.write('\n'.join(hypoDD_inp_text))
    open_file.close()

    # Run hypoDD
    hypoDD_command = hypoDD_dir + 'HYPODD/src/hypoDD/hypoDD' + ' ' + hypoDD_dir + '/INP/hypoDD.inp'
    subprocess.call(hypoDD_command, shell=True)

    # Move hypoDD log and copy hypoDD outputs to relocate catalog output dir
    move_command = 'mv ./hypoDD.log %shypoDD.log' % relocate_catalog_output_dir
    copy_command = 'scp -r %sOUT_hypoDD/ %s' % (hypoDD_dir, relocate_catalog_output_dir)
    copy_command2 = 'scp -r %sOUT_ph2dt/ %s' % (hypoDD_dir, relocate_catalog_output_dir)
    subprocess.call(move_command, shell=True)
    subprocess.call(copy_command, shell=True)
    subprocess.call(copy_command2, shell=True)

    # Write out catalog objects from event.sel, hypoDD.loc and hypoDD.reloc
    print('Now writing catalog objects from hypoDD outputs...')
    hypoDD_loc = loc2cat(hypoDD_dir + 'OUT_hypoDD/hypoDD.loc', catalog, type='loc',
                         depth_correction=depth_correction)
    hypoDD_reloc = loc2cat(hypoDD_dir + 'OUT_hypoDD/hypoDD.reloc', catalog, type='reloc',
                           depth_correction=depth_correction)
    writer(relocate_catalog_output_dir + 'hypoDD_loc.xml', hypoDD_loc)
    writer(relocate_catalog_output_dir + 'hypoDD_reloc.xml', hypoDD_reloc)
    print('Done!')

    return hypoDD_loc, hypoDD_reloc


def plot_hypoDD_results(hypoDD_in,
                        hypoDD_loc,
                        hypoDD_reloc,
                        lat_lims,
                        lon_lims,
                        dep_lims,
                        markersize=5,
                        legend_loc='lower right',
                        export_filepath=None):
    """
    Plot hypoDD results in 3 dimensions, to show differences between the input catalog, relocation candidates, and relocated events
    :param hypoDD_in (:class:`~obspy.core.event.Catalog`): relocatable catalog
    :param hypoDD_loc (:class:`~obspy.core.event.Catalog`): catalog of all relocation candidate events (that have sufficient differential time observations)
    :param hypoDD_reloc (:class:`~obspy.core.event.Catalog`): catalog of all successfully relocated events
    :param lat_lims (list or tuple): minimum and maximum latitude plot limits expressed as a tuple/list of length 2
    :param lon_lims (list or tuple): minimum and maximum longitude plot limits expressed as a tuple/list of length 2
    :param dep_lims (list or tuple): minimum and maximum depth plot limits expressed as a tuple/list of length 2 (km)
    :param markersize (float): marker size for earthquake plots (defaults to 5)
    :param legend_loc (str): string describing position of legend in each subplot (e.g. 'lower right', 'lower left')
    :param export_filepath (str): (str or `None`): If str, exports plotted figure as a '.png' file. If `None`, show figure in interactive python.
    :return: N/A
    """

    from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper
    from matplotlib import pyplot as plt
    from toolbox import remove_catalog_repeats
    print('\nPlotting results...')

    # Get subset of relocation candidates that remain unshifted
    hypoDD_loc_filtered = remove_catalog_repeats(hypoDD_loc, hypoDD_reloc)

    # Extract lat, lon and depth for all catalogs
    hypoDD_in_lats = [event.origins[0].latitude for event in hypoDD_in]
    hypoDD_in_lons = [event.origins[0].longitude for event in hypoDD_in]
    hypoDD_in_deps = [event.origins[0].depth/1000 for event in hypoDD_in]
    hypoDD_loc_lats = [event.origins[0].latitude for event in hypoDD_loc]
    hypoDD_loc_lons = [event.origins[0].longitude for event in hypoDD_loc]
    hypoDD_loc_deps = [event.origins[0].depth/1000 for event in hypoDD_loc]
    hypoDD_loc_filtered_lats = [event.origins[0].latitude for event in hypoDD_loc_filtered]
    hypoDD_loc_filtered_lons = [event.origins[0].longitude for event in hypoDD_loc_filtered]
    hypoDD_loc_filtered_deps = [event.origins[0].depth/1000 for event in hypoDD_loc_filtered]
    hypoDD_reloc_lats = [event.origins[0].latitude for event in hypoDD_reloc]
    hypoDD_reloc_lons = [event.origins[0].longitude for event in hypoDD_reloc]
    hypoDD_reloc_deps = [event.origins[0].depth/1000 for event in hypoDD_reloc]

    # Now do simple plotting
    fig, ax = plt.subplots(2, 3, figsize=(10, 8))
    # [LON VS LAT] Plot selected events among input events (blue against black)
    ax[0, 0].plot(hypoDD_in_lons, hypoDD_in_lats, 'k.', markersize=markersize, label='input events')
    ax[0, 0].plot(hypoDD_loc_lons, hypoDD_loc_lats, 'b.', markersize=markersize, label='relocatable')
    ax[0, 0].set_aspect(2)
    ax[0, 0].axis([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]])
    ax[0, 0].legend(loc=legend_loc, fontsize=7.5)
    ax[0, 0].set_title('LON vs LAT (before hypoDD)')
    ax[0, 0].set_ylabel('Latitude')
    ax[0, 0].set_xlabel('Longitude')
    # [LON VS LAT] Plot relocated events among retained events (red against blue)
    ax[1, 0].plot(hypoDD_loc_filtered_lons, hypoDD_loc_filtered_lats, 'b.', markersize=markersize, label='relocatable')
    ax[1, 0].plot(hypoDD_reloc_lons, hypoDD_reloc_lats, 'r.', markersize=markersize, label='relocated')
    ax[1, 0].set_aspect(2)
    ax[1, 0].axis([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]])
    ax[1, 0].legend(loc=legend_loc, fontsize=7.5)
    ax[1, 0].set_title('LON vs LAT (after hypoDD)')
    ax[1, 0].set_ylabel('Latitude')
    ax[1, 0].set_xlabel('Longitude')
    # [LON VS DEP] Plot selected events among input events (blue against black)
    ax[0, 1].plot(hypoDD_in_lons, hypoDD_in_deps, 'k.', markersize=markersize, label='input events')
    ax[0, 1].plot(hypoDD_loc_lons, hypoDD_loc_deps, 'b.', markersize=markersize, label='relocatable')
    ax[0, 1].axis([lon_lims[0], lon_lims[1], dep_lims[0], dep_lims[1]])
    ax[0, 1].invert_yaxis()
    ax[0, 1].legend(loc=legend_loc, fontsize=7.5)
    ax[0, 1].set_title('LON vs DEP (before hypoDD)')
    ax[0, 1].set_ylabel('Depth (km)')
    ax[0, 1].set_xlabel('Longitude')
    # [LON VS DEP] Plot relocated events among retained events (red against blue)
    ax[1, 1].plot(hypoDD_loc_filtered_lons, hypoDD_loc_filtered_deps, 'b.', markersize=markersize, label='relocatable')
    ax[1, 1].plot(hypoDD_reloc_lons, hypoDD_reloc_deps, 'r.', markersize=markersize, label='relocated')
    ax[1, 1].axis([lon_lims[0], lon_lims[1], dep_lims[0], dep_lims[1]])
    ax[1, 1].invert_yaxis()
    ax[1, 1].legend(loc=legend_loc, fontsize=7.5)
    ax[1, 1].set_title('LON vs DEP (after hypoDD)')
    ax[1, 1].set_ylabel('Depth (km)')
    ax[1, 1].set_xlabel('Longitude')
    # [LAT VS DEP] Plot selected events among input events (blue against black)
    ax[0, 2].plot(hypoDD_in_lats, hypoDD_in_deps, 'k.', markersize=markersize, label='input events')
    ax[0, 2].plot(hypoDD_loc_lats, hypoDD_loc_deps, 'b.', markersize=markersize, label='relocatable')
    ax[0, 2].axis([lat_lims[0], lat_lims[1], dep_lims[0], dep_lims[1]])
    ax[0, 2].invert_yaxis()
    ax[0, 2].legend(loc=legend_loc, fontsize=7.5)
    ax[0, 2].set_title('LAT vs DEP (before hypoDD)')
    ax[0, 2].set_ylabel('Depth (km)')
    ax[0, 2].set_xlabel('Latitude')
    # [LAT VS DEP] Plot relocated events among retained events (red against blue)
    ax[1, 2].plot(hypoDD_loc_filtered_lats, hypoDD_loc_filtered_deps, 'b.', markersize=markersize, label='relocatable')
    ax[1, 2].plot(hypoDD_reloc_lats, hypoDD_reloc_deps, 'r.', markersize=markersize, label='relocated events')
    ax[1, 2].axis([lat_lims[0], lat_lims[1], dep_lims[0], dep_lims[1]])
    ax[1, 2].invert_yaxis()
    ax[1, 2].legend(loc=legend_loc, fontsize=7.5)
    ax[1, 2].set_title('LAT vs DEP (after hypoDD)')
    ax[1, 2].set_ylabel('Depth (km)')
    ax[1, 2].set_xlabel('Latitude')
    plt.tight_layout()

    if export_filepath:
        fig.savefig(export_filepath, bbox_inches='tight')
        plt.close()
    else:
        fig.show()

    print('Done!')