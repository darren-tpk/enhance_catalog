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
    print('Creating subdirectories for workflow outputs...')
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
                    seed_filename = sta + '.' + chan + '.' + trace_datetime_str

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
        print('Data collection for the day complete. Elapsed time: %.2f hours' % ((time_current - time_start) / 3600))

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
              trigalg='classicstalta',
              offset=0,
              plotformat='eqrate,fi,occurrence+occurrencefi,longevity',
              printsta=0,
              minplot=5,
              dybin=1,
              hrbin=1.,
              occurbin=1,
              recbin=1,
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
    :param trigalg (str): chosen STA/LTA algorithm supported by ObsPy's network coincidence trigger (defaults to 'classicstalta')
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
    initialize_command = 'python ./redpy/initialize.py -v -c %s' % (output_destination + cfg_name)
    subprocess.call(initialize_command, shell=True)

    # Call backfill command
    initialize_command = 'python ./redpy/backfill.py -v -c %s -s %s -e %s' % (output_destination + cfg_name,
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
    from toolbox import reader, writer, pull_cores, read_trace

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
            PEC_index = np.argmin(abs(analyst_event_times - detection_time))
            analyst_event = analyst_catalog[PEC_index]

            # Use the closest event's information to fill up event object
            redpy_event = analyst_event.copy()
            full_event_tag = event_tag + ' in analyst catalog, REDPy detection time = %s' % str(detection_time)
            redpy_event.origins[0].comments.append(Comment(text=full_event_tag))

            # Add cluster number to valid clusters list
            associated_cluster_list.append(cluster)

            # Remove avo index from unmatched list if it still exists
            if PEC_index in unmatched_indices_redpy:
                unmatched_indices_redpy.remove(PEC_index)
            if PEC_index in unmatched_indices_core and event_tag.split(' ')[2] == 'core':
                unmatched_indices_core.remove(PEC_index)

        # If event is not part of the PEC, we tag it as so
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
    np.save(convert_redpy_output_dir + 'associated_clusters.npy', unassociated_clusters)
    np.save(convert_redpy_output_dir + 'unmatched_indices_redpy.npy', unassociated_clusters)
    np.save(convert_redpy_output_dir + 'unmatched_indices_core.npy', unassociated_clusters)

    # Also generate unmatched analyst catalogs (that can be used as templates later)
    unmatched_analyst_events_redpy = Catalog() # PEC events that don't match redpy catalog
    unmatched_analyst_events_core = Catalog() # PEC events that don't match redpy cores
    for j, analyst_event in enumerate(analyst_catalog):
        if j in unmatched_indices_redpy:
            unmatched_analyst_events_redpy += analyst_event
        if j in unmatched_indices_core:
            unmatched_analyst_events_core += analyst_event

    # Write out unmatched PEC catalog to .xml file
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
                station_tr = read_trace(data_dir=data_path, station=redpy_station[k], channel=redpy_channel[k],
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
                coin_trigger = coincidence_trigger('classicstalta', trigon, trigoff, stream, nstaC, sta=swin, lta=lwin, details=True)
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
                    station_tr = read_trace(data_dir=data_path, station=campaign_station[k], channel=campaign_channel[k],
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
                    coin_trigger = coincidence_trigger('classicstalta', trigon, trigoff, stream, nstaC, sta=swin, lta=lwin, details=True)
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
                 samprate,
                 fmin,
                 fmax,
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
    :param samprate (float): desired sampling rate for templates
    :param fmin (float): bandpass filter lower bound (Hz)
    :param fmax (float): bandpass filter upper bound (Hz)
    :param prepick (float): time before pick time to start template waveform trim (s)
    :param length (float): time from pre-pick to stop template waveform trim (s)
    :param min_snr (float): minimum signal-to-noise to accept waveform into template
    :param process_len (float): maximum data length (in seconds) to process at once (defaults to 86400 s)
    :param tolerance (float): factor to median tolerance for boxcar removal from data (as a factor to median)
    :param channel_convention (bool): if `True`, enforce strict compliance for P/S picks to be on vertical/horizontal components
    :param use_all_analyst_events (bool): if `True`, disregard REDPy clustering and create templates from all analyst-derived events
    :return:
    """

    # Import all dependencies
    import time
    import numpy as np
    from itertools import compress
    from obspy import UTCDateTime, Catalog
    from eqcorrscan.core.match_filter.tribe import Tribe
    from toolbox import prepare_catalog_stream, reader, writer

    # If template_stations is a comma separated string, convert to list
    if type(template_stations) == list:
        template_stations = template_stations.split(',')

    # Read in, and combine, the picked core catalog and unmatched PEC catalog
    core_catalog_picked = reader(convert_redpy_output_dir + 'core_catalog_picked.xml')
    if use_all_analyst_events:
        unmatched_PEC_events = reader(convert_redpy_output_dir + 'unmatched_PEC_events_core.xml')
    else:
        unmatched_PEC_events = reader(convert_redpy_output_dir + 'unmatched_PEC_events_redpy.xml')
    template_catalog = core_catalog_picked + unmatched_PEC_events

    # Clean catalog to only include picks from our station list
    print('Removing picks on stations not in input station list:')
    for i, event in enumerate(template_catalog):
        for pick in reversed(event.picks):
            if pick.waveform_id.station_code not in template_stations:
                print('Removed ' + pick.waveform_id.station_code + ' pick from template catalog index ' + str(i))
                event.picks.remove(pick)

    # Initialize tribe and tracker for valid events
    tribe = Tribe()
    valid_event = []
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
            sub_catalog = Catalog(list(compress(template_catalog, sub_catalog_index)))
            tracker[np.where(sub_catalog_index==True)] = True

        # Prepare catalog stream, trim to the start and end of the day to enforce daylong processing
        stream = prepare_catalog_stream(data_path, sub_catalog, samprate, tolerance)
        stream = stream.trim(starttime=event_daystart, endtime=event_dayend, pad=True)

        # Do stream checks
        for trace in stream:
            # Check if trace is masked with too many zeros (more than half of samples are masked)
            if hasattr(trace.data, 'mask') and (np.sum(trace.data.mask) / len(trace.data.mask)) > 0.5:
                print('%s.%s got removed due to overly large data gaps.' % (trace.stats.station, trace.stats.channel))
                stream.remove(trace)

        # Construct sub tribe out of sub catalog
        try:
            sub_tribe = Tribe().construct(
                method="from_meta_file", lowcut=fmin, highcut=fmax, samp_rate=samp_rate, length=length,
                filt_order=4, prepick=prepick, meta_file=sub_catalog, st=stream, process=True,
                process_len=process_len, min_snr=min_snr, parallel=True)
            if sub_tribe is not None or len(sub_tribe) > 0:
                tribe += sub_tribe

        # If tribe fails then we add the current sub catalog to the client catalog to try the client method later
        except:
            print('Tribe creation attempt on %s to %s failed.')

    # Conclude process
    time_end = time.time()
    print('\nTemplate creation attempts complete. Time elapsed: %.2f s' % (time_end-time_start))

    print('%d out of %d events in the catalog were converted to templates.' % (len(tribe),len(template_catalog)))
    writer(create_tribe_output_dir + 'tribe.tgz', tribe)