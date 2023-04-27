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
                   'filename='+h5_name,
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


