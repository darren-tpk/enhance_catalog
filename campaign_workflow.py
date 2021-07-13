# import dependencies
import os
import time
from toolbox import reader
from obspy import read
from obspy.signal.trigger import coincidence_trigger, classic_sta_lta, trigger_onset

# set up catalog of picks during campaign period
catalog_inpath = '/home/ptan/attempt_eqcorrscan/output/'
core_catalog = reader(catalog_inpath+'redpy_catalog_picked.xml')
unmatched_avo_catalog = (catalog_inpath+'unmatched_avo_events.xml')
catalog = core_catalog + unmatched_avo_catalog

# prepare for pick-making
campaign_stations = ['RD01','RD02','RD03','RDW']
campaign_channels = ['BHZ'] # try P picks only to start. Maybe use ['BHZ','BHE','BHN'] in the future
campaign_dir = '/home/admin/databases/REDOUBT_2009/wf/data/campaign/2009'
campaign_juldays = [x[0] for x in os.listdir(campaign_dir)]

# make picks on campaign stations
print('\nMaking additional picks on campaign stations for existing events...')
time_start = time.time()
# loop through events in current catalog
for event in catalog:
    # extract event time and convert to julday to point at correct directory in database
    event_time = event.origins[0].time
    event_julday = '%03d' % event_time.julday
    # proceed if this julday falls within campaign duration
    if event_julday in campaign_juldays:
        # now load each stachan data and do picking
        for campaign_station in campaign_stations:
            for campaign_channel in campaign_channels:
                data_filename = 'AV_' + campaign_station + '_' + campaign_channel + '__2009-' + event_julday
                data_filenames.append(data_filename)



for cluster_NA in clusters_NA:
    # find all associated detections
    contribution_list = list(np.where(np.array(redpy_det.Cluster) == cluster_NA)[0])
    # explicitly loop through associated detections
    for contribution_index in contribution_list:
        # determine time difference to offset pick times
        contribution_time = redpy_catalog[contribution_index].origins[0].time
        # load data from local machine in a +/- 12s window
        starttime = contribution_time - 12
        endtime = contribution_time + 12
        # remove RSO after first explosion [HARD CODED]
        redpy_stations = redpy_station_list
        if starttime > UTCDateTime(2009,3,24,0,0,0) and starttime < UTCDateTime(2009,4,16,0,0,0):
            redpy_stations = ['RDN', 'REF']
        # gather waveforms and filter, taper
        stream = Stream()
        for redpy_station in redpy_stations:
            station_tr = read_trace(data_dir=data_dir, station=redpy_station, channel='EHZ',
                                    starttime=starttime , endtime=endtime, tolerance=5e4)
            stream = stream + station_tr
        stream = stream.split()
        stream = stream.filter('bandpass',freqmin=1.0, freqmax=10.0, corners=2, zerophase=True)
        stream = stream.taper(0.05, type='hann', max_length=(0.75*1024/100))  # [HARD CODED]
        stream = stream.merge(method=1, fill_value=0)
        # use coincidence trigger to get pick time estimate
        coin_trigger = coincidence_trigger('classicstalta', 3, 2, stream, 2,
                                      sta=0.7, lta=8, details=True)
        # if there are no coincidence triggers, move to next event
        if coin_trigger == []:
            continue
        # otherwise, continue to single channel STA LTA
        else:
            # extract coincidence trigger time
            coin_trigger_time = coin_trigger[0]["time"]
            # for each channel
            for tr in stream:
                # calculate the value of the characteristic function
                sampling_rate = tr.stats.sampling_rate
                cft = classic_sta_lta(tr.data, int(0.7 * sampling_rate), int(8 * sampling_rate))
                #plot_trigger(tr, cft, 3, 2)
                # obtain trigger limits
                trigger_limits = np.array(trigger_onset(cft, 3, 2))
                # if there exists some trigger limits
                if trigger_limits.size != 0:
                    # convert to UTCDateTime and find the trigger on time closest to the coincidence trigger
                    trigger_on = np.array([tr.stats.starttime + t for t in (trigger_limits[:, 0] / sampling_rate)])
                    pick_time = trigger_on[np.argmin(abs(trigger_on-coin_trigger_time))]
                    # craft waveform stream ID
                    tr_id = tr.id.split('.')
                    waveform_id = WaveformStreamID(tr_id[0],tr_id[1],tr_id[2],tr_id[3])
                    # create pick and arrival objects and add to redpy_catalog
                    add_pick = Pick(time=pick_time,waveform_id=waveform_id,phase_hint='P')
                    add_arrival = Arrival(pick_id=ResourceIdentifier(id=add_pick.resource_id),phase='P',time_weight=0.1)
                    redpy_catalog[contribution_index].picks.append(add_pick)
                    redpy_catalog[contribution_index].origins[0].arrivals.append(add_arrival)
time_stop = time.time()
print('Pick making complete, processing time: %.2f s' % (time_stop-time_start))
writer(catalog_outpath+'redpy_catalog_picked.xml', redpy_catalog)
print('Picked redpy catalog saved.')