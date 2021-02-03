def build_spectrogram(event,data_pad,process_len)

    from obspy.clients.fdsn import Client
    from obspy import Stream
    from obspy import UTCDateTime
    from waveform_collection import gather_waveforms

    num_picks = len(event.picks)

    STARTTIME = event.origins[0].time
    ENDTIME = event.origins[0].time + 30

    for i in range(num_picks):
        print(i)
        pick = event.picks[i]
        network = pick.waveform_id.network_code
        station = pick.waveform_id.station_code
        channel = pick.waveform_id.channel_code
        location = pick.waveform_id.location_code
        st = gather_waveforms(source='AVO', network=network, station=station, channel=channel,
                              location='*', starttime=STARTTIME, endtime=ENDTIME)
