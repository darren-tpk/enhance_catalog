# [calculate_relative_magnitudes] function to calculate magnitudes using CC based on Schaff & Richards (2014)
# Uses template events in PEC as reference magnitudes for relocatable catalog
def calculate_relative_magnitudes(catalog, tribe, data_dir, noise_window, signal_window, min_cc,
                                  shift_len, tolerance,resampling_frequency):

    # Import dependencies
    import numpy as np
    from eqcorrscan.utils.mag_calc import relative_magnitude, _get_signal_and_noise, _get_pick_for_station
    from toolbox import reader, prepare_catalog_stream
    from obspy import Catalog
    from obspy.core.event import Magnitude

    # Derive template names from tribe
    template_names = [template.name for template in tribe]

    # Loop through events in target catalog
    for target_event in catalog:

        # Obtain template event
        template_index = template_names.index(target_event.comments[0].text.split(' ')[1])
        template_event = tribe[template_index].event

        # Obtain day-long streams for both the target event and the template event. (Also detrend and filter)
        target_st = prepare_catalog_stream(data_dir,Catalog() + target_event,resampling_frequency,tolerance)
        target_st = target_st.split()
        target_st = target_st.detrend('simple')
        target_st = target_st.filter(type='bandpass',freqmin=1,freqmax=10,corners=4)
        target_st = target_st.merge()
        template_st = prepare_catalog_stream(data_dir,Catalog() + template_event,resampling_frequency,tolerance)
        template_st = template_st.split()
        template_st = template_st.detrend('simple')
        template_st = template_st.filter(type='bandpass',freqmin=1,freqmax=10,corners=4)
        template_st = template_st.merge()

        # Calculate magnitude differences weighted by cc
        delta_mags, ccs = relative_magnitude(target_st, template_st, target_event, template_event,
                                             noise_window=noise_window, signal_window=signal_window,
                                             min_snr=min_snr, min_cc=min_cc, use_s_picks=False, correlations=None,
                                             shift=shift_len, return_correlations=True, correct_mag_bias=True)

        # Now calculate target magnitude
        delta_mag_values = [delta_mag[1] for delta_mag in delta_mags.items()]
        cc_values =  [cc[1] for cc in ccs.items() if (cc[1] > min_cc)]
        original_mag = template_event.preferred_magnitude().mag
        target_mag = original_mag + np.sum(delta_mag_values) / np.sum(cc_values)

        # Add target magnitude into input catalog
        target_event.magnitudes.append(Magnitude(mag=target_mag))

    # Return catalog
    return catalog
