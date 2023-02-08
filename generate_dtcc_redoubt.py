#%% GENERATE DTCC

# Import all dependencies
import os
import numpy as np
from eqcorrscan.utils.catalog_to_dd import write_correlations, _generate_event_id_mapper
from obspy import Catalog
from toolbox import reader, prepare_stream_dict, adjust_weights, estimate_s_picks

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/' # '/Users/darrentpk/Desktop/GitHub/enhance_catalog/' #
data_dir = '/home/data/redoubt/'
output_dir = main_dir + 'output/redoubt5/'
scan_data_output_dir = output_dir + 'scan_data/'
relocatable_catalog_filename = 'relocatable_catalog_FImag.xml'
by_cluster = False
make_s_picks = False

length_actual = 2.66    # even-numbered decimal so that 50Hz sampling rate won't have incorrect trimming lengths
length_excess = 4       # in excess for stream_dict
pre_pick_actual = 0.15  # shorter than template
pre_pick_excess = 1     # in excess for stream_dict
shift_len = 0.4         # width of search for max_cc
lowcut = 1              # same as template
highcut = 10            # same as template

max_sep = 50             # max separation tolerated (8km)
min_link = 3            # minimum number of matching pick stachans (4)
for min_cc in [0.5, 0.6, 0.7]:
    # min_cc = 0.4            # minimum cc to be considered pick
    weight_func = None      # input function to be applied on raw cc coefficient squared

    #%% Generate a stream dictionaries using prepare_stream_dict and execute cross correlations

    # Load in catalog
    catalog = reader(scan_data_output_dir + relocatable_catalog_filename)
    if make_s_picks:
        catalog = estimate_s_picks(catalog)

    # If we are writing correlations in bulk of the whole catalog, we process the entire catalog at one go
    if not by_cluster:

        # Generate stream dictionary (refer to toolbox.py)
        # Note that we want pre_pick and length to be in excess, since write_correlations trims the data for us
        stream_dict = prepare_stream_dict(catalog,pre_pick=pre_pick_excess,length=length_excess,local=True,data_dir=data_dir)

        # Execute cross correlations and write out a .cc file using write_correlations (refer to EQcorrscan docs)
        # Note this stores a file called "dt.cc" in your current working directory
        event_id_mapper = write_correlations(catalog=catalog, stream_dict=stream_dict, extract_len=length_actual,
                                             pre_pick=pre_pick_actual, shift_len=shift_len, lowcut=lowcut,
                                             highcut=highcut, max_sep=max_sep, min_link=min_link, min_cc=min_cc,
                                             interpolate=False, max_workers=None, parallel_process=True)

        # Define source and target cc filepaths
        original_dt_dir = os.getcwd() + '/dt.cc'
        target_dt_dir = os.getcwd() + '/master_dt_rdt_all_' + str(int(min_cc*10)) + '_S.cc'
        adjust_weights(original_dt_dir, target_dt_dir, dt_cap=shift_len, min_link=min_link, append=False, weight_func=weight_func)

    # If we are correlating by cluster
    else:

        # Start by generating event id mapper for full catalog
        event_id_mapper = _generate_event_id_mapper(catalog, event_id_mapper=None)

        # Obtain a unique list of source templates
        templates = []
        for event in catalog:
            templates.append(event.comments[0].text.split()[1])
        unique_templates = np.unique(templates)

        # Initialize a list to store the catalog index of every template's self-detection
        template_indices = []

        # Define source and target cc filepaths
        original_dt_dir = os.getcwd() + '/dt.cc'
        target_dt_dir = os.getcwd() + '/master_dt_rdt_clust_' + str(int(min_cc*10)) + '_S.cc'

        # Loop through each unique template
        for i, unique_template in enumerate(unique_templates):

            # Find index of catalog events that correspond to this template
            template_detection_indices = [k for k, template in enumerate(templates) if unique_template in template]

            # Craft sub-catalog
            sub_catalog = Catalog()
            detect_vals = []
            for template_detection_index in template_detection_indices:
                template_detection = catalog[template_detection_index]
                detect_val = float(template_detection.comments[2].text.split('=')[1])
                detect_vals.append(detect_val)
                sub_catalog += template_detection

            # Store the index of the template's self-detection, this will be used for our final inter-cluster step
            template_indices.append(template_detection_indices[np.argmax(detect_vals)])

            # Now craft stream dictionary
            stream_dict = prepare_stream_dict(sub_catalog, pre_pick=pre_pick_excess, length=length_excess, local=True, data_dir=data_dir)

            # Execute cross correlations
            _ = write_correlations(catalog=sub_catalog, stream_dict=stream_dict, event_id_mapper=event_id_mapper,
                                   extract_len=length_actual, pre_pick=pre_pick_actual, shift_len=shift_len,
                                   lowcut=lowcut, highcut=highcut, max_sep=max_sep, min_link=min_link,
                                   min_cc=min_cc, interpolate=False, max_workers=None, parallel_process=True)

            # Write/append dt.cc to target cc file
            if i == 0:
                adjust_weights(original_dt_dir, target_dt_dir, dt_cap=shift_len, min_link=min_link, append=False, weight_func=weight_func)
            else:
                adjust_weights(original_dt_dir, target_dt_dir, dt_cap=shift_len, min_link=min_link, append=True, weight_func=weight_func)

        # Now execute inter-template cross correlation by generating a sub-catalog containing template self-detections
        sub_catalog = Catalog()
        for template_index in template_indices:
            sub_catalog += catalog[template_index]

        # Now craft stream dictionary
        stream_dict = prepare_stream_dict(sub_catalog, pre_pick=pre_pick_excess, length=length_excess, local=True,
                                          data_dir=data_dir)

        # Execute cross correlations
        _ = write_correlations(catalog=sub_catalog, stream_dict=stream_dict, event_id_mapper=event_id_mapper,
                               extract_len=length_actual, pre_pick=pre_pick_actual, shift_len=shift_len,
                               lowcut=lowcut, highcut=highcut, max_sep=max_sep, min_link=min_link,
                               min_cc=min_cc, interpolate=False, max_workers=None, parallel_process=False)

        # Append dt.cc to master dt.cc
        adjust_weights(original_dt_dir, target_dt_dir, dt_cap=shift_len, min_link=min_link, append=True, weight_func=None)

    # # Write event_id_mapper
    # with open(relocate_catalog_output_dir + 'event_id_mapper.pkl', 'wb') as evid_pickle:  # Pickling
    #     pickle.dump(event_id_mapper, evid_pickle)

    print('Cross correlations done!')