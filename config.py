#%% CONFIG
# This text file stores all necessary inputs across our modules

#%% Import all dependencies
from obspy import UTCDateTime

#%% Global variables and filepaths
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/redoubt/'  # redoubt data directory on local
output_dir ='/home/ptan/enhance_catalog/output/'
redpy_results_dir = '/home/ptan/enhance_catalog/redpy_results/redoubt/'
PEC_dir = '/home/ptan/enhance_catalog/data/avo/'
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
sitelist_dir = '/home/ptan/enhance_catalog/data/avo/'

convert_redpy_output_dir = output_dir + 'convert_redpy/'
create_tribe_output_dir = output_dir + 'create_tribe/'
scan_data_output_dir = output_dir + 'scan_data/'
relocate_catalog_output_dir = output_dir + 'relocate_catalog/'
plot_hypo_output_dir = output_dir + 'plot_hypo/'

#%% convert_redpy.py
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks
redpy_station_list = ['RDN','REF','RSO']
redpy_channel_list = [['EHZ'],['EHZ'],['EHZ']]
tolerance = 4e4  # tolerance for boxcar removal from data (as a factor to median)
                 # used to be 5e4?

#%% create_tribe.py
tribe_filename = 'tribe_test.tgz'
channel_convention = True  # strict compliance for P/S picks on vertical/horizontal components
resampling_frequency = 50  # for resampling traces prior to final merge
tolerance = 4e4            # tolerance for boxcar removal from data (as a factor to median)
lowcut = 1.0               # filter lowcut (Hz), follows Jeremy's recommendation
highcut = 10.0             # filter highcut (Hz), follows Jeremy's recommendation
samp_rate = 50.0           # new sampling rate (Hz)
length = 10.0              # template length (s), Wech et al. (2018) chose 30s
filt_order = 4             # number of corners for filter
prepick = 1.0              # pre-pick time (s), Wech et al. (2018) chose 5s
process_len = 86400        # length to process data in (s)
min_snr = 3             # minimum SNR, Jeremy's recommendation was 5.0 (same as EQcorrscan tutorial)
local_volcano = 'redoubt'  # for get_local_stations function, since we only take picks from stations near Redoubt
local_radius = 25          # for get_local_stations function; radius around volcano to accept stations

#%% scan_data.py
tribe_filename = 'tribe_test.tgz'
party_filename = 'party_test.tgz'
catalog_filename = 'party_catalog_test.xml'
repicked_catalog_filename = None
min_stations = 3                               # to remove templates that are anchored by too little stations
min_picks = 0                                  # to remove templates that are anchored by too little picks
start_time = UTCDateTime(2009, 1, 1, 0, 0, 0)  # start: UTCDateTime(2009, 1, 1, 0, 0, 0)
end_time = UTCDateTime(2009, 5, 1, 0, 0, 0)    # goal: UTCDateTime(2009, 5, 1, 0, 0, 0)
samp_rate = 50                                 # to resample streams to match tribes
tolerance = 4e4                                # for boxcar removal
threshold_type = 'av_chan_corr'                # also consider 'MAD', Jeremy uses 12
threshold = 0.60                               # used to be 0.74 from sensitivity test. Was too high
trig_int = 10                                  # trigger interval for each template. Also used to remove repeats (decluster)
parallel_process = 'True'                      # parallel process for speed
generate_repicked_catalog = False              # option to use lag_calc to do catalog repicking

#%% relocate_catalog.py
tribe_filename = 'tribe_test.tgz'
catalog_filename = 'party_catalog_test.xml'
relocated_catalog_filename = 'relocated_catalog_test.xml'
raw_station_list_dir = main_dir + 'data/avo/'
raw_station_list_filename = 'station_list.csv'
raw_vzmodel_dir = main_dir + 'data/avo/'
raw_vzmodel_filename = 'redoubt_vzmodel.txt'
growclust_exe = main_dir + '/growclust/SRC/growclust'
length_actual = 10      # same as template
length_excess = 30      # in excess for stream_dict
pre_pick_actual = 2     # same as template
pre_pick_excess = 10    # in excess for stream_dict
shift_len = 3           # width of search for max_cc
lowcut = 1              # same as template
highcut = 10            # same as template
max_sep = 8             # max separation tolerated (8km)
min_link = 3            # minimum number of matching pick stachans (3)
min_cc = 0.4            # minimum cc to be considered pick