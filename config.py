#%% CONFIG
# This text file stores all necessary inputs across our modules

#%% convert_redpy.py
hypoi_file = 'redoubt_20090101_20090501_hypoi.txt'
hypoddpha_file = 'redoubt_20090101_20090501_hypoddpha.txt'
max_dt = 4  # maximum time difference between REDPy detections and AVO events allowed, in seconds
adopt_weight = 0.1  # phase weight for adopted picks
redpy_station_list = ['RDN','REF','RSO']  # note that current code that makes picks look at EHZ only
data_dir = '/home/data/redoubt/'  # redoubt data directory
redpy_dir = '/home/ptan/attempt_eqcorrscan/'
main_dir = ''