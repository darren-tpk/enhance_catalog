# import glob
# from obspy import read
# import os
#
# campaign_dir = '/home/ptan/enhance_catalog/data/redoubt_campaign/'
# campaign_files = glob.glob(campaign_dir + '*')
#
# campaign_output_dir = '/home/data/redoubt/'
#
# for campaign_file in campaign_files:
#
#     # if it is BH, read file
#     st = read(campaign_file)
#     seed_filepath = campaign_output_dir + campaign_file.split('/')[-1]
#     st.write(seed_filepath, format="MSEED")
#     os.remove(campaign_file)
#
