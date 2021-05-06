# generate data products

# import packages and functions that we need
from toolbox import remove_repeats, reader, writer

# define all variables here
time_interval = 5  # number of seconds before and after each detection to check for repeats

# read tribe data
party = party_all
# party = reader('/home/ptan/attempt_eqcorrscan/output/party_marchswarm.tgz')

# derive a party that only has associated cores


# clean party off of repeated detections
party_clean = remove_repeats(party,time_interval)



# save cleaned party
party_clean_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_clean_outpath+'party_clean.tgz',party_clean)