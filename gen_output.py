# generate data products

# import packages and functions that we need
import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from eqcorrscan import Party, Family, Template, Detection
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from toolbox import remove_repeats

from eqcorrscan import Party
party_in = '/home/ptan/attempt_eqcorrscan/output/party.tgz'
party = Party().read(party_in)
# party2_in = '/home/ptan/project/output/party2.tgz'
# party2 = Party().read(party2_in)

# clean party off of repeated detections
time_interval = 30  # number of seconds before and after each detection to check for repeats
party_clean = remove_repeats(party,time_interval)
# save cleaned party
party_clean_out = 'party_clean'
party_clean_outpath = '/home/ptan/attempt_eqcorrscan/output/' + party_clean_out
if os.path.exists(party_clean_outpath+'.tgz'):
    os.remove(party_clean_outpath+'.tgz')
    party_clean.write(party_clean_outpath + '.tgz' , format='tar')
else:
    party_clean.write(party_clean_outpath + '.tgz' , format='tar')