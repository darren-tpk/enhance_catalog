# generate data products

# import packages and functions that we need
import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from eqcorrscan import Party, Family, Template, Detection
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha
from toolbox import remove_repeats, reader, writer

from eqcorrscan import Party
party = reader('/home/ptan/attempt_eqcorrscan/output/party.tgz')
# party2 = reader('/home/ptan/attempt_eqcorrscan/output/party2.tgz')

# clean party off of repeated detections
time_interval = 30  # number of seconds before and after each detection to check for repeats
party_clean = remove_repeats(party,time_interval)
# save cleaned party
party_clean_outpath = '/home/ptan/attempt_eqcorrscan/output/'
writer(party_clean_outpath+'party_clean.tgz')