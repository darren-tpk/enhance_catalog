import os
import subprocess

directory = '/home/ptan/enhance_catalog/data/mammoth/'

old_filenames = os.listdir(directory)

for old_filename in old_filenames:

    old_filepath = directory + old_filename

    pieces = old_filename.split('.')
    new_filename = '.'.join([pieces[1],pieces[0],pieces[2]])
    new_filepath = directory + new_filename

    rename_command = 'mv ' + old_filepath + ' ' + new_filepath
    subprocess.call(rename_command, shell=True)




## Fix Great Sitkin cc file
import pandas as pd
from obspy import Catalog, UTCDateTime
from obspy.core.event import Event, Origin
from toolbox import writer
import subprocess

cc_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/output/relocate_catalog/IN/'
old_cc_filepath = cc_dir + 'xcordata1.txt'
new_cc_filepath = cc_dir + 'xcordata2.txt'
new_cc2_filepath = cc_dir + 'xcordata.txt'
old_cc = pd.read_csv(old_cc_filepath,header=None)
new_cc = open(new_cc_filepath, 'w')
for i in range(len(old_cc.values)):
    text = old_cc.values[i][0]
    if text[0] == '#':
        new_cc.write(text + '\n')
    else:
        sta,dt,wt,pha = text.split()
        if abs(float(dt)) > 5:
            continue
        else:
            new_cc.write(text + '\n')
new_cc.close()

new_cc = pd.read_csv(new_cc_filepath,header=None)
new_cc2 = open(new_cc2_filepath, 'w')
for i in range(len(new_cc.values)-1):
    if new_cc.values[i][0] == '#' and new_cc.values[i+1][0] == '#':
        continue
    else:
        new_cc2.write(new_cc.values[i][0] + '\n')
new_cc2.close()

growclust_command = '/Users/darrentpk/Desktop/Github/enhance_catalog/growclust/SRC/growclust /Users/darrentpk/Desktop/Github/enhance_catalog/output/relocate_catalog/config.inp'
subprocess.call(growclust_command, shell=True)


relocated_event_table = pd.read_csv('/Users/darrentpk/Desktop/Github/enhance_catalog/output/relocate_catalog/OUT/out.growclust_cat', header=None, sep='\s+',
                                    names=['yr','mon','day','hr','mm','sec','evid','lat','lon','dep','mag','qID','cID',
                                           'nbranch','qnpair','qndiffP','qndiffS','rmsP','rmsS','eh','ez','et','latC',
                                           'lonC','depC'])

# Initialize ObsPy catalog object and populate with events on a row-by-row basis
relocated_catalog = Catalog()
for i in range(len(relocated_event_table)):
    time = UTCDateTime(relocated_event_table.yr[i],relocated_event_table.mon[i],relocated_event_table.day[i],relocated_event_table.hr[i],relocated_event_table.mm[i],relocated_event_table.sec[i])
    event = Event(origins=[Origin(time=time,latitude=relocated_event_table.lat[i],longitude=relocated_event_table.lon[i],depth=relocated_event_table.dep[i]*1000)])
    relocated_catalog += event
writer('/Users/darrentpk/Desktop/Github/enhance_catalog/output/relocate_catalog/relocated_catalog_GS.xml', relocated_catalog)
