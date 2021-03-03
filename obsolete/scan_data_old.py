# Scan data to make detections

# Import packages and functions that we need
import os
import pickle
import numpy as np
from eqcorrscan import Party
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Load tribe data
save_dir = '/home/ptan/project/output/'
tribe_file = save_dir + 'tribe.pickle'
tribe_in = open(tribe_file,'rb')
tribe = pickle.load(tribe_in)
print(tribe)

# print station numbers of each template before filter
print('\nBefore filter:')
for i in range(len(tribe)):
     t = tribe[i]
     num_stations = len({tr.stats.station for tr in t.st})
     print('Template ' + str(i+1) + ': ' + str(num_stations) + ' stations')
# filter away templates with < 3 stations
tribe.templates = [t for t in tribe if len({tr.stats.station for tr in t.st}) >= 3]
print('\nAfter removing templates with < 3 stations...')
# print station numbers of each template again
for i in range(len(tribe)):
     t = tribe[i]
     num_stations = len({tr.stats.station for tr in t.st})
     print('Template ' + str(i+1) + ': ' + str(num_stations) + ' stations')

# BRUTE FORCE SCAN
# detect events across desired duration, in time steps to cope with memory limitations
client = Client("IRIS")
run_time = UTCDateTime()
scan_time_start = UTCDateTime(2009,2,15,21,0,0)
scan_time_end = UTCDateTime(2009,3,15,21,0,0)
num_time_steps = 28
step_duration = (scan_time_end - scan_time_start)/num_time_steps
# Initialize detection counter for outfile numbering
resume_boolean = int(input('Start new detection counter (0) or resume detection (1): '))
if resume_boolean == 0:
    detection_counter = np.zeros(len(tribe))
elif resume_boolean == 1:
    detection_counter = np.zeros(len(tribe))
    all_files = os.listdir(save_dir+'examples')
    for i in range(len(tribe)):
        num_detections = []
        for j in range(len(all_files)):
            file_string = all_files[j]
            if ('temp'+str(i+1)+'_') in file_string:
                file_string = all_files[j]
                num_detections.append(int(file_string.split('det')[-1].split('.')[0]))
        detection_counter[i] = max(num_detections)
else:
    raise ValueError('Invalid value for resumption boolean')
# commence loop for time steps
resume_timesteps = int(input('Start from time step number (>0): '))
if resume_timesteps < 1 or resume_timesteps > num_time_steps:
    raise ValueError('Invalid time step number entered')
print('\nTotal number of time steps to be scanned: ' + str(num_time_steps-resume_timesteps+1))
for i in range(resume_timesteps-1,num_time_steps):
    # print statement to keep tabs on time step
    print('\nScanning time step ' + str(i+1) + '...')
    # define t1 and t2 for current time step
    t1 = scan_time_start + (i*step_duration)
    t2 = scan_time_start + ((i+1)*step_duration)
    # detect events in step duration
    party, st = tribe.client_detect(
        client=client, starttime=t1, endtime=t2, threshold=20.,
        threshold_type="MAD", trig_int=30.0, plot=False, return_stream=True)
    # check if length of tribe matches number of families
    if len(tribe) != len(party.families):
        print('WARNING: Length of tribe does not match number of families. timestep = ',str(i+1))
    # commence a nested loop to extract detection value and store templates and detections
    print('Extracting templates and detections...')
    for j in range(len(party.families)):
        # extract family
        family = party[j]
        # plot template and save fig
        outfile = save_dir + 'temp' + str(j+1) + '.png'
        fig = family.template.st.plot(equal_scale=False, size=(800, 600), outfile=outfile)
        # extract streams for detections
        streams = family.extract_streams(stream=st, length=30.0, prepick=10.0)
        # loop through detections related to this template to save figures
        for k in range(len(family.detections)):
            # add one to detection counter
            detection_counter[j] = detection_counter[j] + 1
            # filter stream and save detection figure
            try:
                stream = streams[family.detections[k].id]
                stream.filter('bandpass', freqmin=1, freqmax=10)
                # plot streams and save fig
                outfile = save_dir + 'temp' + str(j+1) + '_det' + str(int(detection_counter[j])) + '.png'
                fig = stream.plot(equal_scale=False, size=(800, 600), outfile=outfile)
            except NotImplementedError:
                continue
        # empty family and streams to free memory
        family = None
        streams = None
    # before terminating the loop, save party and stm then empty to free memory
    party_file = save_dir + 'timestep' + str(i+1) + '_party.pickle'
    party_out = open(party_file,'wb')
    pickle.dump(party,party_out)
    pickle.dump(st,party_out)
    party_out.close()
    party = None
    st = None
print('Time Elapsed: ',UTCDateTime()-run_time)

# WINDOWED SCAN
# Determine plus minus window and initialize party
plus_minus_window = 7*24*60*60
party_dbd = Party()
# Loop through templates to conduct windowed scan for each of them
print('Conducting windowed scan for +/-',str(int(plus_minus_window/86400)),'days\n')
for template in tribe:
    # Extract template time and implement window
    template_time = template.event.origins[0].time
    t1 = template_time - plus_minus_window
    t2 = template_time + plus_minus_window
    # Detect events
    print('Scanning template ',str(template_time),'...')
    tribe_single = Tribe(template)
    party = tribe_single.client_detect(
        client=client, starttime=t1, endtime=t2, threshold=20.,
        threshold_type="MAD", trig_int=30.0, plot=False, return_stream=False)
    # Check if length of tribe matches number of families
    if len(tribe_single) != len(party.families):
        print('WARNING: Length of tribe does not match number of families. Template: ',str(template_time))
    # Append party of detections
    party_dbd = party_dbd + party
    party = None
# save windowed scan party
save_dir = '/home/ptan/project/output/'
tribe_file = save_dir + 'windowed_party.pickle'
tribe_out = open(tribe_file,'wb')
pickle.dump(tribe,tribe_out)
tribe_out.close()
