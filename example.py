## Example script that executes the following in sequence:
# 1. Download data
# 2. Run REDPy
# 3. Convert REDPy to EQcorrscan-compatible outputs
# 4. Create EQcorrscan templates
# 5. Run EQcorrscan matched-filter
# 6. Generate dt.cc file
# 7. Use dt.ct relocation on HypoDD
# 8a. Use dt.cc relocation on HypoDD
# 8b. Use dt.cc relocation on GrowClust (commented out)
# 9. Plot resulting earthquakes in time and space

# Import all dependencies
import os
from obspy import UTCDateTime
from functions import download_data, run_redpy

# (0) Prepare output directory and parse AVO catalog
output_dir = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/example/'
print('Creating subdirectories for workflow outputs...')
output_subdirs = [output_dir, output_dir + 'run_redpy', output_dir + 'convert_redpy/', output_dir + 'create_tribe/',
                  output_dir + 'scan_data/', output_dir + 'relocate_catalog/']
for output_subdir in output_subdirs:
    try:
        os.mkdir(output_subdir)
    except FileExistsError:
        print('%s already exists.' % output_subdir)
print('All output subdirectories created.')

# (1) Download data
data_destination = '/Users/darrentpk/Desktop/example_data/'
starttime = UTCDateTime(2009,2,25,0,0,0)
endtime = UTCDateTime(2009,2,28,0,0,0)
client = 'IRIS'
network = 'AV,AV,AV'
station = 'REF,RDN,RSO'
channel = 'EHZ,EHZ,EHZ'
location = '--,--,--'

# download_data(data_destination=data_destination,
#               starttime=starttime,
#               endtime=endtime,
#               client=client,
#               network=network,
#               station=station,
#               channel=channel,
#               location=location)

# (2) Run REDPy
run_title = 'Redoubt Example'
output_destination = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/output/example/run_redpy/'
data_path = data_destination
stalats=60.5224,60.4888,60.4616  # for teleseism removal
stalons=-152.7401,-152.694,-152.756  # for teleseism removal
samprate=100  # sampling rate for stations -- non-matching traces will be resampled
fmin=1.  # bandpass filter low bound
fmax=10.  # bandpass filter high bound
nstaC=2  # number of stations needed to coincidentally trigger
lwin=8  # STA/LTA long window
swin=0.7  # STA/LTA short window
trigon=3  # STA/LTA trigger on ratio
trigoff=2  # STA/LTA trigger off ratio
winlen=1024  # cross-correlation window length in samples
cmin=0.85  # minimum cross-correlation coefficient value to consider a repeater
ncor=2  # minimum number of stations where cmin must be met to determine a repeater
minorph=0.05  # amount of days to keep orphans in the queue when it just triggers above threshold (> trigon)
maxorph=7  # amount of days to keep orphans in the queue when it triggers way above threshold (> trigon+7)

run_redpy(run_title=run_title,
          output_destination=output_destination,
          data_path=data_path,
          starttime=starttime,
          endtime=endtime,
          network=network,
          station=station,
          channel=channel,
          location=location,
          stalats=stalats,
          stalons=stalons,
          samprate=samprate,
          fmin=fmin,
          fmax=fmax,
          nstaC=nstaC,
          lwin=lwin,
          swin=swin,
          trigon=trigon,
          trigoff=trigoff,
          winlen=winlen,
          cmin=cmin,
          ncor=ncor,
          minorph=minorph,
          maxorph=maxorph)