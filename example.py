## Example script that executes the following in sequence:
# 1. Download data
# 2. Run REDPy
# 3. Convert REDPy into EQcorrscan-compatible ObsPy catalogs
# 4. Create EQcorrscan templates
# 5. Run EQcorrscan matched-filter scan
# 6. Calculate frequency index and relative magnitudes
# 7. Generate dt.cc file
# 8. Relocate catalog candidates
# 9. Plot resulting earthquakes in time and space

# Import all dependencies
from obspy import UTCDateTime
from functions import initialize_run, download_data, run_redpy, convert_redpy, create_tribe, scan_data, rethreshold_results
from toolbox import reader, writer, calculate_catalog_FI, calculate_relative_magnitudes

## (0) Prepare output directory and parse AVO catalog
subdir_name = 'example'
catalog = reader('./redoubt_20080401_20090901.xml')  # Filtered AVO catalog for the 3 day example
initialize_run(subdir_name)

## (1) Download data
data_destination = './data/'+ subdir_name + '/'
starttime = UTCDateTime(2009,2,25,0,0,0)
endtime = UTCDateTime(2009,2,28,0,0,0)
client = 'IRIS'
network = 'AV,AV,AV,AV,AV,AV,AV,AV'
station = 'DFR,NCT,RDJH,RDN,RDT,RDWB,REF,RSO'
channel = 'EHZ,EHZ,BHZ,EHZ,EHZ,BHZ,EHZ,EHZ'
location = '--,--,--,--,--,--,--,--'

download_data(data_destination=data_destination,
              starttime=starttime,
              endtime=endtime,
              client=client,
              network=network,
              station=station,
              channel=channel,
              location=location)

## (2) Run REDPy
run_title = 'Redoubt Example'
redpy_output_destination = './output/' + subdir_name + '/run_redpy/'
data_path = data_destination
redpy_network = 'AV,AV,AV'
redpy_station = 'RDN,REF,RSO'
redpy_channel = 'EHZ,EHZ,EHZ'
redpy_location = '--,--,--'
redpy_stalats = 60.5224,60.4888,60.4616  # for teleseism removal
redpy_stalons = -152.7401,-152.694,-152.756  # for teleseism removal
samprate = 100  # sampling rate for stations -- non-matching traces will be resampled
fmin = 1.  # bandpass filter low bound
fmax = 10.  # bandpass filter high bound
nstaC = 2  # number of stations needed to coincidentally trigger
lwin = 8  # STA/LTA long window
swin = 0.7  # STA/LTA short window
trigon = 3  # STA/LTA trigger on ratio
trigoff = 2  # STA/LTA trigger off ratio
winlen = 1024  # cross-correlation window length in samples
cmin = 0.85  # minimum cross-correlation coefficient value to consider a repeater
ncor = 2  # minimum number of stations where cmin must be met to determine a repeater
minorph = 0.05  # amount of days to keep orphans in the queue when it just triggers above threshold (> trigon)
maxorph = 7  # amount of days to keep orphans in the queue when it triggers way above threshold (> trigon+7)

run_redpy(run_title=run_title,
          output_destination=redpy_output_destination,
          data_path=data_path,
          starttime=starttime,
          endtime=endtime,
          network=redpy_network,
          station=redpy_station,
          channel=redpy_channel,
          location=redpy_location,
          stalats=redpy_stalats,
          stalons=redpy_stalons,
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

## (3) Convert REDPy into separate catalog objects
analyst_catalog = catalog
redpy_output_path = redpy_output_destination + "".join(run_title.split()) + '/'
convert_redpy_output_dir = './output/' + subdir_name + '/convert_redpy/'
add_redpy_pick_to_associated = True  # Add STA/LTA pick to associated catalog events
add_campaign_pick_to_associated = False  # Campaign data not downloaded for test case

convert_redpy(analyst_catalog=analyst_catalog,
              redpy_output_path=redpy_output_path,
              convert_redpy_output_dir=convert_redpy_output_dir,
              data_path=data_path,
              redpy_station=redpy_station,
              redpy_channel=redpy_channel,
              fmin=fmin,
              fmax=fmax,
              nstaC=nstaC,
              lwin=lwin,
              swin=swin,
              trigon=trigon,
              trigoff=trigoff,
              add_redpy_pick_to_associated=add_redpy_pick_to_associated,
              add_campaign_pick_to_associated=add_campaign_pick_to_associated)

## (4) Convert REDPy into EQcorrscan-compatible ObsPy catalogs
convert_redpy_output_dir = convert_redpy_output_dir
create_tribe_output_dir = './output/' + subdir_name + '/create_tribe/'
samprate = 50  # desired sampling rate for templates
prepick = 1  # time before pick time to start template waveform trim (s)
length = 8  # time from pre-pick to stop template waveform trim (s)
min_snr = 1  # minimum signal-to-noise to accept waveform into template

tribe = create_tribe(convert_redpy_output_dir=convert_redpy_output_dir,
                     create_tribe_output_dir=create_tribe_output_dir,
                     data_path=data_path,
                     template_stations=station,
                     resampling_frequency=samprate,
                     lowcut=fmin,
                     highcut=fmax,
                     prepick=prepick,
                     length=length,
                     min_snr=min_snr)

## (5) Run EQcorrscan matched-filter scan
# tribe = reader(create_tribe_output_dir + 'tribe.tgz')
scan_data_output_dir = './output/' + subdir_name + '/scan_data/'
min_stations = 3  # minimum number of stations where each template must be observed on before being used for matched-filter
min_picks = 0  # minimum number of picks that each template must possess before being used for matched-filter
samprate = 50  # standardized sampling rate used for seismic data matched-filter scan (make sure this is equal to template samprate)
threshold_type = 'av_chan_corr'  # EQcorrscan threshold type -- choose between 'MAD', 'absolute', 'av_chan_corr'
threshold = 0.7  # threshold value used for matched-filter detections
trig_int = 8  # minimum trigger interval for individual template (s)
decluster = True  # remove overlapping detections from different templates

party, detected_catalog, relocatable_catalog = scan_data(tribe=tribe,
                                                         scan_data_output_dir=scan_data_output_dir,
                                                         data_path=data_path,
                                                         min_stations=min_stations,
                                                         min_picks=min_picks,
                                                         starttime=starttime,
                                                         endtime=endtime,
                                                         resampling_frequency=samprate,
                                                         threshold_type=threshold_type,
                                                         threshold=threshold,
                                                         trig_int=trig_int,
                                                         decluster=decluster)

# # Optional rethresholding option after conducting a manual sensitivity test
# # Note that it is strongly recommended to run scan_data with decluster=False if rethresholding is desired
# new_threshold = 0.75
# decluster = True
# party, detected_catalog, relocatable_catalog = rethreshold_results(tribe=tribe,
#                                                                    party=party,
#                                                                    threshold_type=threshold_type,
#                                                                    new_threshold=new_threshold,
#                                                                    decluster=decluster,
#                                                                    trig_int=trig_int)
# writer(scan_data_output_dir + 'party.tgz', party)
# writer(scan_data_output_dir + 'detected_catalog.xml', detected_catalog)
# writer(scan_data_output_dir + 'relocatable_catalog.xml', relocatable_catalog)

## (6) Calculate frequency index and relative magnitudes
# party = reader(scan_data_output_dir + 'party.tgz')
# detected_catalog = reader(scan_data_output_dir + 'detected_catalog.xml')

# Settings for FI calculation
reference_station = redpy_station  # stations used to compute the average spectra
reference_channel = redpy_channel  # corresponding channels used to compute the average spectra
min_match = 2  # minimum number of matches with reference channels before FI calculation

detected_catalog_FI = calculate_catalog_FI(catalog=detected_catalog,
                                           data_path=data_path,
                                           reference_station=reference_station,
                                           reference_channel=reference_channel,
                                           min_match=min_match,
                                           resampling_frequency=samprate,
                                           prepick=prepick,
                                           length=length,
                                           lowcut=fmin,
                                           highcut=fmax)

relocatable_catalog_FI = calculate_catalog_FI(catalog=relocatable_catalog,
                                              data_path=data_path,
                                              reference_station=reference_station,
                                              reference_channel=reference_channel,
                                              min_match=min_match,
                                              resampling_frequency=samprate,
                                              prepick=prepick,
                                              length=length,
                                              lowcut=fmin,
                                              highcut=fmax)

# # Settings for magnitude calculation
noise_window = (-20.0, -prepick)  # time window to calculate noise amplitude (s, s)
signal_window = (-prepick, -prepick+length)  # time window to calculate signal amplitude (s, s)
shift_len = 1.5  # length of time shift to find max cc between each event pair (s)
min_cc = 0.7  # minimum cross-correlation coefficient to calculate relative magnitude
min_snr = 1  # minimum signal-to-noise to calculate relative magnitude

relocatable_catalog_FImag = calculate_relative_magnitudes(catalog=relocatable_catalog_FI,
                                                          tribe=tribe,
                                                          data_path=data_path,
                                                          noise_window=noise_window,
                                                          signal_window=signal_window,
                                                          min_cc=min_cc,
                                                          min_snr=min_snr,
                                                          shift_len=shift_len,
                                                          resampling_frequency=samprate,
                                                          lowcut=fmin,
                                                          highcut=fmax)

# Write out new catalogs
writer(scan_data_output_dir + 'detected_catalog_FI.xml', detected_catalog_FI)
writer(scan_data_output_dir + 'relocatable_catalog_FImag.xml', relocatable_catalog_FImag)