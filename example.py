## Example script that executes the following in sequence:
# 1. Download data
# 2. Run REDPy
# 3. Convert REDPy into EQcorrscan-compatible ObsPy catalogs
# 4. Create EQcorrscan templates
# 5. Run EQcorrscan matched-filter scan
# 6. Calculate frequency index and relative magnitudes
# 7. Generate dt.cc file for relocatable catlaog
# 8. Relocate catalog candidates
# 9. Plot relocated earthquakes

# Import all dependencies
from obspy import UTCDateTime
from functions import initialize_run, download_data, run_redpy, convert_redpy, create_tribe, scan_data
from functions import rethreshold_results, generate_dtcc, run_hypoDD, plot_hypoDD_results
from toolbox import reader, writer, calculate_catalog_FI, calculate_relative_magnitudes

## (0) Prepare output directory and parse AVO catalog
subdir_name = 'example'
catalog = reader('./redoubt_20080225_20090227.xml')  # Filtered AVO catalog for the 3 day example
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
dybin = 1/24  # length in days for each bin in the histogram subplot (use 1 hr as our analysis only spans 3 days)

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
          maxorph=maxorph,
          dybin=dybin)

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

## (7) Generate dt.cc file for relocatable catalog
relocate_catalog_output_dir = './output/' + subdir_name + '/relocate_catalog/'
pre_pick_actual = 0.15
pre_pick_excess = 1
length_actual = 2.66
length_excess = 4
shift_len = 0.4
min_link = 3
min_cc = 0.7

generate_dtcc(catalog=relocatable_catalog_FImag,
              relocate_catalog_output_dir=relocate_catalog_output_dir,
              data_path=data_path,
              pre_pick_actual=pre_pick_actual,
              pre_pick_excess=pre_pick_excess,
              length_actual=length_actual,
              length_excess=length_excess,
              shift_len=shift_len,
              lowcut=fmin,
              highcut=fmax,
              min_cc=min_cc,
              min_link=min_link)

## (8) Relocate catalog candidates
ph2dt_inc_dict = {'MEV': 12000,
                  'MSTA': 2600,
                  'MOBS': 1000}
ph2dt_inp_dict = {'MINWGHT': 0,
                  'MAXDIST': 120,
                  'MAXSEP': 10,
                  'MAXNGH': 10,
                  'MINLNK': 3,  # 8
                  'MINOBS': 3,  # 8
                  'MAXOBS': 20}
hypoDD_inc_dict = {'MAXEVE': 1000,
                   'MAXDATA': 20000,
                   'MAXEVE0': 20,
                   'MAXDATA0': 2000,
                   'MAXLAY': 30,
                   'MAXSTA': 100,
                   'MAXCL': 200}
hypoDD_inp_dict = {'IDAT': 3,
                   'IPHA': 3,
                   'DIST': 50,
                   'OBSCC': 0,
                   'OBSCT': 0,
                   'MINDS': -999,
                   'MAXDS': -999,
                   'MAXGAP': -999,
                   'ISTART': 2,
                   'ISOLVE': 2,
                   'IAQ': 0,
                   'NSET': 5,
                   'NITER': [10, 10, 10, 10, 10],
                   'WTCCP': [0.01, 0.01, 0.10, 0.50, 1.00],
                   'WTCCS': [0.01, 0.01, 0.10, 0.50, 1.00],
                   'WRCC': [10, 10, 10, 6, 6],
                   'WDCC': [4, 4, 4, 2, 2],
                   'WTCTP': [1.00, 1.00, 1.00, 0.10, 0.01],
                   'WTCTS': [1.00, 1.00, 1.00, 0.10, 0.01],
                   'WRCT': [12, 6, 6, 6, 6],
                   'WDCT': [10, 5, 2.5, 2.5, 2.5],
                   'DAMP': [100, 100, 100, 100, 100]}

relocatable_catalog_FImag = reader(scan_data_output_dir + 'relocatable_catalog_FImag.xml')
stalats = '60.5913,60.5621,60.5905,60.5224,60.5726,60.4875,60.4888,60.4616'
stalons = '-152.6882,-152.9293,-152.8058,-152.7401,-152.4075,-152.8425,-152.6940,-152.7560'
staelevs = '1090,1136,1414,1400,930,1546,1641,1921'
vzmodel_path = './vz/redoubt_vzmodel.txt'
has_ps_ratio = True  # If the input vz model has ratios in the 3rd column, set to None
correct_depths = True

hypoDD_loc, hypoDD_reloc = run_hypoDD(catalog=relocatable_catalog_FImag,
                                      relocate_catalog_output_dir=relocate_catalog_output_dir,
                                      stations=station,
                                      stalats=stalats,
                                      stalons=stalons,
                                      staelevs=staelevs,
                                      vzmodel_path=vzmodel_path,
                                      has_ps_ratio=has_ps_ratio,
                                      correct_depths=correct_depths,
                                      ph2dt_inc_dict=ph2dt_inc_dict,
                                      ph2dt_inp_dict=ph2dt_inp_dict,
                                      hypoDD_inc_dict=hypoDD_inc_dict,
                                      hypoDD_inp_dict=hypoDD_inp_dict)

## 9. Plot relocated earthquakes
relocatable_catalog_FImag = reader(scan_data_output_dir + 'relocatable_catalog_FImag.xml')
hypoDD_loc = reader(relocate_catalog_output_dir + 'hypoDD_loc.xml')
hypoDD_reloc = reader(relocate_catalog_output_dir + 'hypoDD_reloc.xml')
lat_lims = [60.4, 60.6]
lon_lims = [-152.9, -152.6]
dep_lims = [-3, 10]

plot_hypoDD_results(hypoDD_in=relocatable_catalog_FImag,
                    hypoDD_loc=hypoDD_loc,
                    hypoDD_reloc=hypoDD_reloc,
                    lat_lims=lat_lims,
                    lon_lims=lon_lims,
                    dep_lims=dep_lims)