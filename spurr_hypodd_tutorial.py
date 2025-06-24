### HypoDD Tutorial for Spurr Seismicity (December 2024 to January 2025)
import pandas as pd
from toolbox import reader
from obspy import UTCDateTime
from functions import initialize_run, download_data, generate_dtcc, run_hypoDD, plot_hypoDD_results

### Initialize run (creates sub-directories like "data" and "output")
subdir_name = 'spurr'
initialize_run(subdir_name)

### Download data to local directory
# download_data(data_destination='./data/hypoDD_tutorial/',
#               starttime=UTCDateTime(2024,12,1),
#               endtime=UTCDateTime(2025,2,1),
#               client='IRIS',
#               network='AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV,AV',
#               station='RDT,SPBG,SPBG,SPBG,SPBL,SPBL,SPBL,SPCG,SPCG,SPCG,SPCL,SPCL,SPCL,SPCN,SPCN,SPCN,SPCP,SPCP,SPCP,SPNN,SPNN,SPNN,SPU,SPU,SPU,SPWE,SPWE,SPWE,STLK,STLK,STLK',
#               channel='BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ,BHE,BHN,BHZ',
#               location='--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--,--')

### Generate dt.cc file
relocatable_catalog = reader('./spurr_20241201_20250201_relocatable_catalog.xml')
relocate_catalog_output_dir = './output/spurr/relocate_catalog/'
data_path = './data/spurr/'
pre_pick_actual = 0.15
pre_pick_excess = 1
length_actual = 2.66
length_excess = 4
shift_len = 0.4
resampling_frequency = 100
min_link = 3
min_cc = 0.4
fmin = 1
fmax = 10

# generate_dtcc(catalog=relocatable_catalog,
#               relocate_catalog_output_dir=relocate_catalog_output_dir,
#               data_path=data_path,
#               pre_pick_actual=pre_pick_actual,
#               pre_pick_excess=pre_pick_excess,
#               length_actual=length_actual,
#               length_excess=length_excess,
#               shift_len=shift_len,
#               lowcut=fmin,
#               highcut=fmax,
#               resampling_frequency=resampling_frequency,
#               min_cc=min_cc,
#               min_link=min_link)

## (8) Relocate catalog candidates
D = 160
ph2dt_inc_dict = {'MEV': 12000,
                  'MSTA': 2600,
                  'MOBS': 12000}
ph2dt_inp_dict = {'MINWGHT': 0,
                  'MAXDIST': 200,
                  'MAXSEP': 10,
                  'MAXNGH': 10,
                  'MINLNK': 3,
                  'MINOBS': 3,
                  'MAXOBS': 50}
hypoDD_inc_dict = {'MAXEVE': 12000,
                   'MAXDATA': 5000000,
                   'MAXEVE0': 2,
                   'MAXDATA0': 1,
                   'MAXLAY': 30,
                   'MAXSTA': 100,
                   'MAXCL': 200}
hypoDD_inp_dict = {'IDAT': 3,
                   'IPHA': 3,
                   'DIST': 3000,
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
                   'DAMP': [D, D, D, D, D]}


stations = 'RDT,SPBG,SPBL,SPCG,SPCL,SPCN,SPCP,SPNN,SPU,SPWE,STLK'
avosta_df = pd.read_csv('./avo_stations.csv')
stalats = ','.join([str(avosta_df.Latitude[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
stalons = ','.join([str(avosta_df.Longitude[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
staelevs = ','.join([str(avosta_df.Elevation[list(avosta_df.Station).index(s)]) for s in stations.split(',')])
vzmodel_path = './vz/spurr_vzmodel.txt' # './vz/spurr_vzmodel2.txt'
has_ps_ratio = True  # If the input vz model has p/s ratios in the 3rd column
correct_depths = True

# ph2dt is bugged
hypoDD_loc, hypoDD_reloc = run_hypoDD(catalog=relocatable_catalog,
                                      relocate_catalog_output_dir=relocate_catalog_output_dir,
                                      stations=stations,
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

## (9) Plot relocated earthquakes
hypoDD_loc = reader(relocate_catalog_output_dir + 'hypoDD_loc.xml')
hypoDD_reloc = reader(relocate_catalog_output_dir + 'hypoDD_reloc.xml')
lat_lims = [61.25, 61.375]
lon_lims = [-152.5, -152.15]
dep_lims = [-5, 10]
markersize = 2
figsize = (10,7)

plot_hypoDD_results(hypoDD_in=relocatable_catalog,
                    hypoDD_loc=hypoDD_loc,
                    hypoDD_reloc=hypoDD_reloc,
                    lat_lims=lat_lims,
                    lon_lims=lon_lims,
                    dep_lims=dep_lims,
                    markersize=markersize,
                    figsize=figsize,
                    legend_loc='lower left')
