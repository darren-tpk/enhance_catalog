# Import all dependencies
import os
import tarfile
import subprocess
import pandas as pd
from obspy import Catalog, UTCDateTime
from obspy.core.event import Event, Origin, Comment, Pick, WaveformStreamID, Arrival, ResourceIdentifier, Magnitude
from toolbox import reader, writer
import matplotlib.pyplot as plt
from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper

# Define directories and filenames
main_dir = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/'
output_dir = main_dir + 'output/mammoth3/relocate_catalog/'
hypoddver_dir = main_dir + 'hypoDDver/'
hypodd_tarfile = 'HYPODD_2.1b.tar.gz'
hypodd_dir = hypoddver_dir + 'HYPODD/'
include_dir = hypodd_dir + 'include/'
src_dir = hypodd_dir + 'src/'

# Choose operations
xml_filepath = hypoddver_dir + 'IN/relocatable_catalog.xml'
raw_station_list_filepath = main_dir + 'data/stations/mammoth_station_list.csv'
# Edit problem setup
correct_depths = True  # This shifts the problem down based on the max height of the input velocity model
catalog_loaded = True
# Actual hypoDD operations
run_ncsn2pha = False
run_ph2dt = True
run_hypodd = True
# Make quick plots illustrating results
plot_results = True

# Define configurations
## ph2dt.inc
MEV = 12000  # Maximum number of events; recommended: 12000
MSTA = 2000  # Maximum number of stations; recommended: 2600
MOBS = 1000  # Maximum number of phases (P&S) per event; recommended: 1000
## hyppoDD.inc
MAXEVE = 6000  # Maximum number of events; recommended: large: 6500, small: 1000
MAXDATA = 7000000  # Maximum number of observations; recommended: large: 3000000, small: 20000
MAXEVE0 = 2  # Maximum number of events used for SVD; recommended: large: 215, small: 20, LSQR used: 2
MAXDATA0 = 1  # Maximum number of observations used for SVD; recommended: large: 100000, small: 2000, LSQR used: 1
MAXLAY = 50  # Maximum number of model layers; recommended: large: 50, small: 30
MAXSTA = 50  # Maximum number of stations; recommended: large: 400, small: 100
MAXCL = 100  # Maximum number of clusters allowed; recommended: large: 200, small: 200
## ncsn2pha inputs
arc_file = None  # Not running ncsn2pha for workflow
phase_file = './hypoDDver/IN/phase.dat'
## ph2dt.inp
stlist_file = './hypoDDver/IN/stations.dat'
phase_file = './hypoDDver/IN/phase.dat'
MINWGHT = 0  # Minimum pick weight; out-of-box: 0
MAXDIST = 200  # Maximum distance between event pair and station; out-of-box: 120
MAXSEP = 30  # Maximum hypocentral separation between event pairs; out-of-box: 10
MAXNGH = 10  # Maximum number of neighbors per event; out-of-box: 10
MINLNK = 8  # Minimum number of links required to define a neighbor; out-of-box: 8
MINOBS = 8  # Minimum number of links per pair saved; out-of-box: 1
MAXOBS = 100  # Maximum number of links per pair saved; out-of-box: 100
## hypoDD.inp
#### IN
dtcc_file = './hypoDDver/IN/dt.cc'
dtct_file = './hypoDDver/OUT_ph2dt/dt.ct'
evlist_file = './hypoDDver/OUT_ph2dt/event_hypodd.dat' # './hypoDDver/IN/detection.dat' #
stlist_file = './hypoDDver/IN/stations.dat'
vzmodel_file = './data/vz/mammoth_vzmodel3.txt'  # NOT parsed into hypodd directly
ps_ratio = None  # If the input vz model has ratios in the 3rd column, set to None
#### OUT
loc_file = './hypoDDver/OUT_hypodd/hypoDD.loc'
reloc_file = './hypoDDver/OUT_hypodd/hypoDD.reloc'
sta_file = './hypoDDver/OUT_hypodd/hypoDD.sta'
res_file = './hypoDDver/OUT_hypodd/hypoDD.res'
srcparam_file = './hypoDDver/OUT_hypodd/hypoDD.src'
#### config
###### data type selection
IDAT = 3  # Data Type; 1: cc only, 2: ct only, 3: cc & ct
IPHA = 3  # Phase Type; 1: P only, 2: S only, 3: P & S
DIST = 50  # Maximum distance between cluster centroid and stations
###### event clustering
OBSCC = 0  # Minimum number of cc links for cluster, out-of-box: 0
OBSCT = 0  # Minimum number of ct links for cluster, out-of-box: 8
MINDS = -999  # Min distance between individual event pairs and stations
MAXDS = -999  # Max distance between individual event pairs and stations
MAXGAP = -999  # Max azimuthal gap between individual event pairs and stations. (-999 or 365: not used)
###### solution control
ISTART = 2  # Initial locations; 1: start from centroid, 2: start from catalog locations
ISOLV = 2  # Least squares solution: 1: SVD, 2: LSQR
IAQ = 0  # Keep (0) or Remove (1) air quakes
NSET = 5  # Number of sets of iteration specifications following (if NSET=2, next section should be 2 element lists)
###### data weight and re-weighting
NITER = [10, 10, 10, 10, 10] # Number of iterations for the set of weighting parameters that follow
WTCCP = [0.01, 0.01, 0.10, 0.50, 1.00]  #[0.001, 0.001, 1, 1] #  Weight for cc P wave (-9: omit)
WTCCS = [0.01, 0.01, 0.10, 0.50, 1.00] #[0.001, 0.001, 0.8, 0.8] # Weight for cc S wave (-9: omit)
WRCC = [10, 10, 10, 6, 6] #[10, 10, 10, 6] # Cutoff threshold for outliers located on cc tails (-9: no outlier removed)
WDCC = [4, 4, 4, 2, 2] #[4, 4, 4, 0.5] # Maximum event separation for cc data (-9: not activated)
WTCTP = [1.00, 1.00, 1.00, 0.10, 0.01] #[1, 1, 0.001, 0.001] # Weight for ct P wave (-9: omit)
WTCTS = [1.00, 1.00, 1.00, 0.10, 0.01] #[0.8, 0.8, 0.001, 0.001] # Weight for ct S wave (-9: omit)
WRCT = [12, 6, 6, 6, 6] #[12, 6, 6, 6] # Cutoff threshold for outliers located on ct tails (-9: no outlier removed)
WDCT = [10, 5, 2.5, 2.5, 2.5] #[8, 3, 3, 3] # Maximum event separation for ct data (-9: not activated)
DAMP = [50, 50, 50, 50, 50] #[100, 100, 100, 100] # Damping (only for LSQR, ISOLV=2)
###### 1D velocity model
# Instead of parsing the velocity model via NLAY, RATIO, TOP and VEL, we process a txt file input accordingly (see vzmodel_file)
###### event selection
CID = 0
ID = '*'  # If relocating all events, input '*'


#################
# Define functions that we need

# Converts .loc and .reloc into Catalog objects
def loc2cat(loc_filepath, event_id_mapper=None, input_catalog=None, type='loc', depth_correction=0):
    # Get list of input event resource ids
    if event_id_mapper is not None and input_catalog is not None:
        all_ids = [event.resource_id.id for event in input_catalog]
    # Initialize output catalog
    outcat = Catalog()
    # Open loc file and read lines
    with open(loc_filepath, "r") as loc:
        lines = loc.read()
        lines = lines.split('\n')
        # Loop over lines
        for line in lines[:-1]:
            # Extract information and craft event
            if type == 'loc':
                (evid, lat, lon, dep, _, _, _, _, _, _, yr, mo, dy, hr, mm, ss, mag, cid) = line.split()
            elif type == 'reloc':
                (evid, lat, lon, dep, _, _, _, _, _, _, yr, mo, dy, hr, mm, ss, mag, _, _, _, _, _, _,
                 cid) = line.split()
                if lat == 'NaN':
                    continue
            else:
                raise ValueError('type argument is not loc/reloc!')
            yr = int(yr)
            mo = int(mo)
            dy = int(dy)
            hr = int(hr)
            mm = int(hr)
            ss = float(ss)
            if ss < 60:
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            else:
                ss -= 60
                mm += 1
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            lon = float(lon)
            lat = float(lat)
            dep = (float(dep) - depth_correction)
            mag = float(mag)
            ev = Event(origins=[Origin(time=time, longitude=lon, latitude=lat, depth=dep)],
                       magnitudes=[Magnitude(mag=mag)])
            # Find base event and copy over comments
            old_event_id = [id for id, num in event_id_mapper.items() if num == int(evid)]
            if len(old_event_id) != 1:
                raise ValueError('Multiple matched on event id mapper. Check event ids!')
            else:
                old_event_id = old_event_id[0]
                old_event = input_catalog[all_ids.index(old_event_id)]
                ev.comments = old_event.comments
                evid_origin_comment = 'Event ID: %s' % evid
                ev.origins[0].comments.append(Comment(text=evid_origin_comment))
            # Append event to catalog
            outcat.append(ev)
    # Return catalog
    return outcat


# Converts .dat and .sel to Catalog objects
def dat2cat(dat_filepath, event_id_mapper=None, input_catalog=None, depth_correction=0):
    # Get list of input event resource ids
    if event_id_mapper is not None and input_catalog is not None:
        all_ids = [event.resource_id.id for event in input_catalog]
    # Initialize output catalog
    outcat = Catalog()
    # Open loc file and read lines
    with open(dat_filepath, "r") as dat:
        lines = dat.read()
        lines = lines.split('\n')
        # Loop over lines
        for line in lines[:-1]:
            # Extract information and craft event
            (yrmody, hrmmsscs, lat, lon, dep, mag, _, _, _, evid) = line.split()
            yr = int(yrmody[0:4])
            mo = int(yrmody[4:6])
            dy = int(yrmody[6:8])
            hrmmsscs = '%08d' % int(hrmmsscs)
            hr = int(hrmmsscs[:-6])
            mm = int(hrmmsscs[-6:-4])
            ss = float(hrmmsscs[-4:-2]) + float(hrmmsscs[-2:]) / 100
            if ss < 60:
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            else:
                ss -= 60
                mm += 1
                time = UTCDateTime(yr, mo, dy, hr, mm, ss)
            ev = Event(origins=[Origin(time=time, longitude=float(lon), latitude=float(lat),
                                       depth=(float(dep) - depth_correction))],
                       magnitudes=[Magnitude(mag=float(mag))])
            # Find base event and copy over resource id and comments
            old_event_id = [id for id, num in event_id_mapper.items() if num == int(evid)]
            if len(old_event_id) != 1:
                raise ValueError('Multiple matched on event id mapper. Check event ids!')
            else:
                old_event_id = old_event_id[0]
                old_event = input_catalog[all_ids.index(old_event_id)]
                ev.comments = old_event.comments
                evid_origin_comment = 'Event ID: %s' % evid
                ev.origins[0].comments.append(Comment(text=evid_origin_comment))
            # Append event to catalog
            outcat.append(ev)
    # Return catalog
    return outcat


# Returns a copied catalogA, but repeat events in catalogB are removed
def remove_catalog_repeats(catalogA, catalogB):
    catalogA_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogA]
    catalogB_evids = [int(event.origins[0].comments[0].text.split()[-1]) for event in catalogB]
    catalogC = catalogA.copy()
    for i, catalogA_evid in reversed(list(enumerate(catalogA_evids))):
        if catalogA_evid in catalogB_evids:
            catalogC.events.pop(i)
    return catalogC


#################

# Create directories for intermediary and final outputs
print('Creating sub directories...')
for dir in [hypoddver_dir + 'INP', hypoddver_dir + 'OUT_ph2dt', hypoddver_dir + 'OUT_hypodd']:
    try:
        os.mkdir(dir)
    except FileExistsError:
        print('Subdirectory %s already exists.' % dir)

# If depth_correction = True, we want to shift the entire problem downwards
if correct_depths:
    print('WARNING: correcting elevations and depths for input catalog, vzmodel and stations')
    # Read the input vzmodel and use the
    vzmodel_read = open(vzmodel_file, "r")
    vzmodel_lines = vzmodel_read.readlines()
    model_ceiling = float(vzmodel_lines[0].split()[0])
    depth_correction = -1 * model_ceiling
    vzmodel_read.close()
else:
    depth_correction = 0

# Deconstruct .xml file

print('Deconstructing .xml file into phase.dat and stations.dat ...')

# Output filepaths
used_station_list_filepath = hypoddver_dir + 'IN/stations.dat'
phase_list_filepath = hypoddver_dir + 'IN/phase.dat'
detection_list_filepath = hypoddver_dir + 'IN/detection.dat'

# Read catalog and raw station list
if not catalog_loaded:
    catalog = reader(xml_filepath)
raw_station_list = pd.read_csv(raw_station_list_filepath)

# Loop through catalog event picks to get a unique list of stations used
stations_used = []
for event in catalog:
    for pick in event.picks:
        stations_used.append(pick.waveform_id.station_code)
stations_used = list(set(stations_used))

# Write out station list input
used_station_list_file = open(used_station_list_filepath, 'w')
used_station_list_format = '%s %.4f %.4f %.1f\n'
for station_used in stations_used:
    index = list(raw_station_list.station).index(station_used)
    station_lat = raw_station_list.latitude[index]
    station_lon = raw_station_list.longitude[index]
    if correct_depths:
        station_elev = raw_station_list.elevation[index] - depth_correction * 1000
    else:
        station_elev = raw_station_list.elevation[index]
    used_station_list_file.write(used_station_list_format % (station_used, station_lat, station_lon, station_elev))
used_station_list_file.close()

# Loop through catalog and picks again to populate event and phase text files

# Generate event id mapper
event_id_mapper = _generate_event_id_mapper(catalog, event_id_mapper=None)
# Initialize phase list and define line formats
phase_list_file = open(phase_list_filepath, 'w')
phase_list_format1 = '#  %4d %2d %2d %2d %2d %2d.%02d %8.4f %9.4f %8.3f %5.2f %.2f %.2f %.2f %10d\n'
phase_list_format2 = '%s %.3f %.2f  %2s\n'
# Initialize detection list and define line formats
detection_list_file = open(detection_list_filepath, 'w')
detection_list_format = '%4d%02d%02d  %2d%02d%02d%02d %9.4f %10.4f %10.3f %6.2f %7.2f %7.2f %6.2f %10d\n'
# Loop through events in catalog
for event in catalog:
    # Extract key information
    yr = int(event.origins[0].time.year)
    mo = int(event.origins[0].time.month)
    dy = int(event.origins[0].time.day)
    hh = int(event.origins[0].time.hour)
    mm = int(event.origins[0].time.minute)
    ss = int(event.origins[0].time.second)
    cs = int(event.origins[0].time.microsecond / 1e4)
    lat = event.origins[0].latitude
    lon = event.origins[0].longitude
    if correct_depths:
        dep = event.origins[0].depth + depth_correction
    else:
        dep = event.origins[0].depth
    if event.magnitudes == []:
        mag = -9
    else:
        mag_object = event.preferred_magnitude() or event.magnitudes[0]
        mag = mag_object.mag
    eh = 0  # Change accordingly if event object contains error info
    ez = 0  # Change accordingly if event object contains error info
    rms = 0  # Change accordingly if event object contains error info
    evid = event_id_mapper[event.resource_id.id]
    # Check event time against template to see if it is a self detection
    event_time = event.origins[0].time
    try:
        template_name = event.comments[0].text.split()[1]
        template_time_segments = [int(numstr) for numstr in template_name.replace('t', '_').split('_')]
        template_time = UTCDateTime(*template_time_segments)
        time_diff = abs(event_time - template_time)
    except:
        time_diff = 0
    # If it is a self detection, write info into phase file for ph2dt
    if time_diff < 1:
        phase_list_file.write(phase_list_format1 % (yr, mo, dy, hh, mm, ss, cs, lat, lon, dep, mag, eh, ez, rms, evid))
        for i, pick in enumerate(event.picks):
            sta = pick.waveform_id.station_code
            tt = pick.time - event_time
            if event.origins[0].arrivals == []:
                wght = 1.0  # Implement weighting scheme ?
            else:
                wght = event.origins[0].arrivals[i].time_weight
            pha = pick.phase_hint
            phase_list_file.write(phase_list_format2 % (sta, tt, wght, pha))
    else:
        # Otherwise, it is a detection, and we write it into our detection list file for appending later
        detection_list_file.write(
            detection_list_format % (yr, mo, dy, hh, mm, ss, cs, lat, lon, dep, mag, eh, ez, rms, evid))
# close files
phase_list_file.close()
detection_list_file.close()

# Prepare hypodd inputs

# Extract contents from hypoDD tarfile
print('Extracting hypoDD tarfile contents ...')
hypodd_tar = tarfile.open(hypoddver_dir + hypodd_tarfile, "r:gz")
hypodd_tar.extractall(hypoddver_dir)

if run_ph2dt:
    print('Configuring ph2dt.inc ...')
    # Configure ph2dt inc file (adapted from hypoDDpy, https://github.com/krischer/hypoDDpy)
    ph2dt_inc = """
c ph2dt.inc: Stores parameters that define array dimensions in ph2dt.
c            Modify to fit size of problem and available computer memory.
c Parameter Description:
c MEV:   Max number of events.
c MSTA:  Max number of stations.
c MOBS:  Max number of phases (P&S) per eventer event.

      integer	MEV, MSTA, MOBS

      parameter(MEV=    {MEV},
     &		MSTA=   {MSTA},
     &		MOBS=   {MOBS})""".format(
        MEV=MEV,
        MSTA=MSTA,
        MOBS=MOBS,
    )
    # Remove the leading empty line
    ph2dt_inc = ph2dt_inc[1:]

    # Rename default ph2dt.inc file and write out new ph2dt.inc file
    ph2dt_inc_rename_command = 'mv %s/ph2dt.inc %s/ph2dt_default.inc' % (include_dir, include_dir)
    subprocess.call(ph2dt_inc_rename_command, shell=True)
    with open(include_dir + 'ph2dt.inc', "w") as open_file:
        ph2dt_inc_file = open_file.write(ph2dt_inc)
    open_file.close()

if run_hypodd:
    print('Configuring hypoDD.inc ...')
    # Configure hypoDD inc file (adapted from hypoDDpy, https://github.com/krischer/hypoDDpy)
    hypodd_inc = """
c hypoDD.inc: Stores parameters that define array dimensions in hypoDD.
c             Modify to fit size of problem and available computer memory.
c	      If 3D raytracing is used, also set model parameters in vel3d.inc.
c Parameter Description:
c MAXEVE:   Max number of events (must be at least the size of the number 
c           of events listed in the event file)
c MAXDATA:  Max number of observations (must be at least the size of the 
c           number of observations).  
c MAXEVE0:  Max number of events used for SVD. If only LSQR is used, 
c           MAXEVE0 can be set to 2 to free up memory. 
c MAXDATA0: Max number of observations used for SVD. If only LSQR is used, 
c           MAXDATA0 can be set to 1 to free up memory. 
c MAXLAY:   Max number of model layers.
c MAXSTA:   Max number of stations.
c MAXCL:    Max number of clusters allowed. 
      integer*4 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0
      integer*4 MAXCL

      parameter(MAXEVE   = {MAXEVE},
     &          MAXDATA  = {MAXDATA},
     &          MAXEVE0  = {MAXEVE0},
     &          MAXDATA0 = {MAXDATA0},
     &          MAXLAY   = {MAXLAY},
     &          MAXSTA   = {MAXSTA},
     &          MAXCL    = {MAXCL})""".format(
        MAXEVE=MAXEVE,
        MAXDATA=MAXDATA,
        MAXEVE0=MAXEVE0,
        MAXDATA0=MAXDATA0,
        MAXLAY=MAXLAY,
        MAXSTA=MAXSTA,
        MAXCL=MAXCL,
    )
    # Remove the leading empty line
    hypodd_inc = hypodd_inc[1:]

    # Rename default hypoDD.inc file and write out new hypoDD.inc file
    hypodd_inc_rename_command = 'mv %s/hypoDD.inc %s/hypoDD_default.inc' % (include_dir, include_dir)
    subprocess.call(hypodd_inc_rename_command, shell=True)
    with open(include_dir + 'hypoDD.inc', "w") as open_file:
        hypodd_inc_file = open_file.write(hypodd_inc)
    open_file.close()

# Clean and recompile within src dir
print('Compiling ...')
clean_command = 'make clean -C %s' % src_dir
make_command = 'make -C %s' % src_dir
subprocess.call(clean_command, shell=True)
subprocess.call(make_command, shell=True)

if run_ncsn2pha:
    print('Running ncsn2pha ...')
    ncsn2pha_command = './hypoDDver/HYPODD/src/ncsn2pha/ncsn2pha %s %s' % (arc_file, phase_file)
    subprocess.call(ncsn2pha_command, shell=True)

if run_ph2dt:
    print('Running ph2dt ...')
    # Now configure ph2dt input file
    ph2dt_inp = """
* ph2dt.inp - input control file for program ph2dt
* Input station file:
{stlist_file}
* Input phase file:
{phase_file}
*MINWGHT: min. pick weight allowed [0]
*MAXDIST: max. distance in km between event pair and stations [200]
*MAXSEP: max. hypocentral separation in km [10]
*MAXNGH: max. number of neighbors per event [10]
*MINLNK: min. number of links required to define a neighbor [8]
*MINOBS: min. number of links per pair saved [8]
*MAXOBS: max. number of links per pair saved [20]
*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS
   {MINWGHT}      {MAXDIST}     {MAXSEP}     {MAXNGH}     {MINLNK}      {MINLNK}     {MAXOBS}""".format(
        stlist_file=stlist_file,
        phase_file=phase_file,
        MINWGHT=MINWGHT,
        MAXDIST=MAXDIST,
        MAXSEP=MAXSEP,
        MAXNGH=MAXNGH,
        MINLNK=MINLNK,
        MINOBS=MINOBS,
        MAXOBS=MAXOBS
    )
    # Remove the leading empty line
    ph2dt_inp = ph2dt_inp[1:]

    # Write out ph2dt.inp in inputs directory
    with open(hypoddver_dir + 'INP/ph2dt.inp', "w") as open_file:
        ph2dt_inp_file = open_file.write(ph2dt_inp)
    open_file.close()

    # Run ph2dt
    ph2dt_command = './hypoDDver/HYPODD/src/ph2dt/ph2dt ./hypoDDver/INP/ph2dt.inp'
    subprocess.call(ph2dt_command, shell=True)

    # Move ph2dt outputs to "hypoDDver/OUT_ph2dt/"
    to_move = ['ph2dt.log', 'station.sel', 'event.sel', 'event.dat', 'dt.ct']
    for filename in to_move:
        move_command = 'mv %s%s %sOUT_ph2dt/%s' % (main_dir, filename, hypoddver_dir, filename)
        remove_command = 'rm %s%s' % (main_dir, filename)
        subprocess.call(move_command, shell=True)

    # Stitch event.dat and detection.dat to create event_hypodd.dat
    with open(hypoddver_dir + 'OUT_ph2dt/event.dat') as template_list_file:
        template_list = template_list_file.read()
    with open(detection_list_filepath) as detection_list_file:
        detection_list = detection_list_file.read()
    with open(hypoddver_dir + 'OUT_ph2dt/event_hypodd.dat', 'w') as event_hypodd_file:
        event_hypodd_file.write(template_list + detection_list)

if run_hypodd:
    print('Running hypoDD ...')

    # Prepare velocity model information
    vzmodel = pd.read_csv(vzmodel_file, header=None)
    TOP = []
    VELP = []
    RAT = []
    for i in range(len(vzmodel.values)):
        line = vzmodel.values[i][0]
        if correct_depths:
            TOP.append(str(float(line.split(' ')[0]) + depth_correction))
        else:
            TOP.append(str(line.split(' ')[0]))
        VELP.append(str(line.split(' ')[1]))
        if not ps_ratio:
            RAT.append(str(line.split(' ')[2]))
        else:
            RAT.append(str(ps_ratio))
    TOP = ' '.join(TOP)
    VELP = ' '.join(VELP)
    RAT = ' '.join(RAT)

    # Configure weighting
    weight_params = [NITER, WTCCP, WTCCS, WRCC, WDCC, WTCTP, WTCTS, WRCT, WDCT, DAMP]
    # Check if each parameter has the same number of values
    if not all(len(weight_param) == len(weight_params[0]) for weight_param in weight_params[1:]):
        raise ValueError('Not all weight parameters have the same length!')
    # Now craft rows accordingly
    weight_table = ''
    for i in range(len(weight_params[0])):
        weight_table = weight_table + ' %2d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %5d' % tuple(
            [weight_param[i] for weight_param in weight_params])
        if i != len(weight_params[0]) - 1:
            weight_table = weight_table + '\n'

    # Now configure hypoDD input file
    hypodd_inp = """
hypoDD_2
* Parameter control file hypoDD.inp:
*--- INPUT FILE SELECTION
* filename of cross-corr diff. time input(blank if not available):
{dtcc_file}
* filename of catalog travel time input(blank if not available):
{dtct_file}
* filename of initial hypocenter input:
{evlist_file}
* filename of initial hypocenter input:
{stlist_file}
*
*--- OUTPUT FILE SELECTION
* filename of initial hypocenter output (if blank: output to hypoDD.loc):
{loc_file}
* filename of relocated hypocenter output (if blank: output to hypoDD.reloc):
{reloc_file}
* filename of station residual output (if blank: no output written):
{sta_file}
* filename of data residual output (if blank: no output written):
{res_file}
* filename of takeoff angle output (if blank: no output written):
{srcparam_file}
*
*--- DATA SELECTION:
* IDAT IPHA DIST
    {IDAT}     {IPHA}     {DIST}
*
*--- EVENT CLUSTERING:
* OBSCC OBSCT MINDS MAXDS MAXGAP 
    {OBSCC}    {OBSCT}    {MINDS}    {MAXDS}    {MAXGAP}        
*
*--- SOLUTION CONTROL:
* ISTART ISOLV IAQ NSET
    {ISTART}    {ISOLV}    {IAQ}    {NSET} 
*
*--- DATA WEIGHTING AND REWEIGHTING:
* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP
{weight_table}
*
*--- FORWARD MODEL SPECIFICATIONS:
* IMOD
1
* TOP: depths of top of layer (km) 
{TOP}
* VELP: layer P velocities (km/s)
{VELP}
* RAT: p/s ratio 
{RAT}
*
*--- CLUSTER/EVENT SELECTION:
* CID
    {CID}
* ID
{ID}
* end of hypoDD.inp""".format(
        dtcc_file=dtcc_file,
        dtct_file=dtct_file,
        evlist_file=evlist_file,
        stlist_file=stlist_file,
        loc_file=loc_file,
        reloc_file=reloc_file,
        sta_file=sta_file,
        res_file=res_file,
        srcparam_file=srcparam_file,
        IDAT=IDAT,
        IPHA=IPHA,
        DIST=DIST,
        OBSCC=OBSCC,
        OBSCT=OBSCT,
        MINDS=MINDS,
        MAXDS=MAXDS,
        MAXGAP=MAXGAP,
        ISTART=ISTART,
        ISOLV=ISOLV,
        IAQ=IAQ,
        NSET=NSET,
        weight_table=weight_table,
        TOP=TOP,
        VELP=VELP,
        RAT=RAT,
        CID=CID,
        ID=ID
    )
    # Remove the leading empty line
    hypodd_inp = hypodd_inp[1:]

    # Write out hypoDD.inp in inputs directory
    with open(hypoddver_dir + 'INP/hypoDD.inp', "w") as open_file:
        hypodd_inp_file = open_file.write(hypodd_inp)
    open_file.close()

    # Run hypoDD
    hypodd_command = './hypoDDver/HYPODD/src/hypoDD/hypoDD ./hypoDDver/INP/hypoDD.inp'
    subprocess.call(hypodd_command, shell=True)

    # # Move hypodd outputs to "hypoDDver/OUT_hypodd/"
    # to_move = ['hypodd.log']
    # for filename in to_move:
    #     move_command = 'mv %s%s %sOUT_hypodd/%s' % (main_dir,filename,hypoddver_dir,filename)
    #     remove_command = 'rm %s%s' % (main_dir,filename)
    #     subprocess.call(move_command, shell=True)

    # Write out catalog objects from event.sel, hypodd.loc and hypodd.reloc
    print('Now writing catalog objects from hypodd outputs...')
    hypodd_loc = loc2cat(loc_file, event_id_mapper, catalog, type='loc', depth_correction=depth_correction)
    hypodd_reloc = loc2cat(reloc_file, event_id_mapper, catalog, type='reloc', depth_correction=depth_correction)
    writer(output_dir + 'hypodd_loc.xml', hypodd_loc)
    writer(output_dir + 'hypodd_reloc.xml', hypodd_reloc)

if plot_results:
    # Read in hypodd input events and remove overlaps with hypodd output of unrelocated events
    hypodd_in = dat2cat(evlist_file, event_id_mapper, catalog, depth_correction=depth_correction)
    hypodd_loc_filtered = remove_catalog_repeats(hypodd_loc, hypodd_reloc)
    writer(output_dir + 'hypodd_loc_filtered.xml', hypodd_loc_filtered)

    # Extract lat, lon and depth
    hypodd_in_lats = [event.origins[0].latitude for event in hypodd_in]
    hypodd_in_lons = [event.origins[0].longitude for event in hypodd_in]
    hypodd_in_deps = [event.origins[0].depth for event in hypodd_in]
    hypodd_loc_lats = [event.origins[0].latitude for event in hypodd_loc]
    hypodd_loc_lons = [event.origins[0].longitude for event in hypodd_loc]
    hypodd_loc_deps = [event.origins[0].depth for event in hypodd_loc]
    hypodd_loc_filtered_lats = [event.origins[0].latitude for event in hypodd_loc_filtered]
    hypodd_loc_filtered_lons = [event.origins[0].longitude for event in hypodd_loc_filtered]
    hypodd_loc_filtered_deps = [event.origins[0].depth for event in hypodd_loc_filtered]
    hypodd_reloc_lats = [event.origins[0].latitude for event in hypodd_reloc]
    hypodd_reloc_lons = [event.origins[0].longitude for event in hypodd_reloc]
    hypodd_reloc_deps = [event.origins[0].depth for event in hypodd_reloc]

    # Mammoth
    lat_min = 37.54
    lat_max = 37.71
    lon_min = -119.20
    lon_max = -118.98
    dep_min = -5
    dep_max = 30

    # # Augustine
    # lat_min = 59.3
    # lat_max = 59.42
    # lon_min = -153.6
    # lon_max = -153.2
    # dep_min = -3.1
    # dep_max = 10

    # ## Redoubt
    # lat_min = 60.2852
    # lat_max = 60.6800
    # lon_min = -153.1438
    # lon_max = -152.3438
    # dep_min = -5
    # dep_max = 15

    # Now do simple plotting
    fig, ax = plt.subplots(2, 3, figsize=(8, 8))
    ms = 2
    # [LON VS LAT] Plot selected events among input events (blue against black)
    ax[0, 0].plot(hypodd_in_lons, hypodd_in_lats, 'k.', markersize=ms, label='input events')
    ax[0, 0].plot(hypodd_loc_lons, hypodd_loc_lats, 'b.', markersize=ms, label='relocatable')
    ax[0, 0].axis('square')
    ax[0, 0].axis([lon_min, lon_max, lat_min, lat_max])
    # ax[0,0].legend(loc='upper right')
    ax[0, 0].set_title('LON vs LAT (before hypodd)')
    ax[0, 0].set_ylabel('Latitude')
    ax[0, 0].set_xlabel('Longitude')
    ax[0, 0].axis('square')
    # [LON VS LAT] Plot relocated events among retained events (red against blue)
    # ax[1,0].plot(event_sel_lons,event_sel_lats,'k.',markersize=ms,label='input events')
    ax[1, 0].plot(hypodd_loc_filtered_lons, hypodd_loc_filtered_lats, 'b.', markersize=ms, label='relocatable')
    ax[1, 0].plot(hypodd_reloc_lons, hypodd_reloc_lats, 'r.', markersize=ms, label='relocated')
    ax[1, 0].axis('square')
    ax[1, 0].axis([lon_min, lon_max, lat_min, lat_max])
    # ax[1,0].legend(loc='upper right')
    ax[1, 0].set_title('LON vs LAT (after hypodd)')
    ax[1, 0].set_ylabel('Latitude')
    ax[1, 0].set_xlabel('Longitude')
    # [LON VS DEP] Plot selected events among input events (blue against black)
    ax[0, 1].plot(hypodd_in_lons, hypodd_in_deps, 'k.', markersize=ms, label='input events')
    ax[0, 1].plot(hypodd_loc_lons, hypodd_loc_deps, 'b.', markersize=ms, label='relocatable')
    ax[0, 1].axis([lon_min, lon_max, dep_min, dep_max])
    ax[0, 1].invert_yaxis()
    ax[0, 1].legend(loc='lower right')
    ax[0, 1].set_title('LON vs DEP (before hypodd)')
    ax[0, 1].set_ylabel('Depth (km)')
    ax[0, 1].set_xlabel('Longitude')
    # [LON VS DEP] Plot relocated events among retained events (red against blue)
    # ax[1,1].plot(event_sel_lons, event_sel_deps, 'k.', markersize=ms,label='input events')
    ax[1, 1].plot(hypodd_loc_filtered_lons, hypodd_loc_filtered_deps, 'b.', markersize=ms, label='relocatable')
    ax[1, 1].plot(hypodd_reloc_lons, hypodd_reloc_deps, 'r.', markersize=ms, label='relocated')
    ax[1, 1].axis([lon_min, lon_max, dep_min, dep_max])
    ax[1, 1].invert_yaxis()
    ax[1, 1].legend(loc='lower right')
    ax[1, 1].set_title('LON vs DEP (after hypodd)')
    ax[1, 1].set_ylabel('Depth (km)')
    ax[1, 1].set_xlabel('Longitude')
    # [LAT VS DEP] Plot selected events among input events (blue against black)
    ax[0, 2].plot(hypodd_in_lats, hypodd_in_deps, 'k.', markersize=ms, label='input events')
    ax[0, 2].plot(hypodd_loc_lats, hypodd_loc_deps, 'b.', markersize=ms, label='relocatable')
    ax[0, 2].axis([lat_min, lat_max, dep_min, dep_max])
    ax[0, 2].invert_yaxis()
    ax[0, 2].legend(loc='lower right')
    ax[0, 2].set_title('LAT vs DEP (before hypodd)')
    ax[0, 2].set_ylabel('Depth (km)')
    ax[0, 2].set_xlabel('Latitude')
    # [LAT VS DEP] Plot relocated events among retained events (red against blue)
    # ax[1,2].plot(event_sel_lats, event_sel_deps, 'k.', markersize=ms,label='input events')
    ax[1, 2].plot(hypodd_loc_filtered_lats, hypodd_loc_filtered_deps, 'b.', markersize=ms, label='relocatable')
    ax[1, 2].plot(hypodd_reloc_lats, hypodd_reloc_deps, 'r.', markersize=ms, label='relocated events')
    ax[1, 2].axis([lat_min, lat_max, dep_min, dep_max])
    ax[1, 2].invert_yaxis()
    ax[1, 2].legend(loc='lower right')
    ax[1, 2].set_title('LAT vs DEP (after hypodd)')
    ax[1, 2].set_ylabel('Depth (km)')
    ax[1, 2].set_xlabel('Latitude')
    plt.tight_layout()
    fig.show()