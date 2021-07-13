import numpy as np
import pandas as pd
from obspy import UTCDateTime

# Define function to convert input to float if applicable, string if not, and nan if empty
def formater(input):
    if ' ' * len(input) in input:  # If empty
        output = float("nan")  # We choose nan so that scalings in Class will not be affected
    else:
        try:
            output = int(float(input))  # Convert to int if possible
        except ValueError:
            output = input.strip()  # Keep as string if not, trimming white spaces
    return output

# Define function to print out all fields in Classes EventPhase, Phase and More
# Phase and More are shown as their number of fields in EventPhase
def printer(Class):
    for item in vars(Class).items():
        if 'phases' in item[0]:
            print("{0:40}: {1} phases in total".format(item[0], len(Class.phases)))
        elif 'More' in item[0]:
            print("{0:40}: {1} more details".format(item[0], len(vars(Class.More))))
        else:
            print("{0:40}: {1}".format(item[0], item[1]))

# Define Class called EventPhase (not called Event, to avoid clash with obspy)
# Stores key information, array of Phases, and More
class EventPhase:
    def __init__(self, line):

        # Construct datestring
        datestring = "%s-%s-%s %s:%s:%s.%s" % (
        line[0:4], line[4:6], line[6:8], line[8:10], line[10:12], line[12:14], line[14:16])

        # Construct latitude (convert to decimal degrees)
        if 'S' in line[18]:
            latitude = -1 * (formater(line[16:18]) + (0.01 * formater(line[19:23]) / 60))
        elif ' ' in line[18]:
            latitude = (formater(line[16:18]) + (0.01 * formater(line[19:23]) / 60))
        else:
            print("LATITUDE ERROR")

        # Construct longitude (convert to decimal degrees)
        if 'E' in line[26]:
            longitude = (formater(line[23:26]) + (0.01 * formater(line[27:31]) / 60))
        elif ' ' in line[26]:
            longitude = -1 * (formater(line[23:26]) + (0.01 * formater(line[27:31]) / 60))
        else:
            print("LONGITUDE ERROR")

        # Fill in key information for Event Phase
        self.event_id = formater(line[136:146])
        self.datetime = UTCDateTime(datestring)
        self.latitude = latitude # decimal degrees
        self.longitude = longitude # decimal degrees
        self.depth = 0.01 * formater(line[31:36])  # km
        self.magnitude_label = formater(line[146])  # preferred magnitude label
        self.magnitude = 0.01 * formater(line[147:150])  # preferred magnitude
        self.total_weights = 0.1 * formater(line[150:154])  # no. of. readings
        self.all_Phase = float("nan")  # placeholder
        self.More = float("nan")  # placeholder


# Define Class called More
# Stores all additional information within Event
class More:
    def __init__(self, line):
        # Fill in additional information
        self.amplitude_magnitude = 0.01 * formater(line[36:39])
        self.number_of_weighted_ps = formater(line[39:42])  # with weights > 0.1
        self.max_azimuth_gap = formater(line[42:45])  # deg
        self.distance_nearest_station = formater(line[45:48])  # km
        self.rms_travel_time_residual = 0.01 * formater(line[48:52])
        self.largest_principal_error_azimuth = formater(line[52:55])  # deg e of n
        self.largest_principal_error_dip = formater(line[55:57])  # deg
        self.largest_principal_error_size = 0.01 * formater(line[57:61])  # km
        self.intermediate_principal_error_azimuth = formater(line[61:64])  # deg
        self.intermediate_principal_error_dip = formater(line[64:66])  # deg
        self.intermediate_principal_error_size = 0.01 * formater(line[66:70])  # km
        self.coda_duration_magnitude = 0.01 * formater(line[70:73])  # [empty]
        self.location_remark = formater(line[73:76])  # [empty]
        self.smallest_principal_error_size = 0.01 * formater(line[76:80])
        self.analyst_remark = formater(line[80])  # [empty]
        self.program_remark = formater(line[81])  # [empty]
        self.number_of_s = formater(line[82:85])  # with weights > 0.1
        self.horizontal_error = 0.01 * formater(line[85:89])  # km
        self.vertical_error = 0.01 * formater(line[89:93])  # km
        self.number_of_p = formater(line[93:96])  # first motions
        self.amplitude_magnitude_total_weights = 0.1 * formater(line[96:100])  # no. of readings # [empty]
        self.duration_magnitude_total_weights = 0.1 * formater(line[100:104])  # no. of readings # [empty]
        self.median_abs_diff_amplitude_magnitudes = 0.01 * formater(line[104:107])  # [empty]
        self.median_abs_diff_duration_magnitudes = 0.01 * formater(line[107:110])  # [empty]
        self.crust_delay_model = formater(line[110:113])  # 3 letter code # [empty]
        self.authority_code = formater(line[113])  # network that furnished info # [empty]
        self.most_common_ps_source = formater(line[114])  # see manual # [empty]
        self.most_common_duration_source = formater(line[115])  # [empty]
        self.most_common_amplitude_source = formater(line[116])  # [empty]
        self.primary_coda_duration_magnitude_type = formater(line[117])  # [empty]
        self.number_of_valid_ps = formater(line[118:121])  # weight >= 0
        self.primary_amplitude_magntiude_type = formater(line[121])  # [empty]
        self.external_magnitude_label = formater(line[122])
        self.external_magnitude = 0.01 * formater(line[123:126])
        self.external_magnitude_total_weights = 0.1 * formater(line[126:129])  # no. of readings
        self.alternate_amplitude_magnitude_label = formater(line[129])  # or type code # [empty]
        self.alternate_amplitude_magnitude = 0.01 * formater(line[130:133])  # [empty]
        self.alternate_amplitude_total_weights = 0.1 * formater(line[133:136])  # no. of readings # [empty]
        self.alternate_coda_duration_magnitude_label = formater(line[154])  # [empty]
        self.alternate_coda_duration_magnitude = 0.01 * formater(line[155:158])  # [empty]
        self.alternate_coda_duration_total_weights = 0.1 * formater(line[158:162])  # no. of readings # [empty]


# Define Class called Phase
# Stores all phase-related information for phase lines
class Phase:
    def __init__(self, line):

        # Construct datestring
        datestring = "%s-%s-%s %s:%s" % (line[17:21], line[21:23], line[23:25], line[25:27], line[27:29])

        # For P arrival
        p_seconds = formater(line[29:32])
        if np.isnan(p_seconds) == True:
            p_datetime = float("nan")
        elif p_seconds >= 60:
            p_minutes = formater(line[27:29])
            p_minutes = p_minutes + np.floor(p_seconds / 60)
            p_seconds = p_seconds % 60
            p_datestring = "%s-%s-%s %s:%s:%s.%s" % (
            line[17:21], line[21:23], line[23:25], line[25:27], p_minutes, p_seconds, line[32:34])
            p_datetime = UTCDateTime(p_datestring)
        else:
            p_datestring = "%s-%s-%s %s:%s:%s.%s" % (
            line[17:21], line[21:23], line[23:25], line[25:27], line[27:29], p_seconds, line[32:34])
            p_datetime = UTCDateTime(p_datestring)

        # For S arrival
        s_seconds = formater(line[41:44])
        if np.isnan(s_seconds) == True:
            s_datetime = float("nan")
        elif s_seconds >= 60:
            s_minutes = formater(line[27:29])
            s_minutes = s_minutes + np.floor(s_seconds / 60)
            s_seconds = s_seconds % 60
            s_datestring = "%s-%s-%s %s:%s:%s.%s" % (
            line[17:21], line[21:23], line[23:25], line[25:27], s_minutes, s_seconds, line[32:34])
            s_datetime = UTCDateTime(s_datestring)
        else:
            s_datestring = "%s-%s-%s %s:%s:%s.%s" % (
            line[17:21], line[21:23], line[23:25], line[25:27], line[27:29], s_seconds, line[44:46])
            s_datetime = UTCDateTime(s_datestring)

        # Station Info
        self.datetime = UTCDateTime(datestring) # remove?
        self.station = formater(line[0:5])
        self.network = formater(line[5:7])
        self.component_single = formater(line[8])
        self.component_triple = formater(line[9:12])

        # P arrival info
        self.p_first_motion = formater(line[15])
        self.p_datetime = p_datetime
        self.p_arrival = 0.01 * formater(line[29:34])  # second of P arrival
        self.p_travel_time_residual = 0.01 * formater(line[34:38])
        self.p_delay = 0.01 * formater(line[66:70])
        self.p_weight = 0.01 * formater(line[38:41])
        self.p_weight_code = formater(line[16])
        self.p_remark = formater(line[13:15])

        # S arrival info
        self.s_datetime = s_datetime
        self.s_arrival = 0.01 * formater(line[41:46])  # second of S arrival

        # If we have a longer line, we input more info
        if len(line) > 46:  # Hard coded condition

            # More S info
            self.s_travel_time_residual = 0.01 * formater(line[50:54])
            self.s_delay = 0.01 * formater(line[70:74])
            self.s_weight = 0.01 * formater(line[63:66])
            self.s_weight_code = formater(line[49])
            self.s_remark = formater(line[46:48])

            # Amplitude and coda duration info
            self.amplitude_peak_to_peak = 0.01 * formater(line[54:61])  # normally peak to peak
            self.amplitude_units = formater(line[61:63])  # 0 = PP mm, 1 = 0 to peak mm (UCB), 2 = digital counts
            self.amplitude_magnitude_label = formater(line[110])
            self.amplitude_magnitude = 0.01 * formater(line[97:100])
            self.amplitude_magnitude_weight_code = formater(line[81])
            self.period_of_measurement = 0.01 * formater(line[83:86])
            self.coda_duration = formater(line[87:91])  # s
            self.duration_magnitude_label = formater(line[109])
            self.duration_magnitude = 0.01 * formater(line[94:97])
            self.duration_magnitude_weight_code = formater(line[82])

            # Direction and Location
            self.epicentral_distance = 0.1 * formater(line[74:78])
            self.emergence_angle = formater(line[78:81])  # at source
            self.azimuth_to_station = formater(line[91:94])  # deg E of N

            # Misc
            self.p_importance = 0.001 * formater(line[100:104])
            self.s_importance = 0.001 * formater(line[104:108])
            self.station_location = formater(line[111:113])
            self.station_remark = formater(line[86])
            self.data_source_code = formater(line[108])


# Master function that loops through all lines in our hypoinverse file to construct catalog
def construct_phase_catalog(phase_dir):
    # Phase is in hypoinverse Y2000 summary format for headers, Y2000 (station) archive format for phases
    phase_info = pd.read_csv(phase_dir, header=None, skipfooter=1, engine='python')
    phase_info = phase_info[0]

    # Initialize array for Events
    all_Event = np.empty((0, 0))

    # Loop through every line in data
    for i in range(len(phase_info)):

        # Extract line
        line = phase_info[i]

        # Use line to fill in information accordingly

        # If it is a new event line...
        if '1' in line[0] or '2' in line[0]:
            # Create new EventPhase
            new_Event = EventPhase(line)
            new_Event.More = More(line)
            # Initialize array for Phases
            all_Phase = np.empty((0, 0))

        # If it is a new phase line...
        elif line[0].isalpha() == True:
            new_Phase = Phase(line)
            all_Phase = np.append(all_Phase, new_Phase)

        # Else if it is a blank line, we replace the placeholder in Event with the Phase array
        elif ' ' in line[0]:
            new_Event.all_Phase = all_Phase
            # Lastly we append the completed Event to the Event array
            all_Event = np.append(all_Event, new_Event)

    return all_Event