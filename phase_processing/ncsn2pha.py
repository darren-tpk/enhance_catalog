# Convert phase data file from hypoi format to hypoddpha
def ncsn2pha(hypoi_file, hypoddpha_file, channel_convention=False):

    # Import dependency
    import pandas as pd

    # hypoinverse Y2000 summary format for headers
    # Y2000 (station) archive format for phases

    # Read input file
    data_frame = pd.read_csv(hypoi_file, header=None, skipfooter=1, engine='python')
    input = data_frame[0]

    # Open output file, to sequentially enter quake and phase lines
    # This will be in hypoddpha format
    output = open(hypoddpha_file, 'w')

    # Loop through every line in input file
    for i in range(len(input)):

        # Extract line
        line = input[i]

        # Check what type of line it is

        # If it is a new earthquake location line
        if '19' in line[0:2] or '20' in line[0:2]:

            # Construct latitude (convert to decimal degrees)
            if 'S' in line[18]:
                latitude = -1 * (float(line[16:18]) + (0.01 * float(line[19:23]) / 60))
            elif ' ' in line[18]:
                latitude = (float(line[16:18]) + (0.01 * float(line[19:23]) / 60))
            else:
                print("LATITUDE ERROR")

            # Construct longitude (convert to decimal degrees)
            if 'E' in line[26]:
                longitude = (float(line[23:26]) + (0.01 * float(line[27:31]) / 60))
            elif ' ' in line[26]:
                longitude = -1 * (float(line[23:26]) + (0.01 * float(line[27:31]) / 60))
            else:
                print("LONGITUDE ERROR")

            # Organize earthquake location line
            year = line[0:4]
            month = line[4:6]
            day = line[6:8]
            hour = line[8:10]
            minute = line[10:12]
            second = "%4.2f" % (0.01 * float(line[12:16]))
            latitude = "%8.4f" % latitude
            longitude = "%9.4f" % longitude
            # Some events do not have depth
            if (line[31:36]).isspace():
                depth = None
            else:
                depth = "%5.2f" % (0.01 * float(line[31:36]))
            magnitude = "%3.2f" % (0.01 * float(line[147:150]))
            rms_error = "%4.2f" % (0.01 * float(line[48:52]))  # travel time residual
            horizontal_error = "%4.2f" % (0.01 * float(line[85:89]))
            vertical_error = "%4.2f" % (0.01 * float(line[89:93]))
            event_id = line[136:146]

            # Make sure depth is populated when writing entry
            if depth is not None:
                # Construct entry
                entry = ('# ' + year + ' ' + month + ' ' + day + ' ' + hour + ' ' +
                        minute + ' ' + second + ' ' + latitude + ' ' + longitude +
                        ' ' + depth + ' ' + magnitude + ' ' + horizontal_error + ' '
                        + vertical_error + ' ' + rms_error + ' ' + event_id)
                # Write entry in output
                output.write(entry + '\n')

        # If it is a new phase line with potential S picks:
        elif line[0:2].isalpha() and len(line) > 46:

            # If there is no depth defined for the event, skip all of its phases
            if depth is None:
                continue

            # Extract phase information
            network = line[5:7]
            station = (line[0:5]).strip()  # input is left justified
            phase_minute = int(line[27:29])  # minute for p & s

            # If using channel conventions, check if component and pick comply
            if channel_convention:

                # Check if vertical component has P arrival
                if (line[11] == 'Z') and  (line[14] == 'P'):

                    # Convert P weight
                    p_weight_code = int(line[16])
                    if p_weight_code <= 9:
                        if p_weight_code == 0:
                            p_weight = 1.0
                        elif p_weight_code == 1:
                            p_weight = 0.5
                        elif p_weight_code == 2:
                            p_weight = 0.2
                        elif p_weight_code == 3:
                            p_weight = 0.1
                        else:
                            p_weight = 0.0

                        # Check P importance and adjust P weight
                        try:
                            p_importance = 0.001 * float(line[100:104])
                            if p_importance > 0.5:
                                p_weight = -1 * p_weight
                        except:
                            pass
                        p_weight = "%6.3f" % p_weight

                        # Get P travel time
                        p_second = "%5.2f" % (0.01 * float(line[29:34]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute+60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        p_arrival = (arrival_minute*60) + (float(p_second)-float(second))
                        p_arrival = "%8.3f" % p_arrival

                        # Construct entry
                        entry = network + station + ' ' + p_arrival + ' ' + p_weight + ' P'

                        # Write entry in output
                        output.write(entry + '\n')

                # Check if horizontal component has S arrival
                elif (line[11] in ['E','N','1','2']) and (line[47] == 'S'):

                    # Convert S weight
                    s_weight_code = int(line[49])
                    if s_weight_code <= 9:
                        if s_weight_code == 0:
                            s_weight = 1.0
                        elif s_weight_code == 1:
                            s_weight = 0.5
                        elif s_weight_code == 2:
                            s_weight = 0.2
                        elif s_weight_code == 3:
                            s_weight = 0.1
                        else:
                            s_weight = 0.0

                        # Check S importance and adjust S weight
                        s_importance = 0.001 * float(line[104:108])
                        if s_importance > 0.5:
                            s_weight = -1 * s_weight
                        s_weight = "%6.3f" % s_weight

                        # Get S travel time
                        s_second = "%5.2f" % (0.01 * float(line[41:46]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute + 60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        s_arrival = (arrival_minute * 60) + (float(s_second) - float(second))
                        s_arrival = "%8.3f" % s_arrival

                        # Construct entry
                        entry = network + station + ' ' + s_arrival + ' ' + s_weight + ' S'

                        # Write entry in output
                        output.write(entry + '\n')

            # If we are not using channel conventions, then we can have P and S picks on any component
            else:

                # Check if this line has a P pick
                if line[14] == 'P':

                    # Convert P weight
                    p_weight_code = int(line[16])
                    if p_weight_code <= 9:
                        if p_weight_code == 0:
                            p_weight = 1.0
                        elif p_weight_code == 1:
                            p_weight = 0.5
                        elif p_weight_code == 2:
                            p_weight = 0.2
                        elif p_weight_code == 3:
                            p_weight = 0.1
                        else:
                            p_weight = 0.0

                        # Check P importance and adjust P weight
                        try:
                            p_importance = 0.001 * float(line[100:104])
                            if p_importance > 0.5:
                                p_weight = -1 * p_weight
                        except:
                            pass
                        p_weight = "%6.3f" % p_weight

                        # Get P travel time
                        p_second = "%5.2f" % (0.01 * float(line[29:34]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute+60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        p_arrival = (arrival_minute*60) + (float(p_second)-float(second))
                        p_arrival = "%8.3f" % p_arrival

                        # Construct entry
                        entry = network + station + ' ' + p_arrival + ' ' + p_weight + ' P'

                        # Write entry in output
                        output.write(entry + '\n')

                # Check if this line has an S pick
                if line[47] == 'S':

                    # Convert S weight
                    s_weight_code = int(line[49])
                    if s_weight_code <= 9:
                        if s_weight_code == 0:
                            s_weight = 1.0
                        elif s_weight_code == 1:
                            s_weight = 0.5
                        elif s_weight_code == 2:
                            s_weight = 0.2
                        elif s_weight_code == 3:
                            s_weight = 0.1
                        else:
                            s_weight = 0.0

                        # Check S importance and adjust S weight
                        s_importance = 0.001 * float(line[104:108])
                        if s_importance > 0.5:
                            s_weight = -1 * s_weight
                        s_weight = "%6.3f" % s_weight

                        # Get S travel time
                        s_second = "%5.2f" % (0.01 * float(line[41:46]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute + 60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        s_arrival = (arrival_minute * 60) + (float(s_second) - float(second))
                        s_arrival = "%8.3f" % s_arrival

                        # Construct entry
                        entry = network + station + ' ' + s_arrival + ' ' + s_weight + ' S'

                        # Write entry in output
                        output.write(entry + '\n')

        # If it is a new phase line without potential S pick:
        elif line[0:2].isalpha() and len(line) <= 46:

            print('WARNING, i =',i,'has truncated phase line')

            # If there is no depth defined for the event, skip all of its phases
            if depth is None:
                continue

            # Extract phase information
            network = line[5:7]
            station = (line[0:5]).strip()  # input is left justified
            phase_minute = int(line[27:29])  # minute for p & s

            # If using channel conventions, check if component and pick comply
            if channel_convention:

                # Check if vertical component has P arrival
                if (line[11] == 'Z') and  (line[14] == 'P'):

                    # Convert P weight
                    p_weight_code = int(line[16])
                    if p_weight_code <= 9:
                        if p_weight_code == 0:
                            p_weight = 1.0
                        elif p_weight_code == 1:
                            p_weight = 0.5
                        elif p_weight_code == 2:
                            p_weight = 0.2
                        elif p_weight_code == 3:
                            p_weight = 0.1
                        else:
                            p_weight = 0.0

                        # Adjust P weight
                        # (note there is no P importance to check against)
                        p_weight = "%6.3f" % p_weight

                        # Get P travel time
                        p_second = "%5.2f" % (0.01 * float(line[29:34]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute+60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        p_arrival = (arrival_minute*60) + (float(p_second)-float(second))
                        p_arrival = "%8.3f" % p_arrival

                        # Construct entry
                        entry = network + station + ' ' + p_arrival + ' ' + p_weight + ' P'

                        # Write entry in output
                        output.write(entry + '\n')

            # If we are not using channel conventions, then we can have P picks on any component
            else:

                # Check if this line is a P pick
                if line[14] == 'P':

                    # Process P time if it is valid and weighted
                    p_weight_code = int(line[16])
                    if p_weight_code <= 9:
                        if p_weight_code == 0:
                            p_weight = 1.0
                        elif p_weight_code == 1:
                            p_weight = 0.5
                        elif p_weight_code == 2:
                            p_weight = 0.2
                        elif p_weight_code == 3:
                            p_weight = 0.1
                        else:
                            p_weight = 0.0

                        # Adjust P weight
                        # (note there is no P importance to check against)
                        p_weight = "%6.3f" % p_weight

                        # Get P travel time
                        p_second = "%5.2f" % (0.01 * float(line[29:34]))
                        if phase_minute < int(minute):
                            arrival_minute = phase_minute+60 - int(minute)
                        else:
                            arrival_minute = phase_minute - int(minute)
                        p_arrival = (arrival_minute*60) + (float(p_second)-float(second))
                        p_arrival = "%8.3f" % p_arrival

                        # Construct entry
                        entry = network + station + ' ' + p_arrival + ' ' + p_weight + ' P'

                        # Write entry in output
                        output.write(entry + '\n')

        # If it is shadow line or blank line, skip:
        elif '$' in line [0] or '  ' in line[0:2]:
            continue

    # Close output file at end of loop
    output.close()
