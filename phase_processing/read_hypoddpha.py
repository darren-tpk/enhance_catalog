# Read a hypoddpha catalog with reference to a hypoi source (for component info)
def read_hypoddpha(hypoi_file, hypoddpha_file, channel_convention=False):

    # Import dependency
    import numpy as np
    import pandas as pd
    from obspy import Catalog

    # Print commencement
    print('Reading hypoddpha catalog...')

    # Read hypoddpha file into a python catalog using obspy read_events
    from obspy import read_events
    catalog_raw = read_events(hypoddpha_file, "HYPODDPHA")

    # Extract event id, station, network and component from hypoi file
    # Read hypoi data file
    phase_info = pd.read_csv(hypoi_file, header=None, skipfooter=1, engine='python')
    phase_info = phase_info[0]

    # Initialize list for phase id and phase channel components
    phase_codes = []
    phase_comps = []

    # Loop through hypoi lines
    for i in range(len(phase_info)):

        # Extract current line
        line = phase_info[i]

        # If the current line describes the event, get event id
        if '1' in line[0] or '2' in line[0]:
            event_id = line[136:146].strip()
            # Some events do not have depth - we avoid these events
            if (line[31:36]).isspace():
                has_depth = 0
            else:
                has_depth = 1

        # If the current line describes a phase pick
        elif line[0].isalpha() == True:

            # Do not process phases from events without depth
            if has_depth == 0:
                continue

            # Extract real component / channel
            phase_comp = line[9:12]

            # If using channel conventions, check if component and pick comply
            if channel_convention:

                # Check if vertical component has P arrival
                if (line[11] == 'Z') and  (line[14] == 'P'):

                    # Craft a phase "code" that stores net-sta-chan-pick info for the phase line
                    phase_code = event_id + '-' + line[5:7] + '-' + line[0:5].strip() + '-P'

                    # Append phase "code" and real channel (corresponding indices)
                    phase_codes.append(phase_code)
                    phase_comps.append(phase_comp)

                # Check if horizontal component has S arrival
                elif (line[11] in ['E','N','1','2']) and (line[47] == 'S'):

                    # Craft a phase "code" that stores net-sta-chan-pick info for the phase line
                    phase_code = event_id + '-' + line[5:7] + '-' + line[0:5].strip() + '-S'

                    # Append phase "code" and real channel (corresponding indices)
                    phase_codes.append(phase_code)
                    phase_comps.append(phase_comp)


            # If we are not using channel conventions, then we can have P and S picks on any component
            else:

                # Check for P pick on this component's line
                if line[14] == 'P':

                    # Craft a phase "code" that stores net-sta-chan-pick info for the phase line
                    phase_code = event_id + '-' + line[5:7] + '-' + line[0:5].strip() + '-P'

                    # Append phase "code" and real channel (corresponding indices)
                    phase_codes.append(phase_code)
                    phase_comps.append(phase_comp)

                # Check for S pick on this component's line
                if line[47] == 'S':

                    # Craft a phase "code" that stores net-sta-chan-pick info for the phase line
                    phase_code = event_id + '-' + line[5:7] + '-' + line[0:5].strip() + '-S'

                    # Append phase "code" and real channel (corresponding indices)
                    phase_codes.append(phase_code)
                    phase_comps.append(phase_comp)

        # Move to next event
        elif ' ' in line[0]:
            continue

    # Convert to array for np.argwhere in later part of code
    phase_codes = np.array(phase_codes)
    phase_comps = np.array(phase_comps)

    # Initialize a counter to make sure the indices are referenced accurately
    count = 0

    # Quick loop to fix network and station entries, channel code and remove EQs without picks

    # Initialize catalog
    catalog = Catalog()

    # For each event in our catalog
    for i in range(0,len(catalog_raw)):

        # Get number of picks
        num_picks = len(catalog_raw[i].picks)

        # Loop through each event's pick list
        for j in range(num_picks):

            # Extract info for comparison
            event_id = str(catalog_raw[i].resource_id)[16:]

            # If the network field is empty, we dissect the station code to get net and sta
            if len(catalog_raw[i].picks[j].waveform_id.network_code) == 0:
                network = catalog_raw[i].picks[j].waveform_id.station_code[0:2]
                station = catalog_raw[i].picks[j].waveform_id.station_code[2:]

            # If the network field is populated, we can grab net and sta directly
            else:
                network = catalog_raw[i].picks[j].waveform_id.network_code
                station = catalog_raw[i].picks[j].waveform_id.station_code

            # We extract the orientation of the component
            channel_dir = catalog_raw[i].picks[j].waveform_id.channel_code[-1]

            # ObsPy's read_events for hypoddpha defaults all P picks to Z channel and all S picks to an N channel
            # We craft their respective phase "codes" to compare with our reference hypoi
            if channel_dir == 'Z':
                phase_code = event_id + '-' + network + '-' + station + '-P'
            elif channel_dir == 'N':
                phase_code = event_id + '-' + network + '-' + station + '-S'
            else:
                raise ValueError('WARNING: Unrecognized channel directional component:',channel_dir)

            # Fix network and station entries
            catalog_raw[i].picks[j].waveform_id.network_code = network
            catalog_raw[i].picks[j].waveform_id.station_code = station

            # Fix channel code by checking if the phase code matches
            if phase_codes[count] == phase_code:
                real_channel = phase_comps[count]
                catalog_raw[i].picks[j].waveform_id.channel_code = real_channel
            else:
                raise IndexError('WARNING: Phase code did not match. Event ID:', event_id,'Phase Code:',phase_code)

            # Add one to counter for subsequent loops
            count = count + 1

        # Do not add events without picks
        if num_picks != 0:
            catalog.append(catalog_raw[i])

    # Print out successful check for user
    if len(phase_codes) == count:
        print('Phase codes tally with count. Catalog networks, stations and channels are updated.')

    return catalog
