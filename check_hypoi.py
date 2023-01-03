## Check avo catalog (using the hypoi file) on whether it meets the requirements necessary

augustine_hypoi_file = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/data/avo/augustine_20050401_20060501_hypoi.txt'
redoubt_hypoi_file = '/Users/darrentpk/Desktop/GitHub/enhance_catalog/data/avo/redoubt_20080501_20090901_hypoi.txt'
## Namely we check for
#  Number of P and S phases
#  Number of observations
#  Travel time RMSE
#  Maximum Azimuthal Gap
#  Horizontal Error
#  Vertical Error

def check_hypoi(hypoi_file,
                nsta_thres=4,
                nobs_thres=4,
                rmse_thres=0.3,
                gap_thres=180,
                eh_thres=1,
                ev_thres=1,
                num_P_thres=3,
                num_S_thres=2,
                summary=True):

    # Import dependency
    import pandas as pd
    import numpy as np
    from obspy import UTCDateTime

    # hypoinverse Y2000 summary format for headers
    # Y2000 (station) archive format for phases

    # Initialize checks
    all_datetimes = []
    all_PS = []
    all_nobs = []
    all_nsta = []
    all_rmse = []
    all_gap = []
    all_dist2sta = []
    all_depth = []
    all_eh = []
    all_ev = []
    has_mag = []
    has_depth = []

    # Read input file
    data_frame = pd.read_csv(hypoi_file, header=None, skipfooter=1, engine='python')
    input = data_frame[0]

    # Loop through every line in input file
    for i in range(len(input)):

        # Extract line
        line = input[i]

        # Check what type of line it is

        # If it is a new earthquake location line
        if '19' in line[0:2] or '20' in line[0:2]:

            # construct UTCDateTime
            year = int(line[0:4])
            month = int(line[4:6])
            day = int(line[6:8])
            hour = int(line[8:10])
            minute = int(line[10:12])
            second = 0.01 * float(line[12:16])
            datetime = UTCDateTime(year,month,day,hour,minute,second)

            # extract filter params
            nobs = int(line[39:42])
            num_S = int(line[82:85])
            num_P = nobs - num_S
            rmse = 0.01 * float(line[48:52])  # travel time residual
            gap = int(line[42:45])
            dist2sta = int(line[45:48])
            eh = 0.01 * float(line[85:89])
            ev = 0.01 * float(line[89:93])

            all_datetimes.append(datetime)
            all_dist2sta.append(dist2sta)
            all_PS.append('%02dP%02dS' % (num_P,num_S))
            all_nobs.append(nobs)
            all_rmse.append(rmse)
            all_gap.append(gap)
            all_eh.append(eh)
            all_ev.append(ev)

            if (line[31:36]).isspace():
                has_depth.append(0)
            else:
                has_depth.append(1)
                depth = float(line[31:36])/100
                all_depth.append(depth)
            if (line[147:150]).isspace():
                has_mag.append(0)
            else:
                has_mag.append(1)

            nsta = 0

        # If it is a new phase line with potential S picks:
        elif line[0:2].isalpha():

            nsta += 1

        # If it is shadow line or blank line, skip:
        elif '$' in line[0] or '  ' in line[0:2]:

            all_nsta.append(nsta)

    if summary:
        # Generate report
        print('Within the catalog:')
        print('%d/%d events have nsta >= %d' % (sum(np.array(all_nsta)>=4),len(all_nsta),nsta_thres))
        print('%d/%d events have nobs >= %d' % (sum(np.array(all_nobs)>=4),len(all_nobs),nobs_thres))
        print('%d/%d events have travel time RMSE <= %.2fs' % (sum(np.array(all_rmse)<=0.30),len(all_rmse),rmse_thres))
        print('%d/%d events have a maximum azimuthal gap < %d degrees' % (sum(np.array(all_gap)<180),len(all_gap),gap_thres))
        print('%d/%d events have a horizontal error < %d km' %(sum(np.array(all_eh)<1),len(all_eh),eh_thres))
        print('%d/%d events have a vertical error < %d km' %(sum(np.array(all_ev)<1),len(all_ev),ev_thres))
        print('%d/%d events have magnitude information' %(sum(np.array(has_mag)),len(has_mag)))
        print('%d/%d events have depth information' %(sum(np.array(has_depth)),len(has_depth)))

        num_3P = len([n for n in all_PS if int(n[0:2])>=num_P_thres])
        num_2S = len([n for n in all_PS if int(n[3:5])>=num_S_thres])
        num_3P2S = len([n for n in all_PS if (int(n[0:2])>=num_P_thres and int(n[3:5])>=num_S_thres)])

        print('Summary of phase types:')
        print('%d/%d events have at least %d P picks' % (num_3P,len(all_PS),num_P_thres))
        print('%d/%d events have at least %d S picks' % (num_2S,len(all_PS),num_S_thres))
        print('%d/%d events have at least %dP%dS picks' % (num_3P2S,len(all_PS),num_P_thres,num_S_thres))
    else:
        return all_datetimes, all_PS, all_nobs, all_nsta, all_rmse, all_gap, all_dist2sta, all_depth, all_eh, all_ev, has_mag, has_depth


print('AUGUSTINE VOLCANO\n######')
check_hypoi(augustine_hypoi_file)
print('\n')
print('REDOUBT VOLCANO\n######')
check_hypoi(redoubt_hypoi_file)
#
# plt.figure()
# plt.hist(all_gap,bins=np.arange(60,361,2.5),color='teal',edgecolor='black')
# plt.axvline(x=180,color='red')
# plt.title('Histogram of Redoubt Azimuthal Gaps')
# plt.ylabel('Number of events')
# plt.xlabel('Maximum Azimuthal Gap')
# plt.grid()
# plt.show()
#
# plt.figure()
# plt.hist(all_dist2sta,bins=np.arange(0,10,1),color='orange',edgecolor='black')
# plt.title('Histogram of Augustine Event-Station Distances')
# plt.ylabel('Number of events')
# plt.xlabel('Event-Station Distance')
# plt.grid()
# plt.show()

all_datetimes, all_PS, all_nobs, all_nsta, all_rmse, all_gap, all_dist2sta, all_depth, all_eh, all_ev, has_mag, has_depth = check_hypoi(augustine_hypoi_file,summary=False)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(np.array(all_dist2sta)[np.where(np.array(has_depth)==1)],all_depth,'k.',alpha=0.6)
ax.plot(np.arange(-3,np.max(all_depth)+1)*2,np.arange(-3,np.max(all_depth)+1) - 1.252,'r-')
ax.set_ylabel('Depth (km)')
ax.set_xlabel('Distance to nearest station (km)')
ax.set_title('Augustine Catalog Events')
ax.grid()
ax.invert_yaxis()
# ax.axis('equal')
ax.set_xlim([-0.5,np.max(all_depth)+1])
fig.show()

all_datetimes, all_PS, all_nobs, all_nsta, all_rmse, all_gap, all_dist2sta, all_depth, all_eh, all_ev, has_mag, has_depth = check_hypoi(redoubt_hypoi_file,summary=False)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(np.array(all_dist2sta)[np.where(np.array(has_depth)==1)],all_depth,'k.',alpha=0.6)
ax.plot(np.arange(-3,np.max(all_depth)+1)*2,np.arange(-3,np.max(all_depth)+1) - 3.108,'r-')
ax.set_ylabel('Depth (km)')
ax.set_xlabel('Distance to nearest station (km)')
ax.set_title('Redoubt Catalog Events')
ax.grid()
ax.invert_yaxis()
# ax.axis('equal')
ax.set_xlim([-0.5,np.max(all_depth)+1])
fig.show()