# ### DOWNLOAD COMCAT CATALOG
#
# ## comcat attempts
# from libcomcat.search import search
# from obspy import UTCDateTime
#
# t1=UTCDateTime('2018-01-01 00:00')
# t2=UTCDateTime('2022-04-26 00:00')
# earthquakes = search(starttime=t1.datetime,
# 					 endtime=t2.datetime,
# 					 latitude=57.0509,
#     			 	 longitude=-135.7611,
#     			     maxradiuskm=25,
# 					 reviewstatus='reviewed')
#
# print('Found {:g} earthquakes.'.format(len(earthquakes)))
#
# import requests
# from libcomcat.utils import HEADERS, TIMEOUT
# from obspy.io.quakeml.core import Unpickler
# from obspy import Catalog
#
# catalog = Catalog()
#
# for i,eq in enumerate(earthquakes):
# 	print('Event #{:g}: {}'.format(i+1,eq))
# 	detail = eq.getDetailEvent()
# 	phasedata = detail.getProducts('phase-data')[0]
# 	quakeurl = phasedata.getContentURL('quakeml.xml')
# 	response = requests.get(quakeurl, timeout=TIMEOUT, headers=HEADERS)
# 	data = response.text.encode('utf-8')
# 	unpickler = Unpickler()
# 	catalog += unpickler.loads(data)
#
# catalog.write('/Users/darrentpk/Desktop/GitHub/enhance_catalog/data/comcat/edgecumbe_20180114_20220416.xml','QUAKEML')




# fig, ax = plt.subplots(figsize=(12,5))
# ax.hist(days,bins=bin_days,color='red',edgecolor='black',alpha=0.5)
# ax.plot(cum_days,cum_events,'r-',linewidth=2,label='REDPy Detections')
# ax.plot(cum_days[-1],cum_events[-1],'ro')
# ax2 = ax.twinx()
# ax2.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# ax2.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylabel('Number of REDPy events',color='red',fontsize=14)
# ax2.set_ylim([0,cum_comcat_events[-1]+10])
# ax2.set_ylabel('Cumulative number of ComCat events',color='blue',fontsize=14)
# fig.show()

# # fit trend
# from scipy.optimize import curve_fit
# fit_max_days = (UTCDateTime(2019,7,1)-base_time) / 86400
# def bg(x,m):
#     return m*x
# cum_days_arr = np.array(cum_days)
# cum_events_arr = np.array(cum_events)
# valid_index = np.where(cum_days_arr<fit_max_days)
# mopt,mcov = curve_fit(bg,list(cum_days_arr[valid_index]),list(cum_events_arr[valid_index]))
# cum_events_no_bg = cum_events - cum_days*mopt
# fig, ax = plt.subplots(figsize=(12,5))
# ax.plot(cum_days,cum_events,'r-',linewidth=2,label='REDPy Detections')
# ax.plot(cum_days[-1],cum_events[-1],'ro')
# ax.plot(cum_days,cum_days*mopt,'g-',linewidth=2,label='fit')
# ax2 = ax.twinx()
# ax2.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# ax2.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylim([0,cum_events[-1]+100])
# ax.set_ylabel('Cumulative REDPy Events',color='red',fontsize=14)
# ax2.set_ylim([0,cum_comcat_events[-1]+10])
# ax2.set_ylabel('Cumulative ComCat Events',color='blue',fontsize=14)
# ax.set_title('Cumulative event count (REDPy vs ComCat)',fontsize=17)
# plt.tight_layout()
# fig.show()
# # plot!
# fig, ax = plt.subplots(figsize=(12,5))
# ax.plot(cum_days,cum_events-(cum_days*np.array(1.10)),'r-',linewidth=2,label='REDPy Detections')
# ax.plot(cum_days[-1],cum_events[-1]-(cum_days*mopt)[-1],'ro')
# ax2 = ax.twinx()
# ax2.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# ax2.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# # ax.set_ylim([0,cum_events[-1]+100])
# ax.set_ylabel('Cumulative REDPy Events',color='red',fontsize=14)
# ax2.set_ylim([0,cum_comcat_events[-1]+10])
# ax2.set_ylabel('Cumulative ComCat Events',color='blue',fontsize=14)
# ax.set_title('Cumulative event count (REDPy vs ComCat)',fontsize=17)
# plt.tight_layout()
# fig.show()

# ### Make histogram plot of S-minus-P times
# from toolbox import reader
# catalog = reader('/home/ptan/enhance_catalog/data/comcat/edgecumbe_20180114_20220416.xml')
#
# # Extract SminusP times
# SminusPs = []
# for event in catalog:
#     stachan_list = []
#     for pick in event.picks:
#         sta = pick.waveform_id.station_code
#         chan = pick.waveform_id.channel_code
#         stachan_list.append(sta+'.'+chan)
#     if 'SIT.BHN' in stachan_list and 'SIT.BHZ' in stachan_list:
#         BHN_pick = event.picks[stachan_list.index('SIT.BHN')]
#         BHZ_pick = event.picks[stachan_list.index('SIT.BHZ')]
#         SminusP = BHN_pick.time - BHZ_pick.time
#         SminusPs.append(SminusP)
# print('A total of %d S-minus-P times were extracted from the catalog of %d events' % (len(SminusPs),len(catalog)))
#
# # Create histogram to look at results
# import numpy as np
# import matplotlib.pyplot as plt
# plt.figure()
# bin_lims = np.arange(2,4,0.1)
# plt.hist(SminusPs,bins=bin_lims,edgecolor='k')
# plt.xlabel('Approximate S-minus-P times (s)')
# plt.ylabel('Number of events')
# plt.title('Approximate S-minus-P times for Edgecumbe ComCat Catalog')
# plt.savefig('/home/ptan/enhance_catalog/diffs_comcat.png')
# plt.show()

# Now go into REDPy event list and do computations
import pandas as pd
from obspy import UTCDateTime, read, Stream
redpy_catalog = pd.read_csv('/home/ptan/enhance_catalog/redpy/runs/edgecumbe6/catalog.txt',
                            header=None,names=['cluster','time'],delimiter=' ')
clusters = list(redpy_catalog['cluster'])
times = [UTCDateTime(t) for t in redpy_catalog['time']]

# min_cluster_size = 5  # remove overly small clusters
# clusters_to_remove = range(464,691+1) # active source
# more_clusters_to_remove = [0,14,17,18,19,20,25,26,33,36,47,70,99,100,101,102,115,120,127,156,175,
#                            182,267,285,333,381,384,387,388,391,392,412,438,439]
# clusters_to_keep = [255, 265, 361, 450, 454, 463, 696]

min_cluster_size = 5  # remove overly small clusters
clusters_to_remove = range(879,1104+1) # active source
more_clusters_to_remove = [1,2,3,5,7,11,13,17,21,22,23,25,27,33,39,43,46,48,49,53,55,58,59,62,66,72,73,79,81,82,89,90,92,
                           101,105,106,108,111,114,123,141,146,147,151,152,156,158,160,177,186,189,213,217,252,255,258,
                           259,267,269,274,276,278,281,292,301,305,336,341,348,357,363,370,378,386,409,428,483,485,489,
                           530,532,548,628,629,683,712,722,734,735,746,771,789,803,804,807,809,810,813,814,817,831,855,
                           857,1120,1127,1131]
clusters_to_keep = [129,250,261,419,420,676,677,726,731,768,777,799,827,828,829,830,832,833,834,835,843,844,845,846,
                    848,853,861,866,869,874,878,993,1106,1113,1118,1119,1124,1125,1126,1129,1138,1139,1140,1142,1143,
                    1144,1145,1156,1147]
try_to_remove = [127,137,138,148,149,253,387,415,545]
# try_to_remove = []

index_to_remove = []
for i, cluster in enumerate(clusters):
    count = clusters.count(cluster)
    if cluster in clusters_to_remove or cluster in more_clusters_to_remove or cluster in try_to_remove:
        index_to_remove.append(i)
    elif count <= min_cluster_size:
        index_to_remove.append(i)
for i, cluster in enumerate(clusters):
    if cluster in clusters_to_keep and i in index_to_remove:
        index_to_remove.remove(i)
times = [times[i] for i in range(len(times)) if i not in index_to_remove]
clusters = [clusters[i] for i in range(len(clusters)) if i not in index_to_remove]

# def calc_S_minus_P(time,data_dir='/home/data/edgecumbe/',plot=False,verbose=False,prewin=3,postwin=5):
#     # First read daylong BHN adnd BHZ stream
#     import glob
#     from matplotlib.dates import num2date
#     year = time.year
#     julday = time.julday
#     julday_str = '%03d' % julday
#     BHN_filepaths = glob.glob(data_dir + 'SIT.BHN.' + str(year) + ':' + julday_str + ':*')
#     BHN_stream = Stream()
#     for filepath in BHN_filepaths:
#         contributing_stream = read(filepath)
#         BHN_stream += contributing_stream
#     BHZ_filepaths = glob.glob(data_dir + 'SIT.BHZ.' + str(year) + ':' + julday_str + ':*')
#     BHZ_stream = Stream()
#     for filepath in BHZ_filepaths:
#         contributing_stream = read(filepath)
#         BHZ_stream += contributing_stream
#     # Trim around stream and plot
#     BHN_stream.merge()
#     BHZ_stream.merge()
#     st = BHZ_stream + BHN_stream
#     st.trim(starttime=time-prewin,endtime=time+postwin)
#     st = st.split()
#     st.detrend()
#     st.taper(0.05,type='cosine')
#     st.filter(type='bandpass',freqmin=1,freqmax=10)
#     # place np.nan in difference in case it fails
#     difference = np.nan
#     # do a check then run differecne calculator
#     if len(st) == 2:
#         # extract stream values
#         BHZ_t = st[0].times("matplotlib")
#         BHZ_d = st[0].data
#         BHN_t = st[1].times("matplotlib")
#         BHN_d = st[1].data
#         # Try to find P and S arrival with np.max
#         p_arr_ind = np.argmax(np.abs(BHZ_d[:(len(BHZ_d)//2)]))
#         s_arr_ind = np.argmax(np.abs(BHN_d[(len(BHN_d)//2):])) + len(BHN_d)//2
#         difference = UTCDateTime(num2date(BHN_t[s_arr_ind])) - UTCDateTime(num2date(BHZ_t[p_arr_ind]))
#         if verbose:
#             print('Time = %s gave a S-minus-P difference of %.3fs' % (str(time),difference))
#         # Do plotting
#         if plot:
#             fig, ax = plt.subplots(nrows=2,figsize=(7,4))
#             ax[0].plot(BHZ_t, BHZ_d, "k-",label=st[0].id)
#             ax[0].plot(BHZ_t[p_arr_ind], BHZ_d[p_arr_ind], 'r*', markersize=10)
#             ax[0].axvline(time.matplotlib_date,color='r', linestyle='--', alpha=0.5)
#             ax[0].legend(loc=2)
#             ax[0].xaxis_date()
#             ax[0].set_title('BHN peak - BHZ peak = %.3fs' % difference)
#             ax[1].plot(BHN_t, BHN_d, "k-",label=st[1].id)
#             ax[1].plot(BHN_t[s_arr_ind], BHN_d[s_arr_ind], "b*", markersize=10)
#             ax[1].xaxis_date()
#             ax[1].legend(loc=2)
#             fig.autofmt_xdate()
#             fig.show()
#     else:
#         if verbose:
#             print('Time = %s does not trim and merge to 2 traces' % str(time))
#     return difference
#
# diffs = []
# for time in times:
#     diff = calc_S_minus_P(time,data_dir='/home/data/edgecumbe/',plot=False,verbose=False)
#     diffs.append(diff)
#
# plt.figure()
# bin_lims = np.arange(0,5,0.1)
# plt.hist(diffs,bins=bin_lims,color='#d62728',edgecolor='k',label='REDPy')
# # plt.hist(SminusPs,bins=bin_lims,color='#1f77b4',edgecolor='k',label='ComCat')
# # plt.legend()
# plt.axvline(x=2.3,linestyle='--')
# plt.axvline(x=3.3,linestyle='--')
# plt.xlabel('Approximate S-minus-P times (s)')
# plt.ylabel('Number of detections')
# plt.title('Approximate S-minus-P times for Edgecumbe Detections')
# plt.savefig('/home/ptan/enhance_catalog/diffs_hist.png')
# plt.show()
#
# fig, ax = plt.subplots(figsize=(7,3))
# ax.scatter([t.matplotlib_date for t in times],diffs,s=2,color='r',marker='.')
# ax.xaxis_date()
# ax.set_xlabel('Time of Detection (UTC)')
# ax.set_ylabel('Approx S-minus-P difference (s)')
# fig.autofmt_xdate()
# plt.tight_layout()
# fig.savefig('/home/ptan/enhance_catalog/diffs_scatter.png')
# fig.show()

# diff_array = np.array(diffs)
# anomaly_ind = np.where((diff_array>2.5) & (diff_array<2.7))
# anomaly_times = np.array(times)[anomaly_ind]
# len(anomaly_times)
#
# for time in anomaly_times[60:80]:
#     diff = calc_S_minus_P(time,data_dir='/home/data/edgecumbe/',plot=True,verbose=False)

clean_catalog_txt = open('/home/ptan/enhance_catalog/edgecumbe_clean_catalog_vFINAL.txt', 'w')
# format = '%d %s %.3f\n'
format = '%d %s\n'
for i in range(len(times)):
    # clean_catalog_txt.write(format % (clusters[i],times[i],diffs[i]))
    clean_catalog_txt.write(format % (clusters[i], times[i]))
clean_catalog_txt.close()

redpy_times_txt = open('/home/ptan/enhance_catalog/edgecumbe_redpy_times_vFINAL.txt', 'w')
format = '%s\n'
for i in range(len(times)):
    redpy_times_txt.write(format % times[i])
redpy_times_txt.close()

# comcat_times_txt = open('/home/ptan/enhance_catalog/edgecumbe_comcat_times_v2.txt', 'w')
# format = '%s\n'
# for i in range(len(catalog)):
#     comcat_times_txt.write(format % catalog[i].origins[0].time)
# comcat_times_txt.close()





import numpy as np
import pandas as pd
from obspy import UTCDateTime
import matplotlib.pyplot as plt

# read redpy detections
# textfile = '/home/ptan/enhance_catalog/edgecumbe_clean_catalog_vFINAL.txt'
# df = pd.read_csv(textfile,names=['cluster','time','diff'],delimiter=' ')
textfile = '/home/ptan/enhance_catalog/edgecumbe_redpy_times.txt'
df = pd.read_csv(textfile,names=['time'],delimiter=' ')
times = [UTCDateTime(t) for t in list(df['time'])]
times = list(np.sort(times))
# read comcat
comcatcsv = '/home/ptan/enhance_catalog/data/comcat/edgecumbe_query.csv'
df2 = pd.read_csv(comcatcsv)
comcat_times = [UTCDateTime(t) for t in list(df2.time)]
comcat_times = list(np.sort(comcat_times))
# read wech's detections
wechtxt = '/home/ptan/enhance_catalog/wech.txt'
df3 = pd.read_csv(wechtxt,delimiter='\t')
wech_times = [UTCDateTime(t) for t in list(df3.time) if UTCDateTime(t) < UTCDateTime(2022,4,16)]
# wechtxt = '/home/ptan/enhance_catalog/edgecumbe_matchedfilter_times.txt'
# df3 = pd.read_csv(wechtxt,names=['time'],delimiter=' ')
# wech_times = [UTCDateTime(t) for t in list(df3['time']) if UTCDateTime(t) < UTCDateTime(2022,4,16)]
wech_times = list(np.sort(wech_times))
# sort out bin limits
bin_lims = []
for year in range(2014,2022+1):
    for month in range(1,13):
        if year==2014 and month<10:
            continue
        elif year==2022 and month>5:
            continue
        bin_lims.append(UTCDateTime(year,month,1))
# sort out tick limits
tick_lims = []
for year in range(2014,2022+1):
    for month in range(1,13,3):
        if year==2014 and month<10:
            continue
        if year==2022 and month>5:
            continue
        tick_lims.append(UTCDateTime(year,month,1))
# convert into days
base_time = UTCDateTime(2014,10,17)
days = [(time-base_time)/86400 for time in times]
comcat_days = [(time-base_time)/86400 for time in comcat_times]
wech_days = [(time-base_time)/86400 for time in wech_times]
bin_days = [(bin_lim-base_time)/86400 for bin_lim in bin_lims]
tick_days = [(tick_lim-base_time)/86400 for tick_lim in tick_lims]
tick_labels = [tick_lim.strftime('%b \'%y') for tick_lim in tick_lims]
# general cumulative counts for redpy and for comcat
cum_days = [0, days[0], *days]
cum_events = [0, *range(0,len(days)+1)]
cum_comcat_days = [comcat_days[0], *comcat_days]
cum_comcat_events = range(0,len(comcat_days)+1)
cum_wech_days = [wech_days[0], *wech_days]
cum_wech_events = range(0,len(wech_days)+1)

# plot REDPY + OBSPY-CORRSCAN
fig, ax = plt.subplots(figsize=(12,5))
h1 = ax.hist(days,bins=bin_days,color='maroon',edgecolor='black',alpha=0.6,label='REDPy detections')
h2 = ax.hist(wech_days,bins=bin_days,color='teal',edgecolor='black',alpha=0.4,label='ObsPy correlation detections')
ax.legend(fontsize=14,loc='upper left')
ax2 = ax.twinx()
ax2.plot(cum_days,cum_events,'r-',linewidth=2,label='Cumulative REDPy detections')
ax2.plot(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='Cumulative ObsPy correlation detections')
ax2.plot(cum_days[-1],cum_events[-1],'ro')
ax2.plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
ax2.legend(fontsize=14,bbox_to_anchor=(0.472, 0.83))
# ax2.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# ax2.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# ax2.legend(fontsize=14)
# tidy axes
ax.set_xticks(tick_days)
ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
ax.set_ylabel('Number of detections per month',color='black',fontsize=14)
ax2.set_ylim([0,cum_events[-1]+100])
ax2.set_ylabel('Cumulative number of detections',color='black',fontsize=14)
ax.set_title('Edgecumbe Detections (Oct 17, 2014 - Apr 16, 2022)',fontsize=17)
plt.tight_layout()
# fig.savefig('/home/ptan/enhance_catalog/histogram_and_cumulative_count_dual.pdf')
fig.show()

# # Plot REDPY
# fig, ax = plt.subplots(figsize=(12,5))
# ax.hist(days,bins=bin_days,color='teal',edgecolor='black')
# ax2 = ax.twinx()
# ax2.plot(cum_days,cum_events,'r-',linewidth=2,label='Automated Detections')
# ax2.plot(cum_days[-1],cum_events[-1],'ro')
# # ax2.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# # ax2.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # ax2.legend(fontsize=14)
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylabel('Number of REDPy events per month',color='teal',fontsize=14)
# ax2.set_ylim([0,cum_events[-1]+100])
# ax2.set_ylabel('Cumulative number of REDPy events',color='red',fontsize=14)
# ax.set_title('Edgecumbe REDPy Detections (Jan 2018 - Apr 2022)',fontsize=17)
# plt.tight_layout()
# # fig.savefig('/home/ptan/enhance_catalog/histogram_and_cumulative_count.pdf')
# fig.show()
#
#
# # plot!
fig, ax = plt.subplots(figsize=(12,5))
ax.plot(cum_days,cum_events,'-',color='red',linewidth=2,label='REDPy Detections')
ax.plot(cum_days[-1],cum_events[-1],'o',color='red')
ax.plot(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='ObsPy Correlation Detections')
ax.plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
ax.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Located Events')
ax.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# tidy axes
ax.set_xticks(tick_days)
ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
ax.set_ylim([0,cum_events[-1]+100])
# ax.set_xlim((0,(UTCDateTime(2022,4,1)-base_time)/86400))
# ax.set_ylim([0,2200])
ax.legend(fontsize=15,loc='upper left')
ax.set_ylabel('Cumulative number of events',fontsize=14)
ax.set_title('Cumulative event count comparison (REDPy vs ObsPy Correlation vs ComCat)',fontsize=17)
plt.tight_layout()
# fig.savefig('/home/ptan/enhance_catalog/cumulative_count_comparison_nolog.pdf')
fig.show()
#
# # second plot
# # plot!
# fig, axs = plt.subplots(2,1,figsize=(12,8))
# axs[0].plot(cum_days,cum_events,'-',color='maroon',linewidth=2,label='REDPy Detections')
# axs[0].plot(cum_days[-1],cum_events[-1],'o',color='maroon')
# axs[0].plot(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='EQcorrscan Detections')
# axs[0].plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
# axs[0].plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# axs[0].plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # tidy axes
# axs[0].set_xticks(tick_days)
# axs[0].set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# axs[0].set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# axs[0].set_ylim([0,cum_events[-1]+100])
# axs[0].legend(fontsize=15)
# axs[0].set_ylabel('Cumulative number of events',fontsize=14)
# axs[0].set_title('Cumulative event count comparison (REDPy vs EQcorrscan vs ComCat)',fontsize=17)
# axs[1].semilogy(cum_days,cum_events,'-',color='maroon',linewidth=2,label='REDPy Detections')
# axs[1].plot(cum_days[-1],cum_events[-1],'o',color='maroon')
# axs[1].semilogy(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='EQcorrscan Detections')
# axs[1].plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
# axs[1].semilogy(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# axs[1].plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# # tidy axes
# axs[1].set_xticks(tick_days)
# axs[1].set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# axs[1].set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# axs[1].set_ylim([0,10000])
# axs[1].set_ylabel('Cumulative count in log scale',fontsize=14)
# plt.tight_layout()
# # fig.savefig('/home/ptan/enhance_catalog/cumulative_count_comparison_dual.pdf')
# fig.show()
#
# # second plot
# # plot!
# fig, ax = plt.subplots(figsize=(12,5))
# l1 = ax.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Located Events')
# ax.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# ax2 = ax.twinx()
# l2 = ax2.plot(cum_days,cum_events,'-',color='maroon',linewidth=2,label='REDPy Detections')
# ax2.plot(cum_days[-1],cum_events[-1],'o',color='maroon')
# l3 = ax2.plot(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='EQcorrscan Detections')
# ax2.plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
# objs = l2+l3+l1
# labs = [obj.get_label() for obj in objs]
# ax.legend(objs, labs, loc='upper left',fontsize=15)
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylim([0,cum_comcat_events[-1]+200])
# ax.set_ylabel('Cumulative number of ComCat events',color='blue',fontsize=14)
# ax2.set_ylim([0,cum_events[-1]+100])
# ax2.set_ylabel('Cumulative number of automated detections',color='black',fontsize=14)
# multicolor_ylabel(ax2,('Cumulative number of','REDPy','/','EQcorrscan','detections'),('k','maroon','k','teal','k'),axis='y')
# ax.set_title('Cumulative event count comparison (REDPy & EQcorrscan vs ComCat)',fontsize=17)
# plt.tight_layout()
# # fig.savefig('/home/ptan/enhance_catalog/cumulative_count_comparison.pdf')
# fig.show()
#
# # second plot
# # plot!
# fig, ax = plt.subplots(figsize=(12,5))
# ax.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Events')
# ax.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# ax2 = ax.twinx()
# ax2.plot(cum_days,cum_events,'r-',linewidth=2,label='REDPy Detections')
# ax2.plot(cum_days[-1],cum_events[-1],'ro')
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylim([0,cum_comcat_events[-1]+10])
# ax.set_ylabel('Cumulative number of ComCat Events',color='blue',fontsize=14)
# ax2.set_ylim([0,cum_events[-1]+100])
# ax2.set_ylabel('Cumulative number of REDPy Events',color='red',fontsize=14)
# ax.set_title('Cumulative event count comparison (REDPy vs ComCat)',fontsize=17)
# plt.tight_layout()
# # fig.savefig('/home/ptan/enhance_catalog/cumulative_count_comparison.pdf')
# fig.show()
#
# cumu_redpy_txt = open('/home/ptan/enhance_catalog/edgecumbe_redpy_cumu_v3.txt', 'w')
# format = '%.3f %d\n'
# for i in range(len(cum_days)):
#     cumu_redpy_txt.write(format % (cum_days[i],cum_events[i]))
# cumu_redpy_txt.close()
#
# cumu_eqcs_txt = open('/home/ptan/enhance_catalog/edgecumbe_eqcs_cumu.txt', 'w')
# format = '%.3f %d\n'
# for i in range(len(cum_wech_days)):
#     cumu_eqcs_txt.write(format % (cum_wech_days[i],cum_wech_events[i]))
# cumu_eqcs_txt.close()
#
# cumu_comcat_txt = open('/home/ptan/enhance_catalog/edgecumbe_comcat_cumu.txt', 'w')
# format = '%.3f %d\n'
# for i in range(len(cum_comcat_days)):
#     cumu_comcat_txt.write(format % (cum_comcat_days[i],cum_comcat_events[i]))
# cumu_comcat_txt.close()
#
# redpy_times_txt = open('/home/ptan/enhance_catalog/edgecumbe_redpy_times.txt', 'w')
# format = '%s\n'
# for i in range(len(times)):
#     redpy_times_txt.write(format % times[i])
# redpy_times_txt.close()
#
matchedfilter_times_txt = open('/home/ptan/enhance_catalog/edgecumbe_matchedfilter_times.txt', 'w')
format = '%s\n'
for i in range(len(wech_times)):
    matchedfilter_times_txt.write(format % wech_times[i])
matchedfilter_times_txt.close()
#
#
#
#
#
#
# def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
#     """this function creates axes labels with multiple colors
#     ax specifies the axes object where the labels should be drawn
#     list_of_strings is a list of all of the text items
#     list_if_colors is a corresponding list of colors for the strings
#     axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
#     from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
#
#     # x-axis label
#     if axis=='x' or axis=='both':
#         boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
#                     for text,color in zip(list_of_strings,list_of_colors) ]
#         xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
#         anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
#                                           bbox_transform=ax.transAxes, borderpad=0.)
#         ax.add_artist(anchored_xbox)
#
#     # y-axis label
#     if axis=='y' or axis=='both':
#         boxes = [TextArea(text, textprops=dict(color=color, ha='right',va='bottom',rotation=90,**kw))
#                      for text,color in zip(list_of_strings[::-1],list_of_colors) ]
#         ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
#         anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(1.05, -0.2),
#                                           bbox_transform=ax.transAxes, borderpad=0.)
#         ax.add_artist(anchored_ybox)
#
#
# # plot!
# fig, ax = plt.subplots(figsize=(12,5))
# l1 = ax.plot(cum_comcat_days,cum_comcat_events,'b-',linewidth=2,label='ComCat Located Events')
# ax.plot(cum_comcat_days[-1],cum_comcat_events[-1],'bo')
# ax2 = ax.twinx()
# l2 = ax2.plot(cum_days,cum_events,'-',color='maroon',linewidth=2,label='REDPy Detections')
# ax2.plot(cum_days[-1],cum_events[-1],'o',color='maroon')
# l3 = ax2.plot(cum_wech_days,cum_wech_events,'-',color='teal',linewidth=2,label='EQcorrscan Detections')
# ax2.plot(cum_wech_days[-1],cum_wech_events[-1],'o',color='teal')
# objs = l2+l3+l1
# labs = [obj.get_label() for obj in objs]
# ax.legend(objs, labs, loc='upper left',fontsize=15)
# # tidy axes
# ax.set_xticks(tick_days)
# ax.set_xticklabels(tick_labels,fontsize=13,rotation=30,ha='right')
# ax.set_xlim((0,(UTCDateTime(2022,5,1)-base_time)/86400))
# ax.set_ylim([0,cum_comcat_events[-1]+200])
# ax.set_ylabel('Cumulative no. of ComCat events',color='blue',fontsize=14)
# ax2.set_ylim([0,cum_events[-1]+100])
# # ax2.set_ylabel('Cumulative number of automated detections',color='black',fontsize=14)
# multicolor_ylabel(ax2,('Cumulative no. of','REDPy','/','EQcorrscan','events'),('k','teal','k','maroon','k'),axis='y',fontsize=14)
# ax.set_title('Cumulative event count comparison (REDPy & EQcorrscan vs ComCat)',fontsize=17)
# plt.tight_layout()
# # fig.savefig('/home/ptan/enhance_catalog/cumulative_count_comparison_triple.pdf')
# fig.show()