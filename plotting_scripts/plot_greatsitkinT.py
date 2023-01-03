import numpy as np
from toolbox import reader
from matplotlib import pyplot as plt
from obspy import Catalog, UTCDateTime
from pyproj import Geod
from obspy import UTCDateTime
from pandas import read_csv

main_dir = '/home/ptan/enhance_catalog/'
output_dir = main_dir + 'output/greatsitkin/'
# PEC_events_filepath = output_dir + 'relocate_catalog/PEC_events_reFI.xml'
party_catalog_filepath = output_dir + 'scan_data/party_catalog.xml'
relocated_catalog_filepath = output_dir + 'relocate_catalog/relocated_catalog.xml'
# PEC_events = reader(PEC_events_filepath)
raw_party_catalog = reader(party_catalog_filepath)
raw_relocated_catalog = reader(relocated_catalog_filepath)
gs_lat = 52.076
gs_lon = -176.13
volcano_coord = (gs_lat,gs_lon)

# # Load explosion times
# explosion_table = read_csv(main_dir + 'data/supporting/aug_explosion_list.csv')
# explosion_label = []
# explosion_times = []
# for i in range(len(explosion_table)):
#     explosion_label.append(str(explosion_table.event_number[i]))
#     explosion_times.append(UTCDateTime(explosion_table.year[i],explosion_table.month[i],explosion_table.day[i],explosion_table.hh[i],explosion_table.mm[i]))

###################

start_time = UTCDateTime(2021, 5, 12)
end_time = UTCDateTime(2021, 8, 12)
tick_precision = '7day' # '1hour', '6hour','day','7day','month'
event_count_style = 'linear' # 'log', 'linear'
plot_explosions = False

if tick_precision == '1hour':
    tick_width = 3600
    num_ticks = int((end_time-start_time)/tick_width + 1)
    tick_times = [start_time + num_tick*tick_width for num_tick in range(num_ticks)]
elif tick_precision == '6hour':
    tick_width = 21600
    num_ticks = int((end_time-start_time)/tick_width + 1)
    tick_times = [start_time + num_tick*tick_width for num_tick in range(num_ticks)]
elif tick_precision == 'day':
    tick_width = 86400
    num_ticks = int((end_time-start_time)/tick_width + 1)
    tick_times = [start_time + num_tick*tick_width for num_tick in range(num_ticks)]
elif tick_precision == '7day':
    tick_width = 86400*7
    num_ticks = int((end_time-start_time)/tick_width + 1)
    tick_times = [start_time + num_tick*tick_width for num_tick in range(num_ticks)]
elif tick_precision == 'month':
    num_ticks = int(np.floor((end_time-start_time)/(86400*30)) + 1)
    tick_times = []
    current_month = start_time.month
    for num_tick in range(num_ticks):
        if current_month + num_tick <= 12:
            tick_times.append(start_time.replace(month=start_time.month+num_tick))
        else:
            tick_times.append(start_time.replace(year=start_time.year+1, month=start_time.month+num_tick-12))
else:
    raise ValueError('wrong tick precision setting')

party_catalog = Catalog([event for event in raw_party_catalog if UTCDateTime(event.resource_id.id.split('_')[-1]) > start_time])
party_catalog = Catalog([event for event in party_catalog if UTCDateTime(event.resource_id.id.split('_')[-1]) < end_time])
relocated_catalog = Catalog([event for event in raw_relocated_catalog if event.origins[0].time > start_time])
relocated_catalog = Catalog([event for event in relocated_catalog if event.origins[0].time < end_time])

catalog_mag = Catalog([event for event in relocated_catalog if event.magnitudes != []])
catalog_fi = Catalog([event for event in party_catalog if event.comments[-1].text[0:2] == 'FI'])
catalog_fi = Catalog([event for event in catalog_fi if event.comments[-1].text.split('=')[1] != 'None'])

avcc = [round((float(event.comments[2].text.split('=')[1]) / (len(event.comments[3].text.split(')'))-1)),4) for event in party_catalog]
mag = [event.magnitudes[0].mag for event in catalog_mag]
depth = [event.origins[0].depth / 1000 for event in relocated_catalog]
fi = [float(event.comments[-1].text.split('=')[1]) for event in catalog_fi]
station_coords = [(event.origins[0].latitude,event.origins[0].longitude) for event in relocated_catalog]
dist = []
faz = []
baz = []
geodesic = Geod(ellps='WGS84')
for event in relocated_catalog:
    az12, az21 ,dist12 = geodesic.inv(gs_lon, gs_lat, event.origins[0].longitude, event.origins[0].latitude)
    if az12 < 0:
        faz.append(az12 + 360)
    else:
        faz.append(az12)
    baz.append(az21)
    dist.append(dist12/1000)

all_times = [event.origins[0].time for event in relocated_catalog]
days = [(time-start_time)/86400 for time in all_times]
days_cc = [((UTCDateTime(event.resource_id.id.split('_')[-1])-start_time)/86400) for event in party_catalog]
days_mag = [((event.origins[0].time-start_time)/86400) for event in catalog_mag]
days_fi = [((UTCDateTime(event.resource_id.id.split('_')[-1])-start_time)/86400) for event in catalog_fi]
if plot_explosions:
    days_exp = [(explosion_time - start_time) / 86400 for explosion_time in explosion_times]

num_days = int((tick_times[-1] - tick_times[0])/86400)

days_evcount = []
evcounts = []
all_times = [UTCDateTime(event.resource_id.id.split('_')[-1]) for event in party_catalog]
if tick_precision == 'month' or tick_precision == '7day':
    for k in range(num_days):
        count_start = start_time + k*86400
        count_end = start_time + (k+1)*86400
        days_evcount.append(k+0.5)
        evcount = len([time for time in all_times if (time>count_start) and (time<count_end)])
        evcounts.append(evcount)
else:
    for k in range(num_ticks):
        count_start = start_time + k*tick_width
        count_end = start_time + (k+1)*tick_width
        days_evcount.append((k+0.5)*tick_width/86400)
        evcount = len([time for time in all_times if (time>count_start) and (time<count_end)])
        evcounts.append(evcount)

fig, axs = plt.subplots(7, 1, sharex=True,figsize=(9,10))
fig.subplots_adjust(hspace = 0.25)
if event_count_style == 'log':
    axs[0].semilogy(days_evcount,evcounts,'k-')
    axs[0].set_ylim([1,10000])
    axs[0].set_yticks([1,10,100,1000,10000])
elif event_count_style == 'linear':
    axs[0].plot(days_evcount,evcounts,'k-')
axs[0].tick_params(axis='y',labelsize=13)
if tick_precision == 'month' or tick_precision == 'day' or tick_precision == '7day':
    axs[0].set_ylabel('Events/Day',fontsize=13)
else:
    axs[0].set_ylabel('Events/%s' % tick_precision, fontsize=13)
title = 'Davidof Seismicity Trends, %s to %s' % (tick_times[0].strftime('%d %b %Y'),tick_times[-1].strftime('%d %b %Y'))
axs[0].set_title(title, fontsize=18)
if plot_explosions:
    for ax_number in range(1,7):
        for day_exp in days_exp:
            axs[ax_number].axvline(x=day_exp, color='dimgray', alpha=0.5, zorder=0)
axs[1].scatter(days_cc,avcc,s=18,c='firebrick',marker='.',alpha=0.4)
axs[1].set_ylim([0.7,1.0])
axs[1].set_yticks([0.7,0.8,0.9,1.0])
axs[1].set_ylabel('Avg CC',fontsize=13)
axs[1].tick_params(axis='y',labelsize=13)
axs[2].scatter(days_fi,fi,s=18,c='red',marker='.',alpha=0.4)
axs[2].set_ylim([-1,1.1])
axs[2].set_yticks([-1.0,-0.5,0,0.5,1.0])
axs[2].set_ylabel('FI',fontsize=13)
axs[2].tick_params(axis='y',labelsize=13)
axs[3].scatter(days_mag,mag,s=18,c='black',marker='.',alpha=0.4)
axs[3].set_ylim([-1,5])
axs[3].set_yticks([-1,1,3,5])
axs[3].set_ylabel('Magnitude',fontsize=13)
axs[3].tick_params(axis='y',labelsize=13)
axs[4].scatter(days,depth,s=18,c='darkblue',marker='.',alpha=0.4)
axs[4].set_ylim([-5,15])
axs[4].set_yticks([-5,0,5,10,15])
axs[4].invert_yaxis()
axs[4].set_ylabel('Depth (km)',fontsize=13)
axs[4].tick_params(axis='y',labelsize=13)
axs[5].scatter(days,dist,s=18,c='darkolivegreen',marker='.',alpha=0.4)
axs[5].set_ylim([0,21])
axs[5].set_yticks([0,7,14,21])
axs[5].set_ylabel('Dist (km)',fontsize=13)
axs[5].tick_params(axis='y',labelsize=13)
axs[6].scatter(days,faz,s=18,c='rebeccapurple',marker='.',alpha=0.4)
axs[6].set_ylim([0,360])
axs[6].set_yticks([0,90,180,270,360])
axs[6].set_ylabel('Azimuth (°)',fontsize=13)
axs[6].tick_params(axis='y',labelsize=13)
axs[6].set_xlim([0,num_days])
axs[6].set_xticks((np.array(tick_times) - start_time) / 86400)
if tick_precision == '1hour':
    axs[6].set_xticklabels([tick_time.strftime('%m/%d %H:00') for tick_time in tick_times], fontsize=13, rotation=30, ha='right')
elif tick_precision == '6hour':
    axs[6].set_xticklabels([tick_time.strftime('%m/%d %H:00') for tick_time in tick_times], fontsize=13, rotation=30, ha='right')
elif tick_precision == 'day':
    axs[6].set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times], fontsize=13, rotation=30, ha='right')
elif tick_precision == '7day':
    axs[6].set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times], fontsize=13, rotation=30, ha='right')
elif tick_precision == 'month':
    axs[6].set_xticklabels([tick_time.strftime('%b \'%y') for tick_time in tick_times], fontsize=13, rotation=30, ha='right')
else:
    raise ValueError('wrong tick precision setting')
save_path = output_dir + 'seismic_trends_' + start_time.strftime('%Y%m%d') + '_' + end_time.strftime('%Y%m%d')+ '.png'
plt.savefig(save_path)
fig.show()
