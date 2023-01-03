#%% PLOT REDOUBT TEMPORAL

# Import all dependencies
from obspy import UTCDateTime
import numpy as np
from toolbox import reader, get_color_codes
from pandas import read_csv
import matplotlib.pyplot as plt

# Define general paths
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
catalog_filepath = main_dir + 'output/redoubt5/scan_data/party_FI.xml'
catalog = reader(catalog_filepath)
relocated_catalog_filepath = main_dir + 'output/redoubt5/scan_data/relocatable_catalog_FImag.xml'
relocated_catalog = reader(relocated_catalog_filepath)
for event in relocated_catalog:
    if event.magnitudes == []:
        relocated_catalog.events.remove(event)
original_filepath = main_dir + 'output/redoubt5/scan_data/PEC_FI.xml'
original = reader(original_filepath)
size_by_mag = True
add_labels = True

# event list
utctimes = np.array([UTCDateTime(event.resource_id.id.split('_')[-1]) for event in catalog])
FI_values = []
for event in catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        FI_values.append(np.nan)
    else:
        FI_values.append(float(event.comments[-1].text.split('=')[1]))
FI_values = np.array(FI_values)

# relocated catalog
relocated_utctimes = np.array([event.origins[0].time for event in relocated_catalog])
relocated_FI_values = []
for event in relocated_catalog:
    if event.comments[-1].text.split('=')[1] == 'None':
        relocated_FI_values.append(np.nan)
    else:
        relocated_FI_values.append(float(event.comments[-1].text.split('=')[1]))
relocated_FI_values = np.array(relocated_FI_values)

# remove overlapping
remove_index = []
for i, relocated_utctime in enumerate(relocated_utctimes):
    remove_index.append(np.argmin(np.abs(utctimes-relocated_utctime)))
    print(i)
sorted_remove_index = list(np.sort(remove_index))
utctimes = list(utctimes)
FI_values = list(FI_values)
for i in reversed(sorted_remove_index):
    utctimes.pop(i)
    FI_values.pop(i)
utctimes = np.array(utctimes)
FI_values = np.array(FI_values)

original_utctimes = np.array([event.origins[0].time for event in original])
original_FI_values = []
for event in original:
    if event.comments[-1].text.split('=')[1] == 'None':
        original_FI_values.append(np.nan)
    else:
        original_FI_values.append(float(event.comments[-1].text.split('=')[1]))
original_FI_values = np.array(original_FI_values)

def calibrate_sizes(magnitudes):
    # construct reference sizes
    ref_magnitudes = np.array(range(0, 401))
    min_mag_marker = 1
    ref_sizes = []
    for magnitude in list(ref_magnitudes):
        normalized_magnitude = (magnitude - np.min(min_mag_marker)) / np.ptp(
            ref_magnitudes[ref_magnitudes > min_mag_marker])
        ref_size = 50 + normalized_magnitude * (5*160)
        ref_sizes.append(ref_size)
    ref_magnitudes = list(ref_magnitudes)
    # now calculate sizes
    sizes = []
    for magnitude in magnitudes:
        if magnitude < min_mag_marker:
            size = ref_sizes[0]
        else:
            size = ref_sizes[ref_magnitudes.index(int(magnitude * 100))]
        sizes.append(size)
    sizes = np.array(sizes)
    return sizes

if size_by_mag:
    relocated_magnitudes = np.array([event.magnitudes[0].mag for event in relocated_catalog])
    relocated_sizes = calibrate_sizes(relocated_magnitudes)
    original_magnitudes = np.array([event.magnitudes[0].mag for event in original])
    original_sizes = calibrate_sizes(original_magnitudes)

# Load color code labels
color_codes, color_code_times = get_color_codes('redoubt')

# Load explosion times
explosion_table = read_csv(main_dir + 'data/supporting/rdt_explosion_list.csv')
explosion_label = []
explosion_times = []
for i in range(len(explosion_table)):
    explosion_label.append(str(explosion_table.event_number[i]))
    explosion_times.append(UTCDateTime(explosion_table.year[i],explosion_table.month[i],explosion_table.day[i],explosion_table.hh[i],explosion_table.mm[i]))

# Prepare cumulative count stuff
volc_PEC = original
volc_PARTY = catalog
volc_PEC_times = [e.origins[0].time for e in volc_PEC]
volc_PEC_times = list(np.sort(volc_PEC_times))
volc_PARTY_times = [UTCDateTime(e.resource_id.id.split('_')[-1]) for e in volc_PARTY]
volc_PARTY_times = list(np.sort(volc_PARTY_times))

# plot it

# prepare times

# Initialize figure
fig = plt.figure(figsize=(44, 25))
# Cumulative count plot
ax0 = plt.subplot2grid((8, 3), (0, 0), colspan=3, rowspan=2)
volc_tick_lims = [UTCDateTime(2008, 4, 1),
                  UTCDateTime(2008, 5, 1), UTCDateTime(2008, 6, 1),
                  UTCDateTime(2008, 7, 1), UTCDateTime(2008, 8, 1),
                  UTCDateTime(2008, 9, 1), UTCDateTime(2008, 10, 1),
                  UTCDateTime(2008, 11, 1), UTCDateTime(2008, 12, 1),
                  UTCDateTime(2009, 1, 1), UTCDateTime(2009, 2, 1),
                  UTCDateTime(2009, 3, 1), UTCDateTime(2009, 4, 1),
                  UTCDateTime(2009, 5, 1), UTCDateTime(2009, 6, 1),
                  UTCDateTime(2009, 7, 1), UTCDateTime(2009, 8, 1),
                  UTCDateTime(2009, 9, 1)]
volc_bin_lims = [UTCDateTime(2008, 4, 1)]
while volc_bin_lims[-1] < UTCDateTime(2009, 9, 1):
    volc_bin_lims.append(volc_bin_lims[-1] + (86400*7))
volc_base_time = UTCDateTime(2008,4,1)
volc_PEC_days = [(time-volc_base_time)/86400 for time in volc_PEC_times]
volc_PARTY_days = [(time-volc_base_time)/86400 for time in volc_PARTY_times]
volc_tick_days = [(volc_tick_lim-volc_base_time)/86400 for volc_tick_lim in volc_tick_lims]
volc_tick_labels = [volc_tick_lim.strftime('%b \'%y') for volc_tick_lim in volc_tick_lims]
volc_bin_days = [(bin_lim-volc_base_time)/86400 for bin_lim in volc_bin_lims]
cum_volc_PEC_days = [0, volc_PEC_days[0], *volc_PEC_days]
cum_volc_PEC_events = [0, *range(0,len(volc_PEC_days)+1)]
cum_volc_PARTY_days = [0, volc_PARTY_days[0], *volc_PARTY_days]
cum_volc_PARTY_events = [0, *range(0,len(volc_PARTY_days)+1)]
ax0.hist(volc_PARTY_days,bins=volc_bin_days,color='#505050',edgecolor='black',alpha=0.425,label='Redoubt EQcorrscan')
ax0.hist(volc_PEC_days,bins=volc_bin_days,color='darkred',edgecolor='black',alpha=0.425,label='Redoubt Original Catalog')
# ax0.legend(fontsize=35,loc='upper left')
ax0.legend(fontsize=35,loc='upper left', bbox_to_anchor=(0.021,1.05))
ax0.text(2,1000,'a)',fontsize=58,weight='bold')
ax0.set_xticks(volc_tick_days)
ax0.set_xticklabels(volc_tick_labels,fontsize=35,rotation=30,ha='right')
ax0.set_xlim((0,(UTCDateTime(2009,9,1)-volc_base_time)/86400))
ax0.set_ylabel('Detections/Week',color='black',fontsize=38)
ax0.tick_params(axis='y',labelsize=32)
ax0.set_title('Redoubt Detections (2008-04-01 to 2009-09-01)',fontsize=50,fontweight='bold')
plt.yscale('log', nonpositive='clip')
ax0b = ax0.twinx()
ax0b.text(31, 4200, '(Network\nImpaired)', fontsize=35)
ax0b.plot(cum_volc_PARTY_days,cum_volc_PARTY_events,'-',color='black',linewidth=8,label='Cumulative EQcorrscan detections')
ax0b.plot(cum_volc_PEC_days,cum_volc_PEC_events,'-',color='darkred',linewidth=8,label='Cumulative Original Catalog Count',zorder=9)
ax0b.plot(cum_volc_PARTY_days[-1],cum_volc_PARTY_events[-1],'o',color='#505050')
ax0b.plot(cum_volc_PEC_days[-1],cum_volc_PEC_events[-1],'o',color='darkred')
# ax0b.legend(fontsize=35,loc='lower left')
ax0b.set_ylim([0,35000])
ax0b.set_ylabel('Cumulative Count',color='black',fontsize=38)
ax0b.tick_params(axis='y',labelsize=32)
# Full duration plot
ax1 = plt.subplot2grid((8, 3), (2, 0), colspan=3, rowspan=3)
base_time = UTCDateTime(2008,4,1)
end_time = UTCDateTime(2009,9,1)
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
zoom1_base_time = UTCDateTime(2009,3,14)
zoom1_end_time = UTCDateTime(2009,4,5)
zoom2_base_time = UTCDateTime(2009,5,2)
zoom2_end_time = UTCDateTime(2009,5,9)
tick_times = [UTCDateTime(2008, 4, 1),
              UTCDateTime(2008, 5, 1), UTCDateTime(2008, 6, 1),
              UTCDateTime(2008, 7, 1), UTCDateTime(2008, 8, 1),
              UTCDateTime(2008, 9, 1), UTCDateTime(2008, 10, 1),
              UTCDateTime(2008, 11, 1), UTCDateTime(2008, 12, 1),
              UTCDateTime(2009, 1, 1), UTCDateTime(2009, 2, 1),
              UTCDateTime(2009, 3, 1), UTCDateTime(2009, 4, 1),
              UTCDateTime(2009, 5, 1), UTCDateTime(2009, 6, 1),
              UTCDateTime(2009, 7, 1), UTCDateTime(2009, 8, 1),
              UTCDateTime(2009, 9, 1)]
ax1.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
ax1.scatter(relocated_days, relocated_FI_values, s=relocated_sizes, c='indianred', alpha=0.4)
ax1.scatter(original_days, original_FI_values, s=original_sizes, c='darkred', alpha=0.7)
ax1.scatter(8, -1.3, s=750, marker='s', c='#505050', alpha=0.5)
ax1.scatter(8, -1.1, s=750, marker='s', c='indianred', alpha=0.5)
ax1.scatter(8, -0.9, s=750, marker='s', c='darkred')
ax1.scatter(8, -1.3, s=750, marker='s', c='None', edgecolor='k')
ax1.scatter(8, -1.1, s=750, marker='s', c='None', edgecolor='k')
ax1.scatter(8, -0.9, s=750, marker='s', c='None', edgecolor='k')
ax1.text(10.7, -1.3, ': Enhanced Catalog w/o Magnitudes', fontsize=35, va='center')
ax1.text(10.7, -1.1, ': Enhanced Catalog with Magnitudes', fontsize=35, va='center')
ax1.text(10.7, -0.9, ': Original AVO Catalog', fontsize=35, va='center')
ax1.scatter(19, -0.6, s=50, marker='o', c='None', edgecolor='k')
ax1.scatter(30, -0.6, s=250, marker='o', c='None', edgecolor='k')
ax1.scatter(41, -0.6, s=450, marker='o', c='None', edgecolor='k')
ax1.scatter(52, -0.6, s=650, marker='o', c='None', edgecolor='k')
ax1.text(19, -0.4, 'M0', fontsize=35, va='center', ha='center')
ax1.text(30, -0.4, 'M1', fontsize=35, va='center', ha='center')
ax1.text(41, -0.4, 'M2', fontsize=35, va='center', ha='center')
ax1.text(52, -0.4, 'M3', fontsize=35, va='center', ha='center')
ax1.text(3,0.48,'b)',fontsize=58,weight='bold')
ax1.set_xlim([0,num_days])
ax1.set_ylim([-1.5,1.0])
ax1.set_xticks(np.array([(tt-base_time)/86400 for tt in tick_times], dtype=float))
ax1.set_xticklabels([tick_time.strftime('%b \'%y') for tick_time in tick_times],fontsize=35,rotation=30)
ax1.tick_params(axis='y',labelsize=32)
ax1.tick_params(which='major',length=10,width=2,direction='in')
ax1.set_ylabel('Frequency Index',fontsize=38)
ax1.set_title('Variation of earthquake FI with time (N=%d)' % len(catalog),fontsize=50,fontweight='bold')
explosion_days = [(explosion_time - base_time) / 86400 for explosion_time in explosion_times]
for explosion_day in explosion_days:
    ax1.axvline(x=explosion_day, color='k', alpha=0.6)
color_code_days = [(color_code_time - base_time) / 86400 for color_code_time in color_code_times]
color_code_days[0] = 0 # replace first color code with start of time
for k in range(len(color_code_days) - 1):
    ax1.axhspan(ymin=0.8, ymax=1.00, xmin=color_code_days[k] / num_days, xmax=color_code_days[k + 1] / num_days,
               color=color_codes[k], zorder=10)
    ax1.axhspan(ymin=-1.5, ymax=-1.45, xmin=(zoom1_base_time - base_time) / (86400 * num_days),
                xmax=(zoom1_end_time - base_time) / (86400 * num_days), color='blue', alpha=0.5)
    ax1.axhspan(ymin=-1.5, ymax=-1.45, xmin=(zoom2_base_time - base_time) / (86400 * num_days),
                xmax=(zoom2_end_time - base_time) / (86400 * num_days), color='blue', alpha=0.5)
# Zoom 1 plot
ax2 = plt.subplot2grid((8, 3), (5, 0), colspan=2, rowspan=3)
base_time = zoom1_base_time
end_time = zoom1_end_time
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(0,int(num_days+1),2)]
ax2.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
ax2.scatter(relocated_days,relocated_FI_values,s=relocated_sizes,c='indianred',alpha=0.4)
ax2.scatter(original_days,original_FI_values,s=original_sizes,c='darkred',alpha=0.7)
ax2.text(0.2,0.46,'c)',fontsize=58,weight='bold')
ax2.set_xlim([0,num_days])
ax2.set_ylim([-1.5,1.0])
ax2.set_xticks(np.array([(tt-base_time)/86400 for tt in tick_times], dtype=float))
ax2.set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times],fontsize=35,rotation=30,ha='right')
ax2.tick_params(axis='y',labelsize=32)
ax2.tick_params(which='major',length=10,width=2,direction='in')
ax2.set_ylabel('Frequency Index',fontsize=38)
# ax.set_xlabel('UTC Time',fontsize=15)
num_events = len(utctimes[(utctimes > base_time) & (utctimes < end_time)]) + len(relocated_utctimes[(relocated_utctimes > base_time) & (relocated_utctimes < end_time)])
ax2.set_title('(N=%d)' % num_events,fontsize=50,fontweight='bold')
explosion_days = [(explosion_time - base_time) / 86400 for explosion_time in explosion_times]
for explosion_day in explosion_days:
    ax2.axvline(x=explosion_day, color='k', alpha=0.6)
color_code_days = [(color_code_time - base_time) / 86400 for color_code_time in color_code_times]
for k in range(len(color_code_days) - 1):
    ax2.axhspan(ymin=0.8, ymax=1.00, xmin=color_code_days[k] / num_days, xmax=color_code_days[k + 1] / num_days,
               color=color_codes[k], zorder=10)
# Zoom 2 plot
ax3 = plt.subplot2grid((8, 3), (5, 2), colspan=1, rowspan=3)
base_time = zoom2_base_time
end_time = zoom2_end_time
days = (utctimes-base_time)/86400
relocated_days = (relocated_utctimes-base_time)/86400
original_days = (original_utctimes-base_time)/86400
num_days = (end_time-base_time)/86400
tick_times = [base_time+86400*i for i in range(int(num_days+1))]
ax3.scatter(days, FI_values, s=50, c='#505050', alpha=0.5)
ax3.scatter(relocated_days,relocated_FI_values,s=relocated_sizes,c='indianred',alpha=0.4)
ax3.scatter(original_days,original_FI_values,s=original_sizes,c='darkred',alpha=0.7)
ax3.text(0.1,0.46,'d)',fontsize=58,weight='bold')
ax3.set_xlim([0,num_days])
ax3.set_ylim([-1.5,1.0])
ax3.set_xticks(np.array([(tt-base_time)/86400 for tt in tick_times], dtype=float))
ax3.set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times],fontsize=35,rotation=30,ha='right')
ax3.tick_params(axis='y',labelsize=32)
ax3.tick_params(which='major',length=10,width=2,direction='in')
# ax3.set_ylabel('Frequency Index',fontsize=35)
# ax3.set_xlabel('UTC Time',fontsize=15)
num_events = len(utctimes[(utctimes > base_time) & (utctimes < end_time)]) + len(relocated_utctimes[(relocated_utctimes > base_time) & (relocated_utctimes < end_time)])
ax3.set_title('(N=%d)' % num_events,fontsize=50,fontweight='bold')
explosion_days = [(explosion_time - base_time) / 86400 for explosion_time in explosion_times]
for explosion_day in explosion_days:
    ax3.axvline(x=explosion_day, color='k', alpha=0.6)
color_code_days = [(color_code_time - base_time) / 86400 for color_code_time in color_code_times]
for k in range(len(color_code_days) - 1):
    ax3.axhspan(ymin=0.8, ymax=1.00, xmin=color_code_days[k] / num_days, xmax=color_code_days[k + 1] / num_days,
               color=color_codes[k], zorder=10)
plt.tight_layout()
plt.show()
fig.savefig('/Users/darrentpk/Desktop/figures/paper/new/fig4_redoubt_fi_time_series_v7.pdf',bbox_inches='tight')
fig.savefig('/Users/darrentpk/Desktop/figures/paper/new/fig4_redoubt_fi_time_series_v7.png',bbox_inches='tight')


####################################################################
