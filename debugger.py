

writer(scan_data_output_dir + 'master_catalogM.xml',master_catalog)


keys = []
mags = []
for event in master_catalog:
    keys.append(event.comments[0:4])
    mags.append(event.magnitudes)

for event in relocatable_catalog:
    key = event.comments[0:4]
    try:
        master_index = keys.index(key)
        event.magnitudes = master_catalog[master_index].magnitudes
    except ValueError:
        continue

