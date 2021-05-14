from eqcorrscan.utils.catalog_to_dd import _generate_event_id_mapper, _compute_dt_correlations, _make_sparse_event, \
    _compute_dt, compute_differential_times
from eqcorrscan.utils.clustering import dist_mat_km
from obspy import Catalog
import pickle
from toolbox import reader

catalog_inpath = '/home/ptan/attempt_eqcorrscan/output/'
#catalog = reader(catalog_inpath+'repicked_catalog.xml')
catalog = reader(catalog_inpath+'core_catalog_picked.xml')
with open('/home/ptan/attempt_eqcorrscan/avo_data/redpy/clusters_NA_AVO.txt', 'rb') as cluster_pickle:  # Unpickling
    clusters_NA = pickle.load(cluster_pickle)
new_catalog = Catalog()
for event in catalog:
    if int(event.origins[0].comments[0].text.split(' ')[1]) not in clusters_NA:
        new_catalog += event
catalog = new_catalog

correlation = True
stream_dict = None
event_id_mapper = None
max_sep = 8
min_link = 3
min_cc = None
extract_len = None
pre_pick = 2
shift_len = 3
interpolate = False
max_workers = None

include_master = False
correlation_kwargs = dict(
    min_cc=min_cc, stream_dict=stream_dict, extract_len=extract_len,
    pre_pick=pre_pick, shift_len=shift_len, interpolate=interpolate,
    max_workers=max_workers)
if correlation:
    for arg, name in correlation_kwargs.items():
        assert arg is not None, "{0} is required for correlation".format(name)
# Ensure all events have locations and picks.
event_id_mapper = _generate_event_id_mapper(
    catalog=catalog, event_id_mapper=event_id_mapper)
distances = dist_mat_km(catalog)
distance_filter = distances <= max_sep
if not include_master:
    np.fill_diagonal(distance_filter, 0)
    # Do not match events to themselves - this is the default,
    # only included for testing

additional_args = dict(min_link=min_link, event_id_mapper=event_id_mapper)
if correlation:
    differential_times = {}
    additional_args.update(correlation_kwargs)
    n = len(catalog)
    for i, master in enumerate(catalog):
        master_id = master.resource_id.id
        sub_catalog = [ev for j, ev in enumerate(catalog)
                       if distance_filter[i][j]]
        if master_id not in additional_args["stream_dict"].keys():
            Logger.warning(
                f"{master_id} not in waveforms, skipping")
            continue
        differential_times.update({
            master_id: _compute_dt_correlations(
                sub_catalog, master, **additional_args)})
        Logger.info(
            f"Completed correlations for core event {i} of {n}")
else:
    # Reformat catalog to sparse catalog
    sparse_catalog = [_make_sparse_event(ev) for ev in catalog]

    sub_catalogs = ([ev for i, ev in enumerate(sparse_catalog)
                     if master_filter[i]]
                    for master_filter in distance_filter)
    differential_times = {
        master.resource_id: _compute_dt(
            sub_catalog, master, **additional_args)
        for master, sub_catalog in zip(sparse_catalog, sub_catalogs)}

# Remove Nones
for key, value in differential_times.items():
    differential_times.update({key: [v for v in value if v is not None]})
return differential_times, event_id_mapper