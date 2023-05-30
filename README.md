enhance_catalog
============

This repository stores the codes that make up the workflow described in Tan et al. (2023): 

*Tan, D., Fee, D., Hotovec-Ellis, A. J., Pesicek, J. D., Haney, M. M., Power, J. A. & Girona, T. (2021). Volcanic earthquake catalog enhancement using integrated detection, matched-filtering, and relocation tools. [doi: 10.3389/feart.2023.1158442](https://doi.org/10.3389/feart.2023.1158442)*

The codes can be used to enhance AVO earthquake catalogs using a streamlined integration of open-source catalog-enhancing tools: REDPy, EQcorrscan, HypoDD, and GrowClust. The combination of these tools offers the capability of adding seismic event detections and relocating events in a single workflow. The workflow relies on a combination of standard triggering and cross-correlation clustering (REDPy) to consolidate representative templates used in matched-filtering (EQcorrscan). The templates and their detections are then relocated using the differential time methods provided by HypoDD and/or GrowClust. Additional utilities include the incorporation of campaign and backfilled data at appropriate junctures, relative magnitude calculations, and frequency index calculations for valid events.

The general sequence to be followed is:
1. Determine the timespan of interest and download seismic data
2. Run REDPy to determine families of repeaters within the timespan
3. Read in the analyst-derived catalog and compare it with REDPy's output to consolidate representative templates
4. Create the tribe of representative templates for EQcorrscan
5. Run the matched-filter scan using EQcorrscan to obtain the temporally enhanced event list
6. Calculate frequency indices and relative magnitudes for all valid events
7. Cross-correlate relocation candidates with one another to obtain cross-correlation differential times
8. Relocate valid events to obtain the relocated, enhanced catalog
9. Plot resulting data products in time and in space



**Key references and other suggested citations**

*Hotovec-Ellis, A. J., & Jeffries, C. (2016). Near real‐time detection, clustering, and analysis of repeating earthquakes: Application to Mount St. Helens and Redoubt volcanoes.*

*Chamberlain, C. J., Hopp, C. J., Boese, C. M., Warren‐Smith, E., Chambers, D., Chu, S. X., ... & Townend, J. (2018). EQcorrscan: Repeating and near‐repeating earthquake detection and analysis in Python.*

*Trugman, D. T., & Shearer, P. M. (2017). GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.*

*Klein, F. W. (2002). User's guide to HYPOINVERSE-2000, a Fortran program to solve for earthquake locations and magnitudes (No. 2002-171). US Geological Survey.*

*Waldhauser, F. (2001). hypoDD--A program to compute double-difference hypocenter locations.*

Quickstart
----------

1. Clone and enter directory

```
git clone https://github.com/darren-uaf/enhance_catalog
cd enhance_catalog
```

2. Create conda environment

```
conda env create
conda activate enhance_catalog
```

3. Run example (3 day period from the Redoubt 2009 eruption)

```
python example.py
```

Dependencies
------------

Other repositories:
* [REDPy](https://github.com/ahotovec/REDPy)
* [EQcorrscan](https://github.com/eqcorrscan/EQcorrscan)
* [GrowClust](https://github.com/dttrugman/GrowClust)
* [waveform_collection](https://github.com/uafgeotools/waveform_collection)


<!--stackedit_data:
eyJwcm9wZXJ0aWVzIjoiZXh0ZW5zaW9uczpcbiAgcHJlc2V0Oi
BnZm1cbiAgbWFya2Rvd246XG4gICAgYnJlYWtzOiBmYWxzZVxu
IiwiaGlzdG9yeSI6WzYxMTk4MTkwMCwxOTg3MzQ1MzEwLDQzMD
M3MzM1OSw0MzAzNzMzNTldfQ==
-->
