enhance_catalog
============

This repository stores the codes that can be used to enhance AVO earthquake catalogs using an amalgamation of open-source tools, which are REDPy (Hotovec-Ellis & Jeffries, 2016), EQcorrscan (Chamberlain et al., 2018) and GrowClust (Trugman & Shearer, 2017). The three tools are combined in a single workflow, in order to make the redetection and relocation process as seamless as possible. The workflow starts by tying the REDPy output catalog with the pre-existing AVO catalog, before consolidating REDPy cluster cores and unmatched AVO events as templates for EQcorrscan's matched-filtering. After that, each matched-filter detection is re-correlated with their parent template, in order to obtain CC coefficients and lag times for GrowClust's hiearchical clustering and relocation algorithm.

The general sequence to be followed is:
1. convert_redpy.py
2. create_tribe.py
3. scan_data.py
4. sensitivity_test.py (optional)
5. relocate_catalog.py
6. plot_hypo.py

**General References and Suggested Citations**

*Klein, F. W. (2002). User's guide to HYPOINVERSE-2000, a Fortran program to solve for earthquake locations and magnitudes (No. 2002-171). US Geological Survey.*

*Waldhauser, F. (2001). hypoDD--A program to compute double-difference hypocenter locations.*

*Hotovec-Ellis, A. J., & Jeffries, C. (2016). Near real‐time detection, clustering, and analysis of repeating earthquakes: Application to Mount St. Helens and Redoubt volcanoes.*

*Chamberlain, C. J., Hopp, C. J., Boese, C. M., Warren‐Smith, E., Chambers, D., Chu, S. X., ... & Townend, J. (2018). EQcorrscan: Repeating and near‐repeating earthquake detection and analysis in Python.*

*Trugman, D. T., & Shearer, P. M. (2017). GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.*


Quickstart
----------

1. Clone and enter directory

```
git clone https://github.com/darren-uaf/enhance_catalog
cd enhance_catalog
```

2. Create condo environment

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
* [phase_processing](https://github.com/darren-uaf/phase_processing)
* [waveform_collection](https://github.com/uafgeotools/waveform_collection)


Examples
--------

**coming soon**


<!--stackedit_data:
eyJwcm9wZXJ0aWVzIjoiZXh0ZW5zaW9uczpcbiAgcHJlc2V0Oi
BnZm1cbiAgbWFya2Rvd246XG4gICAgYnJlYWtzOiBmYWxzZVxu
IiwiaGlzdG9yeSI6WzYxMTk4MTkwMCwxOTg3MzQ1MzEwLDQzMD
M3MzM1OSw0MzAzNzMzNTldfQ==
-->
