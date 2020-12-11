# pulsarsurveyscraper

## Modules:
* `pulsarsurveyscraper`

## Scripts:
* `cache_survey`: download a survey table and save to a file (HDF5).  Example: 
```
cache_survey -v -s GBNCC
```
* `search_surveys`: do a cone search among all surveys, loading from cached versions if possible.  Example:
```
search_surveys "188.733 -12.5822"
Searching 5.0deg around 12:34:55.92 -12:34:55.92 = 188.733d,-12.5822d
Found 1 matches:
  PSR        RA       Dec      P       DM    survey            date            Distance
            deg       deg      ms   pc / cm3                                     deg   
-------- --------- --------- ------ -------- ------ -------------------------- --------
J0110+11 15.000000 11.000000 432.09    14.90  AO327 2020-12-11 14:21:36.805394   1.8524
```

## Requirements:
* `astropy`, `json`, `h5py`, `pyyaml`, `beautifulsoup4`, `requests`

## Surveys:
* [AO327](http://www.naic.edu/~deneva/drift-search/index.html)
* [GBNCC](http://astro.phys.wvu.edu/GBNCC/)
* [GBT350](http://astro.phys.wvu.edu/GBTdrift350/)
* [PALFA](http://www2.naic.edu/~palfa/newpulsars/index.html)
* [DMB](http://astro.phys.wvu.edu/dmb)
* [SUPERB](https://sites.google.com/site/publicsuperb/discoveries)
* [HTRU-S Low-latitude](https://sites.google.com/site/htrusouthdeep/home/discoveries)
* [LOTAAS](http://old.astron.nl/lotaas/index.php?sort=1)
* [RRATalog](http://astro.phys.wvu.edu/rratalog)
* [CHIME](http://catalog.chime-frb.ca/galactic)
* [FAST](http://crafts.bao.ac.cn/pulsar/)
* [FAST-GPPS](http://zmtt.bao.ac.cn/GPPS)
* [MWA](https://wiki.mwatelescope.org/display/MP/SMART+survey+candidates)
* [ATNF](http://www.atnf.csiro.au/research/pulsar/psrcat/)
