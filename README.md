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
search_surveys "188.733 -12.5822" -r 10
Searching 10.0deg around 12:34:55.92 -12:34:55.92 = 188.733d,-12.5822d
Found 5 matches:
   PSR         RA        Dec        P       DM    survey       retrieval date       Distance
              deg        deg        ms   pc / cm3                                     deg   
---------- ---------- ---------- ------- -------- ------ -------------------------- --------
J1231-1411 187.797083 -14.195444    3.68     8.09   ATNF 2020-12-11 14:22:10.536850   1.8524
  B1254-10 194.269583 -10.451528  617.31    29.63   ATNF 2020-12-11 14:22:10.536850   5.8281
  J1244-18 191.070833 -18.286667 3425.10    17.00 SUPERB 2020-12-11 14:21:41.961846   6.1330
  B1309-12 197.969167 -12.467111  447.52    36.21   ATNF 2020-12-11 14:22:10.536850   9.0166
  J1327-07 196.750000  -7.000000    2.68    27.90 GBT350 2020-12-11 14:21:39.893771   9.6705
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
* [GalacticMSPs](http://astro.phys.wvu.edu/GalacticMSPs/GalacticMSPs.txt)
* [TRAPUM](http://www.trapum.org/discoveries.html)
* [GHRSS](http://www.ncra.tifr.res.in/~bhaswati/GHRSS.html)
