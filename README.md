# pulsarsurveyscraper
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6390905.svg)](https://doi.org/10.5281/zenodo.6390905)

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

* Also do a CHIME sideline arc search (this picks up the Crab):
```
(scraper) kaplan@plock[~/pythonpackages/pulsarsurveyscraper] (main) % search_surveys -r 5 --sidelobe-arc CHIME 61216 -v 65.404228 24.041539
INFO: Read survey 'ATNF' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/ATNF.hdf5: found 4393 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GalacticMSPs' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GalacticMSPs.hdf5: found 680 pulsars [pulsarsurveyscraper]
INFO: Read survey 'AO327' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/AO327.hdf5: found 105 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GBNCC' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GBNCC.hdf5: found 194 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GBT820' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GBT820.hdf5: found 7 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GBT350' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GBT350.hdf5: found 34 pulsars [pulsarsurveyscraper]
INFO: Read survey 'PALFA' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/PALFA.hdf5: found 208 pulsars [pulsarsurveyscraper]
INFO: Read survey 'DMB' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/DMB.hdf5: found 29 pulsars [pulsarsurveyscraper]
INFO: Read survey 'RRATalog' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/RRATalog.hdf5: found 337 pulsars [pulsarsurveyscraper]
INFO: Cannot find cached survey 'CHIME' at /Users/kaplan/pythonpackages/pulsarsurveyscraper/CHIME.hdf5: loading... [pulsarsurveyscraper]
ERROR: Failed to reach server 'http://catalog.chime-frb.ca/galactic': [Errno 8] nodename nor servname provided, or not known [pulsarsurveyscraper]
INFO: Read survey 'CRAFTS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/CRAFTS.hdf5: found 222 pulsars [pulsarsurveyscraper]
INFO: Read survey 'FAST-GPPS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/FAST-GPPS.hdf5: found 850 pulsars [pulsarsurveyscraper]
INFO: Read survey 'MWA' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/MWA.hdf5: found 28 pulsars [pulsarsurveyscraper]
INFO: Read survey 'TRAPUM' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/TRAPUM.hdf5: found 367 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GHRSS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GHRSS.hdf5: found 31 pulsars [pulsarsurveyscraper]
INFO: Read survey 'VAST' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/VAST.hdf5: found 8 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GCGPS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GCGPS.hdf5: found 7 pulsars [pulsarsurveyscraper]
INFO: Read survey 'PKS70' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/PKS70.hdf5: found 1 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GM-SCOPE' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GM-SCOPE.hdf5: found 13 pulsars [pulsarsurveyscraper]
INFO: Read survey 'CHAMPSS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/CHAMPSS.hdf5: found 60 pulsars [pulsarsurveyscraper]
INFO: Read survey 'SSP' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/SSP.hdf5: found 8 pulsars [pulsarsurveyscraper]
INFO: Read survey 'MGBS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/MGBS.hdf5: found 1 pulsars [pulsarsurveyscraper]
INFO: Read survey 'Puschino' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/Puschino.hdf5: found 89 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GBNCCR' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GBNCCR.hdf5: found 51 pulsars [pulsarsurveyscraper]
INFO: Read survey 'LMCR' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/LMCR.hdf5: found 10 pulsars [pulsarsurveyscraper]
INFO: Read survey 'FAST-GLEANERS' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/FAST-GLEANERS.hdf5: found 30 pulsars [pulsarsurveyscraper]
INFO: Read survey 'GRASP' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/GRASP.hdf5: found 5 pulsars [pulsarsurveyscraper]
INFO: Read survey 'NUMBAT' from /Users/kaplan/pythonpackages/pulsarsurveyscraper/NUMBAT.hdf5: found 3 pulsars [pulsarsurveyscraper]
Searching 5.0deg around RA,Dec = 04:21:37.01472 +24:02:29.5404 = 65.4042d,+24.0415d along the sidelobe arc
INFO: Using transit telescope mode: 49.32070922 deg -119.62367743000004 deg 61216.0 [pulsarsurveyscraper]
Downloading ftp://anonymous:mail%40astropy.org@gdc.cddis.eosdis.nasa.gov/pub/products/iers/finals2000A.all
|===========================================================================================================================================================================================| 3.7M/3.7M (100.00%)         0s
Found 91 matches:
    PSR          RA        Dec         P        DM       survey          retrieval date       Distance Sidelobe Arc Distance
                deg        deg         ms    pc / cm3                                           deg             deg         
------------ ---------- ---------- --------- -------- ------------ -------------------------- -------- ---------------------
  J0846-1049 131.500000 -10.816667    330.10    13.70       CRAFTS 2026-09-04 17:45:38.400849  73.3199                0.0041
    J0457+23  74.275000  23.566667    504.90    59.00         ATNF 2026-09-08 13:31:09.053337   8.1287                0.0317
  J0350+2341  57.512500  23.683333   2421.20    61.00     Puschino 2026-09-04 17:45:54.334642   7.2251                0.0766
  J0349+2340  57.488750  23.681389   2420.77    62.96         ATNF 2026-09-08 13:31:09.053337   7.2469                0.0773
    J0625+17  96.333333  17.500000   2518.41    58.00     RRATalog 2026-09-04 17:45:35.924377  29.5833                0.0844
    B0525+21  82.217500  22.001111   3745.54    50.87         ATNF 2026-09-08 13:31:09.053337  15.5985                0.1601
    B0531+21  83.632917  22.014472     33.39    56.77         ATNF 2026-09-08 13:31:09.053337  16.8860                0.1805
...
```

## Requirements:
* `astropy`, `scipy`, `h5py`, `pyyaml`, `beautifulsoup4`, `requests`, `fpdf2`, `pdf2image`
* `poppler` (via `conda-forge`)

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
