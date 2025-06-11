## Download GFS:

The GFS Forecast data is stored in an Amazon S3 bucket in us-east-1 AWS region at (https://registry.opendata.aws/noaa-gfs-bdp-pds/)[https://registry.opendata.aws/noaa-gfs-bdp-pds/].

- `noaa-gfs-bpd-pds`

### pgrb2

The filename structure for `pgrb2` is as follows:

`noaa-gfs-bdp-pds/gfs.YYYYMMDD/XX/atmos/gfs.tXXz.pgrb2.0p25.fVVV`

- Where `YYYY` = year ,`MM` = month, `DD` = day of the initialization date
- Where `XX` will be (00, 06, 12, 18) corresponding to the four forecasts initialized each day
- Where `VVV` is the forecast lead in hours (000,006, .. , 384)

---

Variables available for `pgrb2` :

https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.0p25.anl.shtml