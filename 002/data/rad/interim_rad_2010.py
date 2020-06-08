#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2010-01-01/2010-02-01/2010-03-01/2010-04-01/2010-05-01/2010-06-01/2010-07-01/2010-08-01/2010-09-01/2010-10-01/2010-11-01/2010-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "176.128/177.128/178.128/179.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_rad_2010.nc",
})
