#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2009-01-01/2009-02-01/2009-03-01/2009-04-01/2009-05-01/2009-06-01/2009-07-01/2009-08-01/2009-09-01/2009-10-01/2009-11-01/2009-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "176.128/177.128/178.128/179.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_rad_2009.nc",
})
