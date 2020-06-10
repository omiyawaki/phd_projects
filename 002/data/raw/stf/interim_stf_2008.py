#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2008-01-01/2008-02-01/2008-03-01/2008-04-01/2008-05-01/2008-06-01/2008-07-01/2008-08-01/2008-09-01/2008-10-01/2008-11-01/2008-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_stf_2008.nc",
})
