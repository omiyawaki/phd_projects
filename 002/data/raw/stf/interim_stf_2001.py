#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2001-01-01/2001-02-01/2001-03-01/2001-04-01/2001-05-01/2001-06-01/2001-07-01/2001-08-01/2001-09-01/2001-10-01/2001-11-01/2001-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_stf_2001.nc",
})
