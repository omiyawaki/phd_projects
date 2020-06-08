#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2003-01-01/2003-02-01/2003-03-01/2003-04-01/2003-05-01/2003-06-01/2003-07-01/2003-08-01/2003-09-01/2003-10-01/2003-11-01/2003-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_tf_2003.nc",
})
