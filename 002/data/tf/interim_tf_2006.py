#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2006-01-01/2006-02-01/2006-03-01/2006-04-01/2006-05-01/2006-06-01/2006-07-01/2006-08-01/2006-09-01/2006-10-01/2006-11-01/2006-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_tf_2006.nc",
})
