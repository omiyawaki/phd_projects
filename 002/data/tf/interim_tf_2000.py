#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2000-01-01/2000-02-01/2000-03-01/2000-04-01/2000-05-01/2000-06-01/2000-07-01/2000-08-01/2000-09-01/2000-10-01/2000-11-01/2000-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_tf_2000.nc",
})
