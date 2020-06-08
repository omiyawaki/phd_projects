#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2011-01-01/2011-02-01/2011-03-01/2011-04-01/2011-05-01/2011-06-01/2011-07-01/2011-08-01/2011-09-01/2011-10-01/2011-11-01/2011-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": "interim_tf_2011.nc",
})
