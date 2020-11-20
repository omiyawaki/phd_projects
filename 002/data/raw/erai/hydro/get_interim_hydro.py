#!/usr/bin/env python
import sys

year = list(sys.argv)[1] # argument after the script name is the year

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": ("{yr}-01-01/{yr}-02-01/{yr}-03-01/{yr}-04-01/{yr}-05-01/{yr}-06-01/{yr}-07-01/{yr}-08-01/{yr}-09-01/{yr}-10-01/{yr}-11-01/{yr}-12-01").format(yr=year),
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "142.128/143.128/182.128",
    "step": "12",
    "time": "00:00:00/12:00:00",
    "stream": "mdfa",
    "type": "fc",
    "format": "netcdf",
    "target": ("interim_hydro_{yr}.nc").format(yr=year),
})
