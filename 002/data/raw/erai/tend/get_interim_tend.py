#!/bin/python
import sys

year = list(sys.argv)[1] # argument after the script name is the year

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": ("{yr}-01-01/to/{yr}-12-31").format(yr=year),
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "62.162",
    "step": "0",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "stream": "oper",
    "type": "an",
    "format": "netcdf",
    "target": ("interim_tend_{yr}.nc").format(yr=year),
})
