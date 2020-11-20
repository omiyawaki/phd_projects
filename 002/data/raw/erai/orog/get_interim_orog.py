#!/usr/bin/env python
import sys

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1989-01-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "129.128",
    "step": "0",
    "time": "12:00:00",
    "stream": "oper",
    "type": "an",
    "format": "netcdf",
    "target": "interim_orog.nc",
})
