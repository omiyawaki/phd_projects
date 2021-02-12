#!/usr/bin/env python
import sys

year = list(sys.argv)[1] # second argument is the year

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": ("{yr}-01-01/to/{yr}-12-31").format(yr = year),
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    "param": "133.128",
    "stream": "oper",
    "step": "0",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "format": "netcdf",
    "target": ("interim_hus_{yr}.nc").format(yr = year),
})
