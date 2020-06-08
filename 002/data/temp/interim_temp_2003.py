#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2003-01-01/2003-02-01/2003-03-01/2003-04-01/2003-05-01/2003-06-01/2003-07-01/2003-08-01/2003-09-01/2003-10-01/2003-11-01/2003-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    "param": "130.128",
    "stream": "moda",
    "type": "an",
    "format": "netcdf",
    "target": "interim_temp_2003.nc",
})
