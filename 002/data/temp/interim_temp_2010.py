#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2010-01-01/2010-02-01/2010-03-01/2010-04-01/2010-05-01/2010-06-01/2010-07-01/2010-08-01/2010-09-01/2010-10-01/2010-11-01/2010-12-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "pl",
    "levelist": "100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    "param": "130.128",
    "stream": "moda",
    "type": "an",
    "format": "netcdf",
    "target": "interim_temp_2010.nc",
})
