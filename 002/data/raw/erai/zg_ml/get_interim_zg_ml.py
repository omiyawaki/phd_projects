#!/usr/bin/env python
import sys

year = list(sys.argv)[1] # second argument is the year

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": ("{yr}-01-01/{yr}-02-01/{yr}-03-01/{yr}-04-01/{yr}-05-01/{yr}-06-01/{yr}-07-01/{yr}-08-01/{yr}-09-01/{yr}-10-01/{yr}-11-01/{yr}-12-01").format(yr = year),
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "ml",
    "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60",
    "param": "129.128",
    "stream": "moda",
    "type": "an",
    "format": "netcdf",
    "target": ("interim_zg_ml_{yr}.nc").format(yr = year),
})
