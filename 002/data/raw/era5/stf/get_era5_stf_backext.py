#!/usr/bin/env python

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        'variable': [
            'surface_latent_heat_flux', 'surface_sensible_heat_flux',
        ],
        'year': [
            '1950', '1951', '1952',
            '1953', '1954', '1955',
            '1956', '1957', '1958',
            '1959', '1960', '1961',
            '1962', '1963', '1964',
            '1965', '1966', '1967',
            '1968', '1969', '1970',
            '1971', '1972', '1973',
            '1974', '1975', '1976',
            '1977', '1978',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
    },
    'era5_stf_1950_1978.nc')
