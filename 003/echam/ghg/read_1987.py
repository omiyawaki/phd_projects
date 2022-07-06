import xarray as xr

ds = xr.open_dataset('./ghg_rcp85_1765-2500_c100203.nc')
print(ds.CO2.sel(time=slice('1987-07-01', '1987-07-01')))
# print(ds.CH4.sel(time=slice('1987-07-01', '1987-07-01')))
# print(ds.N2O.sel(time=slice('1987-07-01', '1987-07-01')))
# print(ds.f11.sel(time=slice('1987-07-01', '1987-07-01')))
# print(ds.f12.sel(time=slice('1987-07-01', '1987-07-01')))

# print(ds.CO2.sel(time=slice('1850-07-01', '1850-07-01')))
