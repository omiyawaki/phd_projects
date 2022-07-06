import xarray as xr

ds = xr.open_dataset('./greenhouse_rcp85.shifttime.nc')
print(ds.CO2.sel(time=[38,39]))
