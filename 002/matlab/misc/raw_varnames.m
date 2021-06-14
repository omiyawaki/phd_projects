par.rea.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.rea.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.rea.vars.hydro = {'pr', 'evspsbl'}; % radiation variables to read
par.rea.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.rea.vars.vert = {'ta'}; % 3d variables to read (removed va)
par.rea.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'}; % surface variables to read
par.rea.vars.tend = {'tendmon', 'tend'}; 
par.rea.vars.tend_txt = {'tend', 'tend'}; % 3d variables to read (tend = total mse tend, int = internal energy tendency, lat = latent energy tendency, pot = potential energy tendency)
par.era.vars.rad = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
par.era.vars.radcs = {'ssrc', 'strc', 'tsrc', 'ttrc'}; % radiation variables to read
par.era.vars.hydro = {'cp', 'lsp', 'e'}; % radiation variables to read
par.era.vars.div = {'p85.162', 'p84.162', 'p83.162'}; % radiation variables to read
par.era.vars.div_txt = {'divg', 'divq', 'divt'}; % radiation variables to read
par.era.vars.stf = {'sshf', 'slhf'}; % surface turbulent flux variables to read
par.era.vars.vert = {'t'}; % 3d variables to read (t = temp)
par.era.vars.srfc = {'sp', 'skt', 't2m', 'd2m', 'zs',}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
par.era.vars.tend = {'tend', 'tend'}; 
par.era.vars.tend_txt = {'tendalt', 'tendmon'}; % 3d variables to read (tend = total mse tend, int = internal energy tendency, lat = latent energy tendency, pot = potential energy tendency)
par.era5c.vars.rad = {'ssr', 'str', 'tsr', 'ttr'}; % radiation variables to read
par.era5c.vars.radcs = {'ssrc', 'strc', 'tsrc', 'ttrc'}; % radiation variables to read
par.era5c.vars.hydro = {'cp', 'lsp', 'e'}; % radiation variables to read
par.era5c.vars.div = {'p85.162', 'p84.162', 'p83.162'}; % radiation variables to read
par.era5c.vars.div_txt = {'divg', 'divq', 'divt'}; % radiation variables to read
par.era5c.vars.stf = {'sshf', 'slhf'}; % surface turbulent flux variables to read
par.era5c.vars.vert = {'t'}; % 3d variables to read (t = temp)
par.era5c.vars.srfc = {'sp', 'skt', 't2m', 'd2m', 'zs',}; % surface variables to read (sp = surface pressure, t2m = 2m temp, d2m = 2m dew point temp)
% par.era.vars.tend = {'p62.162', 'p60.162', 'qL', 'gz'}; 
par.era5c.vars.tend = {'tend', 'tend'}; 
par.era5c.vars.tend_txt = {'tendalt', 'tendmon'}; % 3d variables to read (tend = total mse tend, int = internal energy tendency, lat = latent energy tendency, pot = potential energy tendency)
% par.era.vars.tend_txt = {'tend', 'int', 'lat', 'pot'}; % 3d variables to read (tend = total mse tend, int = internal energy tendency, lat = latent energy tendency, pot = potential energy tendency)
par.hahn.vars.rad = {'FLNT', 'FLNS', 'FSNT', 'FSNS'}; % radiation variables to read
par.hahn.vars.hydro = {'PRECC', 'PRECL', 'PRECSC', 'PRECSL'}; % hydrology variables
par.hahn.vars.stf = {'SHFLX', 'LHFLX'};
par.hahn.vars.vert = {'T'};
par.hahn.vars.srfc = {'PS', 'TREFHT', 'TS', 'huss', 'zs'};
par.merra2c.vars.rad = {'SWTNT', 'SWGNT', 'LWTUP', 'LWGNT'}; % radiation variables to read
par.merra2c.vars.radcs = {'SWTNTCLR', 'SWGNTCLR', 'LWTUPCLR', 'LWGNTCLR'}; % radiation variables to read
par.merra2c.vars.hydro = {'PRECTOT', 'PRECCON', 'EVAP'}; % hydrology variables
par.merra2c.vars.stf = {'HFLUX', 'EFLUX'};
par.merra2c.vars.vert = {'T'};
par.merra2c.vars.srfc = {'PS', 'TS', 'T2M', 'QV2M', 'zs'};
par.merra2c.vars.tend = {'tend'}; 
par.merra2c.vars.tend_txt = {'tendmon'}; 
% par.jra55.vars.rad = {'dswrf', 'uswrf', 'dlwrf', 'ulwrf'}; % radiation variables to read
par.jra55.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.jra55.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.jra55.vars.hydro = {'pr', 'prc', 'evspsbl'}; % hydrology variables
% par.jra55.vars.stf = {'lhtfl', 'shtfl'};
par.jra55.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.jra55.vars.vert = {'ta'};
par.jra55.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'};
par.jra55.vars.tend = {'tend'}; 
par.jra55.vars.tend_txt = {'tendmon'}; 
par.gcm.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.gcm.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.gcm.vars.hydro = {'pr', 'evspsbl'}; % radiation variables to read
par.gcm.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.gcm.vars.vert = {'ta'}; % 3d variables to read (removed va)
par.gcm.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'}; % surface variables to read
par.gcm.vars.tend = {'tendmon', 'tend'}; 
par.gcm.vars.tend_txt = {'tend', 'tend'}; % 3d variables to read (tend = total mse tend, int = internal energy tendency, lat = latent energy tendency, pot = potential energy tendency)
par.echam.vars.rad = {'srad0d', 'srads', 'trads', 'srad0', 'trad0', 'tradsu', 'sradsu'}; % radiation variables to read
par.echam.vars.radcs = {'srafs', 'trafs', 'sraf0', 'traf0'}; % radiation variables to read
par.echam.vars.hydro = {'aprc', 'aprl', 'evap'}; % radiation variables to read
par.echam.vars.stf = {'ahfl', 'ahfs'}; % surface turbulent flux variables to read
par.echam.vars.vert = {'t', 'v'}; % 3d variables to read
par.echam.vars.srfc = {'aps', 'tsurf', 'temp2', 'dew2'}; % surface variables to read
par.echam.vars.tend = {'tendmon'}; 
par.echam.vars.tend_txt = {'tend'}; 
par.ceres.vars.rad = {'sfc_net_sw_all_mon', 'sfc_net_lw_all_mon', 'toa_sw_all_mon', 'solar_mon', 'toa_lw_all_mon'}; % radiation variables to read
par.ceres.vars.rad_txt = {'ssr', 'str', 'tsur', 'tsdr', 'ttr'}; % radiation variables to read