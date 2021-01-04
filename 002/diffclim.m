% makes files for the difference between two climates (clim2 - clim1)

addpath(genpath('/project2/tas1/miyawaki/matlab'));

par.clim1 = 'piControl';
par.clim2 = 'historical';

% make sure that the models chosen here is available for both climates
par.gcm.clim = 'hist-pi';
par.lat_interp = 'native'; % native: native model grid, don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
% par.gcm_models = {'ACCESS1-0'};
par.gcm_models = {'ACCESS1-0', 'ACCESS1-3',... % GCM model names
                  'bcc-csm1-1', 'bcc-csm1-1-m',...
                  'BNU-ESM', 'CanESM2', 'CCSM4',...
                  'CNRM-CM5', 'CNRM-CM5-2',...
                  'CSIRO-Mk3-6-0', 'FGOALS-g2', 'FGOALS-s2',...
                  'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',...
                  'GISS-E2-H', 'GISS-E2-R',...
                  'HadGEM2-ES', 'inmcm4',...
                  'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',...
                  'MIROC5', 'MIROC-ESM',...
                  'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P',...
                  'MRI-CGCM3', 'NorESM1-M'};

par.gcm.vars.rad = {'rsus', 'rsds', 'rlus', 'rlds', 'rsdt', 'rsut', 'rlut'}; % radiation variables to read
par.gcm.vars.radcs = {'rsuscs', 'rsdscs', 'rldscs', 'rsutcs', 'rlutcs'}; % radiation variables to read
par.gcm.vars.hydro = {'pr', 'evspsbl'}; % radiation variables to read
par.gcm.vars.stf = {'hfss', 'hfls'}; % surface turbulent flux variables to read
par.gcm.vars.vert = {'ta'}; % 3d variables to read (removed va)
par.gcm.vars.srfc = {'ps', 'ts', 'tas', 'hurs', 'zs'}; % surface variables to read
par.gcm.fw = {'mse', 'dse'};

%% call functions
for k=1:length(par.gcm_models); par.model=par.gcm_models{k};
    type='gcm';
    disp(par.model)
    run_func(type, par);
end

function run_func(type, par)
    % read_grid(type, par) % grid, i.e. lon, lat, plev
    % read_lfrac(type, par) % land fraction
    % diff_struct(type, par) % take difference of data that are in the form of depth 1 structs
    % diff_array(type, par) % take difference of data that are in the form of arrays
    % diff_flux_z(type, par)
    % diff_flux_t(type, par)
    % diff_flux_zt(type, par)
    diff_vh(type, par)
    diff_vh_mon(type, par)
end

function read_grid(type, par)
    filename = 'grid.mat';
    % read data net SW and LW radiation data downloaded from Era5
    % first read lon and lat vectors since this is different from the Donohoe grid
    if any(strcmp(type, {'era5', 'erai'}))
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'longitude'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'latitude'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 10^2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'level')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif any(strcmp(type, {'era5c'}))
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lat'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 10^2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'level')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'merra2')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_rad_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lat'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 1e2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_temp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lev')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'jra55')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_dswrf_%s.ymonmean.nc', type, type, par.(type).yr_span), 'g0_lon_2'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/rad/%s_dswrf_%s.ymonmean.nc', type, type, par.(type).yr_span), 'g0_lat_1'));
        grid.dim3 = grid.dim2;
        grid.dim3.plev = 1e2*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/temp/%s_tmp_%s.ymonmean.nc', type, type, par.(type).yr_span), 'lv_ISBL1')); % multiply by 100 to convert hPa to Pa
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid')
    elseif strcmp(type, 'gcm')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'tas', par.model, par.clim1));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.nc', par.model, 'ta', par.model, par.clim1));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'plev'));
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam')
        if contains(par.echam.clim, 'rp000')
            file.dim2=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/ATM_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        else
            file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/BOT*_%s_*.ymonmean.nc', par.echam.clim));
            file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam/ATM*_%s_*.ymonmean.nc', par.echam.clim));
        end
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        grid.dim3.z = par.z;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_ml')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.a=double(ncread(fullpath.dim3, 'hyam'));
        grid.dim3.b=double(ncread(fullpath.dim3, 'hybm'));
        grid.dim3.z = par.z;
        grid.dim3.plev = par.pa;
        % grid.dim3.si = 1e-5*([grid.dim3.a+grid.dim3.b*1e5; 1e5]);
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'echam_pl')
        file.dim2=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_*.ymonmean.nc'));
        file.dim3=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
        fullpath.dim2=sprintf('%s/%s', file.dim2.folder, file.dim2.name);
        fullpath.dim3=sprintf('%s/%s', file.dim3.folder, file.dim3.name);
        grid.dim2.lon=double(ncread(fullpath.dim2, 'lon'));
        grid.dim3.lon=double(ncread(fullpath.dim3, 'lon'));
        grid.dim2.lat=double(ncread(fullpath.dim2, 'lat'));
        grid.dim3.lat=double(ncread(fullpath.dim3, 'lat'));
        grid.dim3.plev=double(ncread(fullpath.dim3, 'lev'));
        grid.dim3.z = par.z;
        % grid.dim3.si = 1e-5*grid.dim3.plev;
        grid.dim3.si = 1e-5*par.pa;
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_pl');
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='grid.mat';
        save(sprintf('%s/%s', newdir, filename), 'grid');
    elseif strcmp(type, 'ceres')
        grid.dim2.lon = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lon'));
        grid.dim2.lat = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/CERES_EBAF_Ed4.1_Subset_%s.ymonmean.nc', type, par.(type).yr_span), 'lat'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type), 'grid')
    end

    % save grid
end

function read_lfrac(type, par)
% land fraction
    if strcmp(type, 'erai')
        sftlf = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/lmask/interim_lmask.nc', type), 'lsm'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/sftlf.mat', type), 'sftlf');
    elseif strcmp(type, 'gcm')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', par.model));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'sftlf'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    elseif strcmp(type, 'echam')
        file=dir(sprintf('/project2/tas1/CMIP5_piControl/%s/sftlf_*.nc', 'MPI-ESM-LR'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        sftlf=double(ncread(fullpath, 'sftlf'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='sftlf.mat';
        save(sprintf('%s/%s', newdir, filename), 'sftlf');
    end
end

%% define main functions
function diff_struct(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    prefix2 = make_prefix(type, par, par.clim2);

    for vt = {'rad', 'hydro', 'stf', 'srfc', 'ta_si'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix1, vtype)); var1=tmp.(vtype);
        if ~strcmp(vtype, 'ta_si')
            var_vars=tmp.(sprintf('%s_vars',vtype));
        end
        clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix2, vtype)); var2=tmp.(vtype); clear tmp;

        for fn=fieldnames(var1)'; fname=fn{1};
            var.(fname) = var2.(fname)-var1.(fname);
        end

        savedir = make_savedir(type, par);

        if strcmp(vtype, 'rad');       rad = var;   rad_vars=var_vars;   save(sprintf('%s/rad.mat', savedir), 'rad', 'rad_vars');
        elseif strcmp(vtype, 'hydro'); hydro = var; hydro_vars=var_vars; save(sprintf('%s/hydro.mat', savedir), 'hydro', 'hydro_vars');
        elseif strcmp(vtype, 'stf');   stf = var;   stf_vars=var_vars;   save(sprintf('%s/stf.mat', savedir), 'stf', 'stf_vars');
        elseif strcmp(vtype, 'srfc');  srfc = var;  srfc_vars=var_vars;  save(sprintf('%s/srfc.mat', savedir), 'srfc', 'srfc_vars');
        elseif strcmp(vtype, 'ta_si');  ta_si = var; save(sprintf('%s/ta_si.mat', savedir), 'ta_si');
        end

    end
end

function diff_array(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    prefix2 = make_prefix(type, par, par.clim2);

    % for vt = {'zgsi', 'pa_si', 'tai', 'dtdzsi', 'malrsi', 'ma_si_0.95'}; vtype = vt{1}; disp(vtype);
    for vt = {'zgsi', 'pa_si'}; vtype = vt{1}; disp(vtype);
    % for vt = {'ma_si_0.95'}; vtype = vt{1}; disp(vtype);
        if strcmp(vtype, 'malrsi')
            tmp=load(sprintf('%s/%s.mat', prefix1, vtype)); var1=tmp.dtmdzsi; clear tmp;
            tmp=load(sprintf('%s/%s.mat', prefix2, vtype)); var2=tmp.dtmdzsi; clear tmp;
        elseif contains(vtype, 'ma_si')
            tmp=load(sprintf('%s/%s.mat', prefix1, vtype)); var1=tmp.ma_si; clear tmp;
            tmp=load(sprintf('%s/%s.mat', prefix2, vtype)); var2=tmp.ma_si; clear tmp;
        else
            tmp=load(sprintf('%s/%s.mat', prefix1, vtype)); var1=tmp.(vtype); clear tmp;
            tmp=load(sprintf('%s/%s.mat', prefix2, vtype)); var2=tmp.(vtype); clear tmp;
        end

        var = var2-var1;

        savedir = make_savedir(type, par);

        if strcmp(vtype, 'zgsi');  zgsi = var; save(sprintf('%s/zgsi.mat', savedir), 'zgsi');
        elseif strcmp(vtype, 'pa_si');  pa_si = var; save(sprintf('%s/pa_si.mat', savedir), 'pa_si');
        elseif strcmp(vtype, 'tai');  tai = var; save(sprintf('%s/tai.mat', savedir), 'tai');
        elseif strcmp(vtype, 'dtdzsi');  dtdzsi = var; save(sprintf('%s/dtdzsi.mat', savedir), 'dtdzsi');
        elseif strcmp(vtype, 'malrsi');  dtmdzsi = var; save(sprintf('%s/malrsi.mat', savedir), 'dtmdzsi');
        elseif strcmp(vtype, 'ma_si_0.95');  ma_si = var; save(sprintf('%s/ma_si_0.95.mat', savedir), 'ma_si');
        end

    end
end

function diff_flux_z(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    load(sprintf('%s/grid.mat', prefix1)) % read grid data
    lat = grid.dim3.lat;

    prefix_proc1 = make_prefix_proc(type, par, par.clim1);
    prefix_proc2 = make_prefix_proc(type, par, par.clim2);

    for vt = {'flux_z'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix_proc1, vtype)); var1=tmp.(vtype); clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix_proc2, vtype)); var2=tmp.(vtype); clear tmp;

        for l = {'lo'}; land=l{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_z.(land).(fname) = var2.(land).(fname)-var1.(land).(fname);
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_z.(land).(fname).(fw) = var2.(land).(fname).(fw)-var1.(land).(fname).(fw);
                end
            end
        end

        savedir = make_savedir_proc(type, par);

        save(sprintf('%s/flux_z.mat', savedir), 'flux_z', 'lat');

    end
end

function diff_flux_t(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    load(sprintf('%s/grid.mat', prefix1)) % read grid data
    lat = grid.dim3.lat;

    prefix_proc1 = make_prefix_proc(type, par, par.clim1);
    prefix_proc2 = make_prefix_proc(type, par, par.clim2);

    for vt = {'flux_t'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix_proc1, vtype)); var1=tmp.(vtype); clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix_proc2, vtype)); var2=tmp.(vtype); clear tmp;

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_t.(land).(time).(fname) = var2.(land).(time).(fname) - var1.(land).(time).(fname);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_t.(land).(time).(fname).(fw) = var2.(land).(time).(fname).(fw) - var1.(land).(time).(fname).(fw);
                    end
                end
            end
        end

        savedir = make_savedir_proc(type, par);

        save(sprintf('%s/flux_t.mat', savedir), 'flux_t', 'lat');

    end
end

function diff_flux_zt(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    load(sprintf('%s/grid.mat', prefix1)) % read grid data
    lat = grid.dim3.lat;

    prefix_proc1 = make_prefix_proc(type, par, par.clim1);
    prefix_proc2 = make_prefix_proc(type, par, par.clim2);

    for vt = {'flux_zt'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix_proc1, vtype)); var1=tmp.(vtype); clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix_proc2, vtype)); var2=tmp.(vtype); clear tmp;

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_zt.(land).(time).(fname) = var2.(land).(time).(fname) - var1.(land).(time).(fname);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_zt.(land).(time).(fname).(fw) = var2.(land).(time).(fname).(fw) - var1.(land).(time).(fname).(fw);
                    end
                end
            end
        end

        savedir = make_savedir_proc(type, par);

        save(sprintf('%s/flux_zt.mat', savedir), 'flux_zt', 'lat');

    end
end

function diff_vh(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    load(sprintf('%s/grid.mat', prefix1)) % read grid data
    lat = grid.dim3.lat;

    prefix_proc1 = make_prefix_proc(type, par, par.clim1);
    prefix_proc2 = make_prefix_proc(type, par, par.clim2);

    for vt = {'vh'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix_proc1, vtype)); var1=tmp.(vtype); clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix_proc2, vtype)); var2=tmp.(vtype); clear tmp;

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    vh.(land).(time).(fw) = var2.(land).(time).(fw) - var1.(land).(time).(fw);
                end
            end
        end

        savedir = make_savedir_proc(type, par);

        save(sprintf('%s/vh.mat', savedir), 'vh', 'lat');

    end
end

function diff_vh_mon(type, par)
    prefix1 = make_prefix(type, par, par.clim1);
    load(sprintf('%s/grid.mat', prefix1)) % read grid data
    lat = grid.dim3.lat;

    prefix_proc1 = make_prefix_proc(type, par, par.clim1);
    prefix_proc2 = make_prefix_proc(type, par, par.clim2);

    for vt = {'vh_mon'}; vtype = vt{1}; disp(vtype);
        tmp=load(sprintf('%s/%s.mat', prefix_proc1, vtype)); var1=tmp.(vtype); clear tmp;
        tmp=load(sprintf('%s/%s.mat', prefix_proc2, vtype)); var2=tmp.(vtype); clear tmp;

        for l = {'lo'}; land=l{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                vh_mon.(land).(fw) = var2.(land).(fw) - var1.(land).(fw);
            end
        end

        savedir = make_savedir_proc(type, par);

        save(sprintf('%s/vh_mon.mat', savedir), 'vh_mon', 'lat');

    end
end

%% helper functions
% function prefix = make_prefix(type, par, clim)
%     if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55', 'echam_ml', 'echam_pl'}))
%         prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
%     elseif strcmp(type, 'gcm')
%         prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, clim);
%     elseif strcmp(type, 'echam')
%         prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type, clim);
%     end
% end

% function prefix_proc = make_prefix_proc(type, par, clim)
%     if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55', 'echam_ml', 'echam_pl'}))
%         prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.lat_interp);
%     elseif strcmp(type, 'gcm')
%         prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s', type, par.model, clim, par.lat_interp);
%     elseif strcmp(type, 'echam')
%         prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, clim, par.lat_interp);
%     end
% end

% function savedir = make_savedir(type, par)
%     if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
%     elseif strcmp(type, 'gcm')
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
%     elseif strcmp(type, 'echam')
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
%     end
% end

% function savedir = make_savedir_proc(type, par)
%     if any(strcmp(type, {'era5', 'era5c', 'erai', 'merra2', 'jra55'}))
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
%     elseif strcmp(type, 'gcm')
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
%     elseif strcmp(type, 'echam')
%         savedir = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
%     end
%     if ~exist(savedir, 'dir'); mkdir(savedir); end
% end
