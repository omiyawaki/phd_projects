% This script imports all CMIP data and creates the multi-model mean
clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

gcm_info

%% set parameters
if 1
% par.erai.yr_span = '2000_2012'; % spanning years for ERA-Interim
par.gcm.clim = 'hist-pi'; % choose either piControl or abrupt4xCO2
par.outname = 'mmm'; % use mmm for standard multimodel mean (all models available) and mmm_subset where subset=name of climate where you are taking the mmm of the subset of models in common with the subset climate
par.lat_interp = '1.00'; % latitudinal grid spacing to interpolate to (deg)
par.lat = -90:str2num(par.lat_interp):90; % define standard latitude grid for 'std' interpolation
par.lon_interp = '1.00'; % longitudinal grid spacing to interpolone to (deg)
par.lon = -90:str2num(par.lon_interp):90; % define standard lonitude grid for 'std' interpolonion
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold for RCE definition. RCE is defined as where abs(R1) < ep
par.ga_swp = 0.9; % optional threshold for RAE. If undefined, the default value is 1-par.ep
par.ma_init = 0.95; % initial starting level for moist adiabat ('surf' = start with dry adiabat until LCL or enter sigma level for starting saturated adiabat)
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.si = 1e-5*par.pa; % high resolution vertical grid to interpolate to
par.si_eval = [0.8 0.85 0.9]; % sigma level for evaluating inversion strength (T(si_eval) - T(surface))
par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up = 0.4; % sigma level for upper boundary of vertical average for close to moist adiabatic
par.z = [0:500:40e3]';
par.gcm.fw = {'mse', 'dse'};
par.cpd = 1005.7; par.cpv = 1870; par.cpl = 4186; par.cpi = 2108; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
end

type = 'gcm';
make_grid(type, par);
% choose_mmm_lat_mon_silev(type, par); % make mmm of mon x lat x lev data
choose_mmm_lat_mon(type, par); % make mmm of mon x lat data
choose_mmm_lon_lat(type, par); % make mmm of lon x lat data
choose_mmm_lat(type, par); % make mmm of lat data

function make_grid(type, par)
    filename = 'grid.mat';
    grid.dim2.lon = par.lon;
    grid.dim3.lon = par.lon;
    grid.dim2.lat = par.lat;
    grid.dim3.lat = par.lat;
    grid.dim3.plev = par.pa;
    grid.dim3.z = par.z;
    grid.dim3.si = 1e-5*par.pa;
    newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.outname, par.gcm.clim);
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    save(sprintf('%s/%s', newdir, filename), 'grid');
end

%% Choose functions to run
function choose_mmm_lat_mon_silev(type, par)
    mmm_ta_mon_lat(type, par);
    mmm_ma_mon_lat(type, par);
end

function choose_mmm_lat_mon(type, par)
    % mmm_flux_z(type, par);
    mmm_vh_mon(type, par);
end

function choose_mmm_lon_lat(type, par)
    % mmm_flux_t(type, par);
end

function choose_mmm_lat(type, par)
    % mmm_flux_zt(type, par);
    mmm_vh(type, par);
end

%% Mon x lat x lev functions
function mmm_ta_mon_lat(type, par)

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        tasi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        tasi0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, 'native'));

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            % interpolate native grid data to standard grid
            tasi0i = interp1(grid0.grid.dim3.lat, tasi0.tasi.(land), grid.dim3.lat);
            tasi_mmm.(land) = nanmean(cat(4,tasi0i,tasi_mmm.(land)),4);
        end % land

    end % models

    tasi = tasi_mmm;

    printname = [foldername 'ta_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'tasi', 'lat', '-v7.3');

end
function mmm_ma_mon_lat(type, par)

    lat = par.lat;
    for l = {'lo'}; land=l{1};
        masi_mmm.(land) = nan(length(par.lat), 12, length(par.si));
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, 'native');
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        masi0=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_mon_lat.mat', type, par.model, par.gcm.clim, 'native'));

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            % interpolate native grid data to smandard grid
            masi0i = interp1(grid0.grid.dim3.lat, masi0.masi.(land), grid.dim3.lat);
            masi_mmm.(land) = nanmean(cat(4,masi0i,masi_mmm.(land)),4);
        end % land

    end % models

    masi = masi_mmm;

    printname = [foldername 'ma_mon_lat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'masi', 'lat', '-v7.3');

end

%% Mon x lat functions

function mmm_flux_z(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        for fn = var_vec; fname = fn{1};
            flux_z_mmm.(land).(fname) = nan(length(par.lat), 12);
        end
        for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                flux_z_mmm.(land).(fname).(fw) = nan(length(par.lat), 12);
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        flux_z0 = load(sprintf('%s/%s/flux_z.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_z0i.flux_z.(land).(fname) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname), grid.dim3.lat);
                flux_z_mmm.(land).(fname) = nanmean(cat(3, flux_z0i.flux_z.(land).(fname), flux_z_mmm.(land).(fname)), 3);
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_z0i.flux_z.(land).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_z0.flux_z.(land).(fname).(fw), grid.dim3.lat);
                    flux_z_mmm.(land).(fname).(fw) = nanmean(cat(3, flux_z0i.flux_z.(land).(fname).(fw), flux_z_mmm.(land).(fname).(fw)), 3);
                end
            end
        end
    end

    flux_z = flux_z_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_z'), 'flux_z', 'lat', '-v7.3');

end

function mmm_vh_mon(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        f_vec = par.gcm.fw;
        for f = f_vec; fw = f{1};
            vh_mon_mmm.(land).(fw) = nan(length(par.lat), 12);
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        vh_mon0 = load(sprintf('%s/%s/vh_mon.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                vh_mon0i.vh_mon.(land).(fw) = interp1(grid0.grid.dim3.lat, vh_mon0.vh_mon.(land).(fw), grid.dim3.lat);
                vh_mon_mmm.(land).(fw) = nanmean(cat(3, vh_mon0i.vh_mon.(land).(fw), vh_mon_mmm.(land).(fw)), 3);
            end
        end % land
    end % models

    vh_mon = vh_mon_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh_mon'), 'vh_mon', 'lat', '-v7.3');

end

%% Lon x lat functions
function mmm_flux_t(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_t_mmm.(land).(time).(fname) = nan(length(par.lat), length(par.lon));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_t_mmm.(land).(time).(fname).(fw) = nan(length(par.lat), length(par.lon));
                end
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        flux_t0 = load(sprintf('%s/%s/flux_t.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_t0i.flux_t.(land).(time).(fname) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname), grid.dim3.lon);
                    flux_t0i.flux_t.(land).(time).(fname) = permute(flux_t0i.flux_t.(land).(time).(fname), [2 1]);
                    flux_t0i.flux_t.(land).(time).(fname) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname), grid.dim3.lat);
                    flux_t0i.flux_t.(land).(time).(fname) = permute(flux_t0i.flux_t.(land).(time).(fname), [2 1]);
                    flux_t_mmm.(land).(time).(fname) = nanmean(cat(3, flux_t0i.flux_t.(land).(time).(fname), flux_t_mmm.(land).(time).(fname)), 3);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lon, flux_t0.flux_t.(land).(time).(fname).(fw), grid.dim3.lon);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_t0i.flux_t.(land).(time).(fname).(fw), grid.dim3.lat);
                        flux_t0i.flux_t.(land).(time).(fname).(fw) = permute(flux_t0i.flux_t.(land).(time).(fname).(fw), [2 1]);
                        flux_t_mmm.(land).(time).(fname).(fw) = nanmean(cat(3, flux_t0i.flux_t.(land).(time).(fname).(fw), flux_t_mmm.(land).(time).(fname).(fw)), 3);
                    end
                end
            end % time
        end % land
    end % models

    flux_t = flux_t_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_t'), 'flux_t', 'lat', '-v7.3');

end

%% Lat functions
function mmm_flux_zt(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
            for fn = var_vec; fname = fn{1};
                flux_zt_mmm.(land).(time).(fname) = nan(1, length(par.lat));
            end
            for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    flux_zt_mmm.(land).(time).(fname).(fw) = nan(1, length(par.lat));
                end
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        flux_zt0 = load(sprintf('%s/%s/flux_zt.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
                for fn = var_vec; fname = fn{1};
                    flux_zt0i.flux_zt.(land).(time).(fname) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname), grid.dim3.lat);
                    flux_zt_mmm.(land).(time).(fname) = nanmean(cat(1, flux_zt0i.flux_zt.(land).(time).(fname), flux_zt_mmm.(land).(time).(fname)), 1);
                end
                for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
                    f_vec = par.gcm.fw;
                    for f = f_vec; fw = f{1};
                        flux_zt0i.flux_zt.(land).(time).(fname).(fw) = interp1(grid0.grid.dim3.lat, flux_zt0.flux_zt.(land).(time).(fname).(fw), grid.dim3.lat);
                        flux_zt_mmm.(land).(time).(fname).(fw) = nanmean(cat(1, flux_zt0i.flux_zt.(land).(time).(fname).(fw), flux_zt_mmm.(land).(time).(fname).(fw)), 1);
                    end
                end
            end % time
        end % land
    end % models

    flux_zt = flux_zt_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'flux_zt'), 'flux_zt', 'lat', '-v7.3');

end

function mmm_vh(type, par)
    lat = par.lat;

    for l = {'lo'}; land=l{1};
        for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
            f_vec = par.gcm.fw;
            for f = f_vec; fw = f{1};
                vh_mmm.(land).(time).(fw) = nan(1,length(par.lat));
            end
        end
    end

    pb = CmdLineProgressBar("Creating the multi-model mean...");
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        pb.print(k, length(par.gcm_models)); % output progress of moist adiabat calculonion

        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        grid0 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        vh0 = load(sprintf('%s/%s/vh.mat', prefix_proc, 'native')); % load lat x mon RCAE data

        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.outname, par.gcm.clim, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.outname, par.gcm.clim));

        for l = {'lo'}; land=l{1};
            for t = {'ann', 'djf', 'mam', 'jja', 'son'}; time = t{1};
                f_vec = par.gcm.fw;
                for f = f_vec; fw = f{1};
                    vh0i.vh.(land).(time).(fw) = interp1(grid0.grid.dim3.lat, vh0.vh.(land).(time).(fw), grid.dim3.lat);
                    vh_mmm.(land).(time).(fw) = nanmean(cat(1, vh0i.vh.(land).(time).(fw), vh_mmm.(land).(time).(fw)), 1);
                end
            end % time
        end % land
    end % models

    vh = vh_mmm;

    % save energy flux data into mat file
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(sprintf('%s%s', foldername, 'vh'), 'vh', 'lat', '-v7.3');

end
