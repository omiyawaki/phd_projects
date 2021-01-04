function proc_malr_mon_lat(type, par)
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/%s/ta_mon_lat.mat', prefix_proc, par.lat_interp)); % load taerature data
    load(sprintf('%s/%s/ta_lon_lat.mat', prefix_proc, par.lat_interp)); % load taerature data
    % load(sprintf('%s/%s/masks.mat', prefix_proc, par.lat_interp)); % load land and ocean masks

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    pb = CmdLineProgressBar("Calculating moist adiabatic lapse rate...");
    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        for ilat = 1:length(lat);
            pb.print(ilat, length(lat)); % output progress of moist adiabat calculation
            for imon = 1:12;
                tai = ta.(land);
                pa = 1e-2*grid.dim3.plev; pa = repmat(pa, [1 size(tai, 1) size(tai, 2)]); pa = permute(pa, [2 3 1]);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

                tai = taz.(land);
                pa = 1e-2*paz.(land);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz_z.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

                tai = tasi.(land);
                pa = 1e-2*pasi.(land);
                L = 2510 - 2.38*(tai-T0);
                es = es0*exp(eps*L/R.*(1/T0-1./tai));
                desdT = eps*L.*es./(R*tai.^2);
                dtmdz_si.(land) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                clear tai pa L es desdT;

            end
            for ilon = 1:length(grid.dim3.lon);
                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};

                    tai = ta_t.(land).(time);
                    pa = 1e-2*grid.dim3.plev; pa = repmat(pa, [1 size(tai, 1) size(tai, 2)]); pa = permute(pa, [2 3 1]);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                    tai = taz_t.(land).(time);
                    pa = 1e-2*paz_t.(land).(time);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_z_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                    tai = tasi_t.(land).(time);
                    pa = 1e-2*pasi_t.(land).(time);
                    L = 2510 - 2.38*(tai-T0);
                    es = es0*exp(eps*L/R.*(1/T0-1./tai));
                    desdT = eps*L.*es./(R*tai.^2);
                    dtmdz_si_t.(land).(time) = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));
                    clear tai pa L es desdT;

                end
            end
        end
    end

    % save data into mat file
    printname = [foldername 'malr_mon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtmdz', 'dtmdz_z', 'dtmdz_si');

    printname = [foldername 'malr_lon_lat.mat'];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(printname, 'dtmdz_t', 'dtmdz_z_t', 'dtmdz_si_t');

end % compute mon x lat moist adiabat field
