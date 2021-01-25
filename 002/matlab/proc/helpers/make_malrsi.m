function make_malrsi(type, par)
% compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2')
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
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
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/tai.mat', prefix)); clear zgi; % read interpolated temperature

    pa = repmat(par.pa', [1 length(grid.dim3.lon) length(grid.dim3.lat) 12]);
    pa = permute(pa, [2 3 1 4]);
    pa = pa/1e2; % convert to mb

    dtmdz = nan([length(grid.dim3.lon) length(grid.dim3.lat) length(par.pa) 12]);

    % calculate MALR following equations 3-5 in Stone and Carlson (1979)
    dalr = 9.8; % K/km
    eps = 0.622;
    R = 0.287; % J/g/K
    cp = 1.005; % J/g/K
    es0 = 6.11; % mb
    T0 = 273; % K
    L = 2510 - 2.38*(tai-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tai));
    desdT = eps*L.*es./(R*tai.^2);

    dtmdz = dalr * (1+eps*L.*es./(pa.*R.*tai))./(1+(eps.*L./(cp.*pa)).*(desdT));

    clear L es desdT;
    % 2 m dtmdz
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        tas = srfc.t2m;
        ps = srfc.sp;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
        ps = srfc.PS;
    elseif strcmp(type, 'gcm')
        tas = srfc.tas;
        ps = srfc.ps;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
        ps = srfc.aps;
    end
    L = 2510 - 2.38*(tas-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tas));
    desdT = eps*L.*es./(R*tas.^2);
    dtasmdz = dalr * (1+eps*L.*es./(1e-2*ps.*R.*tas))./(1+(eps.*L./(cp.*1e-2*ps)).*(desdT));

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.sp)]), [2 3 1 4]));
    elseif strcmp(type, 'merra2')
        ps_vert = repmat(srfc.PS, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(par.pa', [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.ps)]), [2 3 1 4]);
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 size(dtmdz, 3)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(par.pa', [1 size(srfc.aps)]), [2 3 1 4]);
    end

    pa = permute(pa, [3 1 2 4]); % bring height to 1st
    dtmdz = permute(dtmdz, [3 1 2 4]); % bring height to 1st
    pb = CmdLineProgressBar("Sorting and interpolating dtmdz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)
                dtmdzsi(:,lo,la,mo) = interp1(pa(:,lo,la,mo)/ps_vert(lo,la,1,mo), dtmdz(:,lo,la,mo), grid.dim3.si);
            end
        end
    end
    dtmdzsi(1,:,:,:) = dtasmdz;
    dtmdzsi = permute(dtmdzsi, [2 3 1 4]); % bring height to 3rd

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/gcm/%s/%s', par.model, par.gcm.clim);
    elseif strcmp(type, 'echam'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'})); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='malrsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzsi', '-v7.3');

end
