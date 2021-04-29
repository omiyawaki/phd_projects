function make_malrsi(type, par)
% compute moist adiabatic lapse rate at every lat, time
% calculate moist adiabats
    
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    newdir = make_savedir(type, par);
    
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface variable data
    load(sprintf('%s/tai_simp.mat', prefix)); clear zgi; % read interpolated temperature

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
    tas = rename_tas(type, srfc);
    ps = rename_ps(type, srfc);
    
    L = 2510 - 2.38*(tas-T0);
    es = es0*exp(eps*L/R.*(1/T0-1./tas));
    desdT = eps*L.*es./(R*tas.^2);
    dtasmdz = dalr * (1+eps*L.*es./(1e-2*ps.*R.*tas))./(1+(eps.*L./(cp.*1e-2*ps)).*(desdT));

    [ps_vert, pa] = make_surfmask_vars(grid, type, srfc, par);

    pa = permute(pa, [3 1 2 4]); % bring height to 1st
    dtmdz = permute(dtmdz, [3 1 2 4]); % bring height to 1st
    pb = CmdLineProgressBar("Sorting and interpolating dtmdz to new standard grid...");
    for lo=1:size(pa,2)
        pb.print(lo, size(pa,2));
        for la=1:size(pa,3)
            for mo=1:size(pa,4)

                local_si = squeeze(pa(:,lo,la,mo)/ps_vert(lo,la,1,mo));
                local_lr = squeeze(dtmdz(:,lo,la,mo));

                % remove missing data
                idx_nan = isnan(local_si) | isnan(local_lr);
                local_si(idx_nan) = [];
                local_lr(idx_nan) = [];

                % insert surface lapse rate
                if ~any(ismember(local_si, 1))
                    local_si = [1; local_si];
                    local_lr = [squeeze(dtasmdz(lo,la,mo)); local_lr];
                end

                % interpolate to standard sigma grid
                dtmdzsi(:,lo,la,mo) = interp1(local_si, local_lr, grid.dim3.si);

                clear local_si local_lr idx_nan

            end
        end
    end

    dtmdzsi = permute(dtmdzsi, [2 3 1 4]); % bring height to 3rd

    filename='malrsi.mat';
    save(sprintf('%s/%s', newdir, filename), 'dtmdzsi', '-v7.3');

end
