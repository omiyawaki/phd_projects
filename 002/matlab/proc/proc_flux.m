function proc_flux(type, par)
    if any(strcmp(type, {'era5', 'era5c' 'erai'}))
        don = load(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology')); % read donohoe data
        prefix_ceres=sprintf('/project2/tas1/miyawaki/projects/002/data/read/ceres'); % prefix for CERES data
        % load(sprintf('%s/div.mat', prefix)) % read divergence data
        % load(sprintf('%s/tend.mat', prefix)) % read tendency data
        prefix_don=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'erai');
        if any(strcmp(type, {'erai', 'era5', 'era5c'})) & strcmp(par.(type).yr_span, '1979_2018')
            don79 = load(sprintf('%s/dondiv79.mat', prefix_don)); % read Donohoe data 1979--2018
        elseif any(strcmp(type, {'erai', 'era5', 'era5c'})) & strcmp(par.(type).yr_span, '2000_2018')
            don79 = load(sprintf('%s/dondiv00.mat', prefix_don)); % read Donohoe data 1979--2018
        end
    end

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);

    load(sprintf('%s/grid.mat', prefix)) % read grid data
    load(sprintf('%s/rad.mat', prefix)) % read radiation data
    load(sprintf('%s/hydro.mat', prefix)) % read hydrology data
    load(sprintf('%s/stf.mat', prefix)) % read surface turbulent flux data
    % load(sprintf('%s/masks.mat', prefix_proc)); % load land and ocean masks

    if strcmp(par.lat_interp, 'std')
        lat = par.lat_std;
    else
        lat = grid.dim3.lat;
    end

    % interpolate onto std lat x lon grid
    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        ceres = struct();
        tmp = load(sprintf('%s/grid.mat', prefix_ceres)); ceres.grid = tmp.grid; % read grid data
        tmp = load(sprintf('%s/rad.mat', prefix_ceres)); ceres.rad = tmp.rad; % read CERES radiation data
        for fn = {'TEDIV', 'TETEN'}; fname = fn{1};
            flux.(fname) = permute(don.(fname), [2 3 1]); % bring lat forward
            flux.(fname) = interp1(don.lat, flux.(fname), lat, 'linear');
            flux.(fname) = permute(flux.(fname), [2 1 3]); % bring lon forward
            flux.(fname) = interp1(don.lon, flux.(fname), grid.dim2.lon, 'linear');
        end
        % flux.don79div = interp1(don79.lat, don79.dondiv, lat, 'linear'); % interpolate to standard lat
        % flux.don79div = repmat(flux.don79div, [1 1 length(grid.dim2.lon)]); % repeat to longitude because don79 data is already zonally averaged
        % flux.don79div = permute(flux.don79div, [3 1 2]); % bring lon to 1st
        for fn = fieldnames(ceres.rad)'; fname = fn{1}; % interpolate to std lat and ERA lon
            ceres.(fname) = interp1(ceres.grid.dim2.lon, ceres.rad.(fname), grid.dim2.lon, 'linear');
            ceres.(fname) = permute(ceres.(fname), [2 1 3]); % bring lat forward
            ceres.(fname) = interp1(ceres.grid.dim2.lat, ceres.(fname), lat, 'linear');
            ceres.(fname) = permute(ceres.(fname), [2 1 3]); % bring lon forward
        end
        ceres.tsr = ceres.tsdr - ceres.tsur; ceres.net = ceres.tsr - ceres.ttr;
        ceres.str = -ceres.str;
        ceres.ra = ceres.tsr - ceres.ssr + ceres.str - ceres.ttr;
    end
    for fn = rad_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(rad.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    for fn = hydro_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(hydro.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    for fn = stf_vars; fname = fn{1}; % interpolate to std lat
        flux.(fname) = permute(stf.(fname), [2 1 3]);
        flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
        flux.(fname) = permute(flux.(fname), [2 1 3]);
    end
    if contains(type, 'era') || any(strcmp(type, {'gcm'}))
        load(sprintf('%s/tend.mat', prefix)) % read surface turbulent flux data
        flux.tend = permute(tend.tend, [2 1 3]);
        flux.tend = interp1(grid.dim2.lat, flux.tend, lat, 'linear');
        flux.tend = permute(flux.tend, [2 1 3]);
    end
    % if any(strcmp(type, {'era5', 'era5c' 'erai'}))
    %     for fn = tend_vars_txt; fname = fn{1}; % interpolate to std lat
    %         flux.(fname) = permute(tend.(fname), [2 1 3]);
    %         flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
    %         flux.(fname) = permute(flux.(fname), [2 1 3]);
    %     end; clear tend
    %     for fn = div_vars_txt; fname = fn{1}; % interpolate to std lat
    %         flux.(fname) = permute(div.(fname), [2 1 3]);
    %         flux.(fname) = interp1(grid.dim2.lat, flux.(fname), lat, 'linear');
    %         flux.(fname) = permute(flux.(fname), [2 1 3]);
    %     end; clear div
    % end
    % flux.w500 = permute(w500, [2 1 3]);
    % flux.w500 = interp1(grid.dim3.lat, flux.w500, lat, 'linear');
    % flux.w500 = permute(flux.w500, [2 1 3]);
    % flux.vas = permute(vas, [2 1 3]);
    % flux.vas = interp1(grid.dim3.lat, flux.vas, lat, 'linear');
    % flux.vas = permute(flux.vas, [2 1 3]);

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        % compute surface turbulent fluxes directly from INTP data
        % multiply by negative to define flux from surface to atmosphere as positive
        flux.stf.mse = -( flux.sshf + flux.slhf ); flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*(flux.cp+flux.lsp) - flux.sshf;
    elseif strcmp(type, 'merra2')
        flux.stf.mse = flux.HFLUX + flux.EFLUX;
        flux.stf.dse = par.L*flux.PRECTOT + flux.HFLUX;
    elseif strcmp(type, 'hahn')
        flux.stf.mse = flux.LHFLX + flux.SHFLX;
        flux.stf.dse = par.L*(flux.PRECC+flux.PRECL+flux.PRECSC+flux.PRECSL) + flux.SHFLX;
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        flux.stf.mse = flux.hfls + flux.hfss; flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*flux.pr + flux.hfss;
    elseif strcmp(type, 'echam')
        flux.stf.mse = -(flux.ahfl + flux.ahfs); flux.stf.mse2 = flux.stf.mse;
        flux.stf.dse = par.L*(flux.aprc+flux.aprl) - flux.ahfs;
    end
    flux.stf.mse_old = flux.stf.mse;
    flux.stf.mse_ac = flux.stf.mse;
    flux.stf.mse_sc = flux.stf.mse;
    flux.stf.mse_ac_ra = flux.stf.mse;
    flux.stf.mse_sc_ra = flux.stf.mse;
    flux.stf.dse_old = flux.stf.dse;

    f_vec = assign_fw(type, par);
    for f = f_vec; fw = f{1};
        if any(strcmp(type, {'era5', 'era5c', 'erai'}));
            flux.rtoa = flux.tsr + flux.ttr; % net flux at TOA
            flux.olr = flux.ttr;
            flux.swsfc = -flux.ssr;
            flux.lwsfc = -flux.str;
            flux.lw = flux.ttr - flux.str; flux.sw = flux.tsr-flux.ssr; % compute net shortwave and longwave flux through atmosphere
            if contains(fw, 'ceresrad'); flux.ra.(fw) = ceres.ra;  % compute net radiative cooling from radiative fluxes
            else; flux.ra.(fw) = flux.tsr - flux.ssr + flux.ttr - flux.str; end % use radiative cooling from CERES data
        elseif strcmp(type, 'merra2');
            flux.rtoa = flux.SWTNT - flux.LWTUP; % net flux at TOA
            flux.olr = -flux.LWTUP;
            flux.swsfc = -flux.SWGNT;
            flux.lwsfc = -flux.LWGNT;
            flux.lw = -flux.LWTUP - flux.LWGNT; flux.sw = flux.SWTNT - flux.SWGNT;
            flux.ra.(fw) = flux.SWTNT - flux.LWTUP - flux.SWGNT - flux.LWGNT; % radiative cooling
        elseif strcmp(type, 'hahn');
            flux.rtoa = flux.FSNT - flux.FLNT; % net flux at TOA
            flux.olr = -flux.FLNT;
            flux.swsfc = -flux.FSNS;
            flux.lwsfc = flux.FLNS;
            flux.lw = -flux.FLNT + flux.FLNS; flux.sw = flux.FSNT - flux.FSNS;
            flux.ra.(fw) = flux.FSNT - flux.FLNT - flux.FSNS + flux.FLNS; % radiative cooling
        elseif any(strcmp(type, {'gcm', 'jra55'}));
            flux.rtoa = flux.rsdt - flux.rsut - flux.rlut; % net flux at TOA
            flux.olr = -flux.rlut;
            flux.swsfc = flux.rsus - flux.rsds;
            flux.lwsfc = flux.rlus - flux.rlds;
            flux.lw = flux.rlus - flux.rlds - flux.rlut; flux.sw = flux.rsdt - flux.rsut + flux.rsus - flux.rsds;
            flux.ra.(fw) = flux.rsdt - flux.rsut + flux.rsus - flux.rsds + flux.rlus - flux.rlds - flux.rlut;
        elseif strcmp(type, 'echam');
            flux.rtoa = flux.trad0 + flux.srad0; % net flux at TOA
            flux.olr = flux.trad0;
            flux.swsfc = -flux.srads;
            flux.lwsfc = -flux.trads;
            flux.lw = flux.trad0 - flux.trads; flux.sw = flux.srad0 - flux.srads;
            flux.ra.(fw) = flux.lw + flux.sw;
        end % calculate atmospheric radiative cooling
        flux.rsfc = flux.swsfc + flux.lwsfc;

        if any(strcmp(fw, {'mse_old', 'dse_old'}))
            flux.res.(fw) = flux.ra.(fw) + flux.stf.(fw); % infer MSE tendency and flux divergence as residuals
        elseif any(strcmp(fw, {'mse', 'mse_ac', 'mse_sc', 'mse_ac_ra', 'mse_sc_ra', 'dse'}))
            flux.res.(fw) = flux.ra.(fw) + flux.stf.(fw) - flux.tend; % infer MSE tendency and flux divergence as residuals
        elseif any(strcmp(fw, {'mse2'}))
            flux.res.(fw) = flux.lw + flux.stf.(fw);
        elseif strcmp(fw, 'db13')
            flux.res.(fw) = flux.TEDIV + flux.TETEN; % use MSE tendency and flux divergence from DB13
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'db13s')
            flux.res.(fw) = flux.TEDIV; % use MSE flux divergence from DB13, ignore MSE tendency term
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'db13t')
            flux.res.(fw) = flux.TEDIV - flux.tend; % use MSE flux divergence from DB13, use MSE tendency from ERA
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div')
            flux.res.(fw) = flux.divt + flux.divg + flux.divq*par.L; % use MSE tendency and flux divergence from ERA5 output
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'divt')
            flux.res.(fw) = flux.divt + flux.divg + flux.divq*par.L - flux.tend; % use MSE tendency and flux divergence from ERA5 output
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div79')
            flux.res.(fw) = flux.don79div - flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div00')
            flux.res.(fw) = flux.don79div - flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - ceres.ra; % infer the surface turbulent fluxes
        elseif strcmp(fw, 'div00erarad')
            flux.res.(fw) = flux.don79div - flux.tend; % use MSE tendency from ERA output and MSE flux divergence from Donohoe MSE flux transport data
            flux.stf.(fw) = flux.res.(fw) - flux.ra.(fw); % infer the surface turbulent fluxes
        end

        flux.shf.(fw) = flux.lwsfc + flux.stf.(fw); % surface LW and surface turbulent fluxes
        flux.sfc.(fw) = flux.rsfc + flux.stf.(fw); % net flux at surface

        if strcmp(fw, 'mse2')
            flux.r1.(fw) = (flux.res.(fw))./flux.lw; % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.lw; % calculate nondimensional number R2 disregarding MSE budget closure
        else
            flux.r1.(fw) = (flux.res.(fw))./flux.ra.(fw); % calculate nondimensional number R1 disregarding MSE budget closure
            flux.r2.(fw) = flux.stf.(fw)./flux.ra.(fw); % calculate nondimensional number R2 disregarding MSE budget closure
        end
        if any(strcmp(type, {'era5', 'era5c', 'erai'}));
            flux.ftoa.(fw) = flux.tsr + flux.ttr; flux.fsfc.(fw) = -flux.ssr - flux.str + flux.stf.(fw);
        elseif strcmp(type, 'hahn')
            flux.ftoa.(fw) = flux.FSNT - flux.FLNT;
            flux.fsfc.(fw) = -flux.FSNS + flux.FLNS + flux.stf.(fw);
        elseif strcmp(type, 'merra2')
            flux.ftoa.(fw) = flux.SWTNT - flux.LWTUP;
            flux.fsfc.(fw) = -flux.SWGNT - flux.LWGNT + flux.stf.(fw);
        elseif any(strcmp(type, {'gcm', 'jra55'}));
            flux.ftoa.(fw) = flux.rsdt - flux.rsut - flux.rlut;
            flux.fsfc.(fw) = flux.rsus - flux.rsds + flux.rlus - flux.rlds + flux.stf.(fw);
        elseif strcmp(type, 'echam');
            flux.ftoa.(fw) = flux.trad0 + flux.srad0;
            flux.fsfc.(fw) = -flux.trads - flux.srads + flux.stf.(fw);
        end

        % linear decomposition of R1 seasonality
        flux.comp1.(fw) = (flux.res.(fw)-nanmean(flux.res.(fw),3))./nanmean(flux.ra.(fw),3);
        flux.comp2.(fw) = - nanmean(flux.res.(fw),3)./nanmean(flux.ra.(fw),3).^2 .* (flux.ra.(fw)-nanmean(flux.ra.(fw),3));
        % flux.comp2.(fw) = - nanmean(flux.res.(fw)./flux.ra.(fw).^2, 3) .* (flux.ra.(fw)-nanmean(flux.ra.(fw),3));

    end


    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c');
        % var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend', 'divt', 'divg', 'divq', 'TETEN', 'TEDIV', 'don79div'};
        var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        % foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
    elseif strcmp(type, 'hahn')
        var_vec = {'LHFLX', 'SHFLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    elseif strcmp(type, 'merra2')
        var_vec = {'EFLUX', 'HFLUX', 'PRECCON', 'PRECTOT', 'EVAP', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        % foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
    elseif any(strcmp(type, {'jra55'}))
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        % foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.(type).yr_span, par.lat_interp);
    elseif any(strcmp(type, {'gcm'}))
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        % foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        var_vec = {'ahfl', 'ahfs', 'aprc', 'aprl', 'evap', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        % foldername = sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/', type, par.echam.clim, par.lat_interp);
    end

    for fn = var_vec; fname = fn{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); flux_n.(land).(fname) = flux.(fname);
            elseif strcmp(land, 'l'); flux_n.(land).(fname) = flux.(fname) .*mask.ocean;
            elseif strcmp(land, 'o'); flux_n.(land).(fname) = flux.(fname) .*mask.land;
            end

            % take zonal averages
            flux_z.(land).(fname) = squeeze(nanmean(flux_n.(land).(fname), 1));

            % take time averages
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                if strcmp(time, 'ann')
                    flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname), 3);
                elseif strcmp(time, 'djf')
                    flux_shift.(land).(fname) = circshift(flux_n.(land).(fname), 1, 3);
                    flux_t.(land).(time).(fname) = nanmean(flux_shift.(land).(fname)(:,:,1:3), 3);
                elseif strcmp(time, 'jja')
                    flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,6:8), 3);
                elseif strcmp(time, 'mam')
                    flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,3:5), 3);
                elseif strcmp(time, 'son')
                    flux_t.(land).(time).(fname) = nanmean(flux_n.(land).(fname)(:,:,9:11), 3);
                end
                flux_zt.(land).(time).(fname) = squeeze(nanmean(flux_t.(land).(time).(fname), 1));
            end
        end
    end

    for fn = {'ra', 'stf', 'res', 'r1', 'r2', 'ftoa', 'fsfc', 'sfc', 'shf', 'comp1', 'comp2'}; fname = fn{1};
        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                if strcmp(land, 'lo'); flux_n.(land).(fname).(fw) = flux.(fname).(fw);
                elseif strcmp(land, 'l'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.ocean;
                elseif strcmp(land, 'o'); flux_n.(land).(fname).(fw) = flux.(fname).(fw) .*mask.land;
                end

                if strcmp(fname, 'res')

                    % compute northward MSE transport using the residual data
                    tediv_0 = fillmissing(flux_n.(land).res.(fw), 'constant', 0); % replace nans with 0s so missing data doesn't influence transport
                    % tediv_z = squeeze(trapz(deg2rad(grid.dim2.lon), tediv_0, 1)); % zonally integrate
                    tediv_z = squeeze(nanmean(tediv_0, 1)); % zonally average

                    % integrate from south pole
                    if grid.dim3.lat(1) > 0
                        vh_mon.(land).(fw) = cumtrapz(flip(deg2rad(lat),1), flip(2*pi*par.a^2*cosd(lat).*tediv_z,1), 1); % cumulatively integrate
                        vh_mon.(land).(fw) = flip(vh_mon.(land).(fw),1);
                    else
                        vh_mon.(land).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_z, 1); % cumulatively integrate
                    end

                    % remove annual mean global mean inbalance
                    if any(strcmp(fw, {'mse_ac', 'mse_ac_ra', 'mse_sc', 'mse_sc_ra'}))
                        if contains(fw, 'mse_sc')
                            % remove global mean from flux divergence
                            resmean = permute(vh_mon.(land).(fw)(end,:), [2 1]);
                            resmean = repmat(resmean, [1 length(grid.dim3.lon) length(grid.dim3.lat)]);
                            resmean = permute(resmean, [2 3 1]);
                        elseif contains(fw, 'mse_ac')
                            vh_ann = nanmean(vh_mon.(land).(fw),2);
                            disp(sprintf('Annual average imbalance = %g Wm^-2', vh_ann(end)/(4*pi*par.a^2)))
                            resmean = vh_ann(end)*ones(length(grid.dim3.lon), length(grid.dim3.lat), 12);
                        end

                        if contains(fw, '_ra')
                            flux_n.(land).res.(fw) = flux_n.(land).res.(fw) - resmean./(4*pi*par.a^2);
                            flux_n.(land).ra.(fw) = flux_n.(land).ra.(fw) - resmean./(4*pi*par.a^2);
                            [flux_z, flux_t, flux_zt] = comp_flux(flux_z, flux_t, flux_zt, flux_n, land, 'ra', fw);
                        else
                            flux_n.(land).res.(fw) = flux_n.(land).res.(fw) - resmean./(4*pi*par.a^2);
                        end

                        % re-compute northward MSE transport using the residual data
                        tediv_0 = fillmissing(flux_n.(land).res.(fw), 'constant', 0); % replace nans with 0s so missing data doesn't influence transport
                        % tediv_z = squeeze(trapz(deg2rad(grid.dim2.lon), tediv_0, 1)); % zonally integrate
                        tediv_z = squeeze(nanmean(tediv_0, 1)); % zonally average

                        if grid.dim3.lat(1) > 0
                            vh_mon.(land).(fw) = cumtrapz(flip(deg2rad(lat),1), flip(2*pi*par.a^2*cosd(lat).*tediv_z,1), 1); % cumulatively integrate
                            vh_mon.(land).(fw) = flip(vh_mon.(land).(fw),1);
                        else
                            vh_mon.(land).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_z, 1); % cumulatively integrate
                        end
                    end

                    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                        if strcmp(time, 'ann')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw), 3)); % take time average
                        elseif strcmp(time, 'djf')
                            flux_n_shift.(land).res = circshift(flux_n.(land).res.(fw), 1, 3); % shift month by 1 to allow for djf average
                            tediv_t = squeeze(nanmean(flux_n_shift.(land).res(:,:,1:3), 3)); % take time average
                        elseif strcmp(time, 'jja')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,6:8), 3)); % take time average
                        elseif strcmp(time, 'mam')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,3:5), 3)); % take time average
                        elseif strcmp(time, 'son')
                            tediv_t = squeeze(nanmean(flux_n.(land).res.(fw)(:,:,9:11), 3)); % take time average
                        end
                        % tediv_tz = trapz(deg2rad(grid.dim2.lon), tediv_t, 1); % zonally integrate
                        tediv_tz = squeeze(nanmean(tediv_t, 1)); % zonally integrate
                        tediv_tz = fillmissing(tediv_tz, 'constant', 0);

                        if grid.dim3.lat(1) > 0
                            vh.(land).(time).(fw) = cumtrapz(flip(deg2rad(lat),1), flip(2*pi*par.a^2*cosd(lat).*tediv_tz',1),1); % cumulatively integrate in latitude
                            vh.(land).(time).(fw) = flip(vh.(land).(time).(fw),1);
                        else
                            vh.(land).(time).(fw) = cumtrapz(deg2rad(lat), 2*pi*par.a^2*cosd(lat).*tediv_tz',1); % cumulatively integrate in latitude
                        end
                    end
                end

                [flux_z, flux_t, flux_zt] = comp_flux(flux_z, flux_t, flux_zt, flux_n, land, fname, fw);

            end
        end
    end

    % save energy flux data into mat file
    foldername = make_savedir_proc(type, par);
    % if ~exist(foldername, 'dir')
    %     mkdir(foldername)
    % end
    for v = {'flux', 'flux_t', 'flux_z', 'flux_zt', 'vh', 'vh_mon'}; varname = v{1}; % removed flux and flux_t because it takes forever
        save(sprintf('%s%s', foldername, varname), varname, 'lat', '-v7.3');
    end

end % process radiative fluxes into one struct

function [flux_z, flux_t, flux_zt] = comp_flux(flux_z, flux_t, flux_zt, flux_n, land, fname, fw)
    % take zonal average
    flux_z.(land).(fname).(fw) = squeeze(nanmean(flux_n.(land).(fname).(fw), 1));

    % take time averages
    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        if strcmp(time, 'ann')
            flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw), 3);
        elseif strcmp(time, 'djf')
            flux_shift.(land).(fname).(fw) = circshift(flux_n.(land).(fname).(fw), 1, 3);
            flux_t.(land).(time).(fname).(fw) = nanmean(flux_shift.(land).(fname).(fw)(:,:,1:3), 3);
        elseif strcmp(time, 'jja')
            flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,6:8), 3);
        elseif strcmp(time, 'mam')
            flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,3:5), 3);
        elseif strcmp(time, 'son')
            flux_t.(land).(time).(fname).(fw) = nanmean(flux_n.(land).(fname).(fw)(:,:,9:11), 3);
        end
        flux_zt.(land).(time).(fname).(fw) = squeeze(nanmean(flux_t.(land).(time).(fname).(fw), 1));
    end
end
