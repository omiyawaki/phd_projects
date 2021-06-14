function plot_dlh_polar_line_asymdecomp(type, par)
        make_dirs(type, par)

        prefix = make_prefix(type, par);
        prefix_proc = make_prefix_proc(type, par);
        plotdir = make_plotdir(type, par);

        load(sprintf('%s/grid.mat', prefix)); % read grid data
        % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
        load(sprintf('%s/srfc.mat', prefix)); % read melting of ice
        load(sprintf('%s/sfcWind.mat', prefix)); % read melting of ice
        load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
        % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

        % sftlf = nanmean(sftlf, 1); % zonal average
        % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

        tas = rename_tas(type, srfc);
        ts = rename_ts(type, srfc);
        ps = rename_ps(type, srfc);

        % compute saturation vapor pressure w.r.t. water
        esatw = 100 * (1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));

        % compute saturation vapor pressure w.r.t. ice
        esati = 100 * (1.0003 + 4.18e-6*ps/100).*6.1115.*exp(22.452*(tas-273.15)./(272.55+(tas-273.15)));

        sw_ylim_lo = 0;
        sw_ylim_up = 20;
        rh_ylim_lo = 50;
        rh_ylim_up = 120;
        dlh_ylim_lo = -10;
        dlh_ylim_up = 10;

        % compute RH w.r.t. ice
        if any(strcmp(type, {'era5', 'erai', 'era5c'}))
            % e = calc_esat(srfc.d2m, par.frz);
            e = 100 * (1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(srfc.d2m-273.15)./(247.15+(srfc.d2m-273.15)));
            C_e = 2e-3; % water transfer coefficient
            u_eff = 1e-3; % effective velocity
            ylim_lo = -1;
            ylim_up = 15;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'rea', 'gcm'}))
            e = esatw .* srfc.hurs/100;
            C_e = 2e-3; % water transfer coefficient
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -1;
            ylim_up = 15;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'merra2', 'merra2c'}))
            e = calc_e(srfc.PS, srfc.QV2M, par);
            C_e = 2e-3; % water transfer coefficient
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 230;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'jra55'}))
            e = esatw .* srfc.hurs/100;
            C_e = 2e-3; % water transfer coefficient
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 230;
            ts_ylim_up = 280;
        elseif strcmp(type, 'hahn')
            e = calc_e(srfc.PS, srfc.huss, par);
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 230;
            ts_ylim_up = 280;
        elseif strcmp(type, 'echam')
            % e = calc_esat(srfc.dew2, par.frz);
            e = 100 * (1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(srfc.dew2-273.15)./(247.15+(srfc.dew2-273.15)));
            u_eff = 1.4e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 230;
            ts_ylim_up = 280;
        end
        rhw = e./esatw;
        rhi = e./esati;

        % scale estimate of LH (following Andreas 2002)
        u10 = 3; % typical surface wind speed m/s
        ce10 = 2e-3; % bulk transfer coefficient unitless
        rh10 = 0.9; % relative humidity where evaporation takes place
        % u_eff = u10*ce10*(1-rh10); % effective velocity
        lh_est = par.Ls * sfcWind.*C_e.*(1-rhi).*esati./(par.Rv*tas);
        lh_estw = par.L * sfcWind.*C_e.*(1-rhw).*esatw./(par.Rv*tas);

        esati_tas = squeeze(nanmean(esati./tas, 1));
        rhi_def = squeeze(nanmean(1-rhi, 1));

        tas = squeeze(nanmean(tas, 1));
        ts = squeeze(nanmean(ts, 1));
        e = squeeze(nanmean(e));
        esatw = squeeze(nanmean(esatw, 1));
        esati = squeeze(nanmean(esati, 1));
        rhw = squeeze(nanmean(rhw, 1));
        rhi = squeeze(nanmean(rhi, 1));
        lh_est = squeeze(nanmean(lh_est, 1));
        lh_estw = squeeze(nanmean(lh_estw, 1));
        sfcWind = squeeze(nanmean(sfcWind, 1));

        lat_bound_list = [-80 80];

        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam');
                if strcmp(par.echam.clim, '20170908')
                    land_text = 'Snowball';
                else
                    land_text = par.echam.(par.echam.clim);
                end
            end;
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

            f_vec = assign_fw(type, par);
            for f = f_vec; fw = f{1};
                for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                    dlat = 0.25; % step size for standard lat grid
                    if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby=0; monlabel=par.monlabelnh;
                    else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby=6; monlabel=par.monlabelsh; end;
                    clat = cosd(lat); % cosine of latitude for cosine weighting
                    clat_mon = repmat(clat', [1 12]);

                    folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g_asym', plotdir, fw, land, lat_bound);
                    if ~exist(folder, 'dir'); mkdir(folder); end;

                    [lh, sh] = rename_stf(type, flux_z, land);
                    lh_lat = interp1(grid.dim3.lat, lh, lat);
                    lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                    tas_lat = interp1(grid.dim2.lat, tas, lat);
                    tas_lat = nansum(tas_lat.*clat_mon)/nansum(clat);
                    ts_lat = interp1(grid.dim2.lat, ts, lat);
                    ts_lat = nansum(ts_lat.*clat_mon)/nansum(clat);
                    rhw_lat = interp1(grid.dim2.lat, rhw, lat);
                    rhw_lat = nansum(rhw_lat.*clat_mon)/nansum(clat);
                    rhi_lat = interp1(grid.dim2.lat, rhi, lat);
                    rhi_lat = nansum(rhi_lat.*clat_mon)/nansum(clat);
                    rhi_def_lat = interp1(grid.dim2.lat, rhi_def, lat);
                    rhi_def_lat = nansum(rhi_def_lat.*clat_mon)/nansum(clat);
                    esati_lat = interp1(grid.dim2.lat, esati, lat);
                    esati_lat = nansum(esati_lat.*clat_mon)/nansum(clat);
                    esati_tas_lat = interp1(grid.dim2.lat, esati_tas, lat);
                    esati_tas_lat = nansum(esati_tas_lat.*clat_mon)/nansum(clat);
                    lh_est_lat = interp1(grid.dim2.lat, lh_est, lat);
                    lh_est_lat = nansum(lh_est_lat.*clat_mon)/nansum(clat);
                    lh_estw_lat = interp1(grid.dim2.lat, lh_estw, lat);
                    lh_estw_lat = nansum(lh_estw_lat.*clat_mon)/nansum(clat);
                    sfcWind_lat = interp1(grid.dim2.lat, sfcWind, lat);
                    sfcWind_lat = nansum(sfcWind_lat.*clat_mon)/nansum(clat);

                    % summer months
                    if lat_bound > 0
                        % idx_nh = [5,6,7];
                        idx_nh = [6];

                        % summer (DJF for SH, JJA for NH) means and deviations
                        sfcWind_lat_nh = nanmean(sfcWind_lat(idx_nh));
                        esati_tas_lat_nh = nanmean(esati_tas_lat(idx_nh));
                        rhi_def_lat_nh = nanmean(rhi_def_lat(idx_nh));
                        lh_est_lat_nh = nanmean(lh_est_lat(idx_nh));
                        lh_lat_nh = nanmean(lh_lat(idx_nh));
                    else
                        % idx_sh = [11,12,1];
                        idx_sh = [12];

                        % summer (DJF for SH, JJA for NH) means and deviations
                        sfcWind_lat_sh = nanmean(sfcWind_lat(idx_sh));
                        esati_tas_lat_sh = nanmean(esati_tas_lat(idx_sh));
                        rhi_def_lat_sh = nanmean(rhi_def_lat(idx_sh));
                        lh_est_lat_sh = nanmean(lh_est_lat(idx_sh));
                        lh_lat_sh = nanmean(lh_lat(idx_sh));
                    end

                end

                sfcWind_lat_dev = sfcWind_lat_nh - sfcWind_lat_sh;
                esati_tas_lat_dev = esati_tas_lat_nh - esati_tas_lat_sh;
                rhi_def_lat_dev = rhi_def_lat_nh - rhi_def_lat_sh;
                lh_est_lat_dev = lh_est_lat_nh - lh_est_lat_sh;
                lh_lat_dev = lh_lat_nh - lh_lat_sh;

                lh_comp_sfcWind = par.Ls * sfcWind_lat_dev * C_e * rhi_def_lat_sh * esati_tas_lat_sh / par.Rv;
                lh_comp_rhi_def = par.Ls * sfcWind_lat_sh * C_e * rhi_def_lat_dev * esati_tas_lat_sh / par.Rv;
                lh_comp_esati_tas = par.Ls * sfcWind_lat_sh * C_e * rhi_def_lat_sh * esati_tas_lat_dev / par.Rv;

                lh_res = lh_est_lat_dev - (lh_comp_sfcWind + lh_comp_rhi_def + lh_comp_esati_tas);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                bar([lh_lat_dev lh_est_lat_dev lh_comp_esati_tas lh_comp_rhi_def lh_comp_sfcWind lh_res]);
                % make_title_type_lat(type, lat_bound, lat_pole, par);
                % xlabel('Month');
                ylabel(sprintf('$\\Delta\\mathrm{LH}$ (Wm$^{-2}$)'));
                % legend([dlhe, comp_ts, comp_sw, comp_rh, res], '$\Delta \mathrm{LH_{pred}}$', '$L_sC R_v^{-1}\overline{U} \Delta(e_{\mathrm{sat,\,i}}T^{-1})\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \Delta U \overline{e_{\mathrm{sat,\,i}}T^{-1}}\,\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \overline{U} \overline{e_{\mathrm{sat,\,i}}T^{-1}}\Delta(1-\mathrm{RH})$', 'Residual', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq_large)
                set(gca, 'xtick', 1:6, 'xticklabel', {'$\Delta \mathrm{LH}$', '$\Delta \mathrm{LH_{pred}}$', '$L_sC R_v^{-1}\overline{U} \Delta(e_{\mathrm{sat,\,i}}T^{-1})\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \overline{U} \overline{e_{\mathrm{sat,\,i}}T^{-1}}\Delta(1-\mathrm{RH})$', '$L_sCR_v^{-1} \Delta U \overline{e_{\mathrm{sat,\,i}}T^{-1}}\,\overline{(1-\mathrm{RH})}$', 'Residual'}, 'yminortick', 'on', 'tickdir', 'out');
                xtickangle(45);
                print(sprintf('%s/0_mon_dlh_decomp', folder), '-dpng', '-r300');
                close;

            end % for mse dse
        end % for land

end % for function
