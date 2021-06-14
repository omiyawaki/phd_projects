function plot_dlh_polar_line_decomp(type, par)
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
            load(sprintf('%s/hus_ml0.mat', prefix)); % read melting of ice
            hus_ml0 = interp1(hus_ml0_lon, hus_ml0, grid.dim3.lon);
            hus_ml0 = permute(hus_ml0, [2 1 3]);
            hus_ml0 = interp1(hus_ml0_lat, hus_ml0, grid.dim3.lat);
            hus_ml0 = permute(hus_ml0, [2 1 3]);
            e = calc_e(srfc.ps, hus_ml0, par);
            % e = esatw .* srfc.hurs/100;
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

        esatw_tas = squeeze(nanmean(esatw./tas, 1));
        esati_tas = squeeze(nanmean(esati./tas, 1));
        rhw_def = squeeze(nanmean(1-rhw, 1));
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

                    folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g', plotdir, fw, land, lat_bound);
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
                    rhw_def_lat = interp1(grid.dim2.lat, rhw_def, lat);
                    rhw_def_lat = nansum(rhw_def_lat.*clat_mon)/nansum(clat);
                    rhi_lat = interp1(grid.dim2.lat, rhi, lat);
                    rhi_lat = nansum(rhi_lat.*clat_mon)/nansum(clat);
                    rhi_def_lat = interp1(grid.dim2.lat, rhi_def, lat);
                    rhi_def_lat = nansum(rhi_def_lat.*clat_mon)/nansum(clat);
                    esatw_lat = interp1(grid.dim2.lat, esatw, lat);
                    esatw_lat = nansum(esatw_lat.*clat_mon)/nansum(clat);
                    esatw_tas_lat = interp1(grid.dim2.lat, esatw_tas, lat);
                    esatw_tas_lat = nansum(esatw_tas_lat.*clat_mon)/nansum(clat);
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

                    % annual means and deviations
                    sfcWind_lat_ann = nanmean(sfcWind_lat);
                    esatw_tas_lat_ann = nanmean(esatw_tas_lat);
                    esati_tas_lat_ann = nanmean(esati_tas_lat);
                    rhw_def_lat_ann = nanmean(rhw_def_lat);
                    rhi_def_lat_ann = nanmean(rhi_def_lat);
                    lh_est_lat_ann = nanmean(lh_est_lat);
                    lh_estw_lat_ann = nanmean(lh_estw_lat);

                    sfcWind_lat_dev = sfcWind_lat - sfcWind_lat_ann;
                    esatw_tas_lat_dev = esati_tas_lat - esati_tas_lat_ann;
                    esati_tas_lat_dev = esati_tas_lat - esati_tas_lat_ann;
                    rhw_def_lat_dev = rhw_def_lat - rhw_def_lat_ann;
                    rhi_def_lat_dev = rhi_def_lat - rhi_def_lat_ann;
                    lh_est_lat_dev = lh_est_lat - lh_est_lat_ann;
                    lh_estw_lat_dev = lh_estw_lat - lh_estw_lat_ann;

                    lh_comp_sfcWind = par.Ls * sfcWind_lat_dev * C_e * rhi_def_lat_ann * esati_tas_lat_ann / par.Rv;
                    lh_comp_rhi_def = par.Ls * sfcWind_lat_ann * C_e * rhi_def_lat_dev * esati_tas_lat_ann / par.Rv;
                    lh_comp_esati_tas = par.Ls * sfcWind_lat_ann * C_e * rhi_def_lat_ann * esati_tas_lat_dev / par.Rv;

                    lhw_comp_sfcWind = par.L * sfcWind_lat_dev * C_e * rhw_def_lat_ann * esatw_tas_lat_ann / par.Rv;
                    lhw_comp_rhw_def = par.L * sfcWind_lat_ann * C_e * rhw_def_lat_dev * esatw_tas_lat_ann / par.Rv;
                    lhw_comp_esatw_tas = par.L * sfcWind_lat_ann * C_e * rhw_def_lat_ann * esatw_tas_lat_dev / par.Rv;

                    lh_res = lh_est_lat_dev - (lh_comp_sfcWind + lh_comp_rhi_def + lh_comp_esati_tas);
                    lhw_res = lh_estw_lat_dev - (lhw_comp_sfcWind + lhw_comp_rhw_def + lhw_comp_esatw_tas);

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    psw=plot([1:12], circshift(sfcWind_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$U_{10\,\mathrm{m}}$ (ms$^{-1}$)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [sw_ylim_lo sw_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_sfcWind', folder), '-dpng', '-r300');
                    close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % prhw=plot([1:12], 100 * circshift(rhw_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('$\mathrm{RH_{w,\,2\,m}}$ (\%)');
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'ylim', [rh_ylim_lo rh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % print(sprintf('%s/0_mon_rhw', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % prhi=plot([1:12], 100 * circshift(rhi_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('$\mathrm{RH_{i,\,2\,m}}$ (\%)');
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'ylim', [rh_ylim_lo rh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % print(sprintf('%s/0_mon_rhi', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % yyaxis right
                    % ptas=plot([1:12], circshift(tas_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('$T_{2\,\mathrm{m}}$ (K)');
                    % yyaxis left
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % ax = gca;
                    % ax.YAxis(2).Color = 'k';
                    % ax.YAxis(1).Color = par.blue;
                    % print(sprintf('%s/0_mon_tas_lh', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % yyaxis right
                    % pts=plot([1:12], circshift(ts_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('$T_{s}$ (K)');
                    % yyaxis left
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % ax = gca;
                    % ax.YAxis(2).Color = 'k';
                    % ax.YAxis(1).Color = par.blue;
                    % print(sprintf('%s/0_mon_ts_lh', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % yyaxis right
                    % ei=plot([1:12], circshift(esati_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('$e_{\mathrm{esat,\, i}}$ (hPa)');
                    % yyaxis left
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % ax = gca;
                    % ax.YAxis(2).Color = 'k';
                    % ax.YAxis(1).Color = par.blue;
                    % print(sprintf('%s/0_mon_esati_lh', folder), '-dpng', '-r300');
                    % close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_est_full_lh', folder), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhe=plot([1:12], circshift(lh_estw_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_estw_full_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    dlhe=plot([1:12], circshift(lh_est_lat_dev,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    comp_ts=plot([1:12], circshift(lh_comp_esati_tas,shiftby,2), '-', 'color', 'r', 'linewidth', 1.2);
                    comp_sw=plot([1:12], circshift(lh_comp_sfcWind,shiftby,2), '-', 'color', 'b', 'linewidth', 1.2);
                    comp_rh=plot([1:12], circshift(lh_comp_rhi_def,shiftby,2), '-', 'color', 'g', 'linewidth', 1.2);
                    res=plot([1:12], circshift(lh_res,shiftby,2), '-', 'color', 0.7*[1 1 1], 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('$\\Delta\\mathrm{LH}$ (Wm$^{-2}$)'));
                    legend([dlhe, comp_ts, comp_sw, comp_rh, res], '$\Delta \mathrm{LH_{pred}}$', '$L_sC R_v^{-1}\overline{U} \Delta(e_{\mathrm{sat,\,i}}T^{-1})\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \Delta U \overline{e_{\mathrm{sat,\,i}}T^{-1}}\,\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \overline{U} \overline{e_{\mathrm{sat,\,i}}T^{-1}}\Delta(1-\mathrm{RH})$', 'Residual', 'location', 'eastoutside');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                    set(gca, 'ylim', [dlh_ylim_lo dlh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dlh_decomp', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    dlhe=plot([1:12], circshift(lh_est_lat_dev,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    comp_ts=plot([1:12], circshift(lh_comp_esati_tas,shiftby,2), '-', 'color', 'r', 'linewidth', 1.2);
                    comp_sw=plot([1:12], circshift(lh_comp_sfcWind,shiftby,2), '-', 'color', 'b', 'linewidth', 1.2);
                    comp_rh=plot([1:12], circshift(lh_comp_rhi_def,shiftby,2), '-', 'color', 'g', 'linewidth', 1.2);
                    res=plot([1:12], circshift(lh_res,shiftby,2), '-', 'color', 0.7*[1 1 1], 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('$\\Delta\\mathrm{LH}$ (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [dlh_ylim_lo dlh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dlh_decomp_noleg', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    dlhe=plot([1:12], circshift(lh_est_lat_dev,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    comp_ts=plot([1:12], circshift(lh_comp_esati_tas,shiftby,2), '-', 'color', 'r', 'linewidth', 1.2);
                    comp_sw=plot([1:12], circshift(lh_comp_sfcWind,shiftby,2), '-', 'color', 'b', 'linewidth', 1.2);
                    comp_rh=plot([1:12], circshift(lh_comp_rhi_def,shiftby,2), '-', 'color', 'g', 'linewidth', 1.2);
                    res=plot([1:12], circshift(lh_res,shiftby,2), '-', 'color', 0.7*[1 1 1], 'linewidth', 1.2);
                    % xlabel('Month');
                    ylabel(sprintf('$\\Delta\\mathrm{LH}$ (Wm$^{-2}$)'));
                    axis off
                    axis([10,11,10,11])
                    legend([dlhe, comp_ts, comp_sw, comp_rh, res], '$\Delta \mathrm{LH_{pred}}$', '$L_sC R_v^{-1}\overline{U} \Delta(e_{\mathrm{sat,\,i}}T^{-1})\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \Delta U \overline{e_{\mathrm{sat,\,i}}T^{-1}}\,\overline{(1-\mathrm{RH})}$', '$L_sCR_v^{-1} \overline{U} \overline{e_{\mathrm{sat,\,i}}T^{-1}}\Delta(1-\mathrm{RH})$', 'Residual', 'location', 'northwest');
                    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 3.3 1.2])
                    print(sprintf('%s/0_mon_dlh_decomp_legonly', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    dlhwe=plot([1:12], circshift(lh_estw_lat_dev,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    comp_ts=plot([1:12], circshift(lhw_comp_esatw_tas,shiftby,2), '-', 'color', 'r', 'linewidth', 1.2);
                    comp_sw=plot([1:12], circshift(lhw_comp_sfcWind,shiftby,2), '-', 'color', 'b', 'linewidth', 1.2);
                    comp_rh=plot([1:12], circshift(lhw_comp_rhw_def,shiftby,2), '-', 'color', 'g', 'linewidth', 1.2);
                    res=plot([1:12], circshift(lhw_res,shiftby,2), '-', 'color', 0.7*[1 1 1], 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('$\\Delta\\mathrm{LH}$ (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [dlh_ylim_lo dlh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dlhw_decomp_noleg', folder), '-dpng', '-r300');
                    close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % lhe=plot([1:12], circshift(lh_estw_lat,shiftby,2), '-.', 'color', par.blue, 'linewidth', 1.2);
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % print(sprintf('%s/0_mon_lh_estw_lh', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % axis off;
                    % axis([10,11,10,11])
                    % legend([lhf lhe], 'LH', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, i}}T^{-1}$', 'location', 'northwest', 'orientation', 'horizontal');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 3.3 0.5])
                    % print(sprintf('%s/0_mon_lh_est_lh_legonly', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    % lhew=plot([1:12], circshift(lh_estw_lat,shiftby,2), '-.', 'color', par.blue, 'linewidth', 1.2);
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % axis off;
                    % axis([10,11,10,11])
                    % legend([lhf lhe lhew], 'LH', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, i}}T^{-1}$', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, w}}T^{-1}$', 'location', 'northwest', 'orientation', 'horizontal');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 5.5 0.5])
                    % print(sprintf('%s/0_mon_lh_estw_lh_legonly', folder), '-dpng', '-r300');
                    % close;

                end

            end % for mse dse
        end % for land

end % for function
