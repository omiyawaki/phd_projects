function plot_dlh_polar_line_basic(type, par)
        make_dirs(type, par)

        prefix = make_prefix(type, par);
        prefix_proc = make_prefix_proc(type, par);
        plotdir = make_plotdir(type, par);

        load(sprintf('%s/grid.mat', prefix)); % read grid data
        % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
        load(sprintf('%s/srfc.mat', prefix)); % read melting of ice
        load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
        % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

        % sftlf = nanmean(sftlf, 1); % zonal average
        % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

        % take zonal mean of LH components
        tas = rename_tas(type, srfc);
        ts = rename_ts(type, srfc);
        ps = rename_ps(type, srfc);

        % compute saturation vapor pressure w.r.t. water
        % esatw = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));
        esatw = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(ts-273.15)./(247.15+(ts-273.15)));

        ps_line = 1e5;
        tas_line = 230:270;
        esatw_line = 1e2*(1.0007 + 3.46e-6*ps_line/100).*6.1121.*exp(17.966*(tas_line-273.15)./(247.15+(tas_line-273.15)));

        % compute saturation vapor pressure w.r.t. ice
        % esati = 1e2*(1.0003 + 4.18e-6*ps/100).*6.1115.*exp(22.452*(tas-273.15)./(272.55+(tas-273.15)));
        esati = 1e2*(1.0003 + 4.18e-6*ps/100).*6.1115.*exp(22.452*(ts-273.15)./(272.55+(ts-273.15)));
        esati_line = 1e2*(1.0003 + 4.18e-6*ps_line/100).*6.1115.*exp(22.452*(tas_line-273.15)./(272.55+(tas_line-273.15)));

        % ALL lat x mon dependence of RCE and RAE
        figure(); clf; hold all; box on;
        % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
        pew=plot(tas_line, esatw_line, '-k', 'linewidth', 1.2);
        pei=plot(tas_line, esati_line, '--k', 'linewidth', 1.2);
        xlabel('T (K)');
        ylabel('$e_{\mathrm{sat}}$ (Pa)');
        % make_title_type_lat(type, lat_bound, lat_pole, par);
        legend('water', 'ice', 'location', 'northwest');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'xminortick', 'on', 'yminortick', 'on');
        folder=sprintf('/project2/tas1/miyawaki/projects/002/figures/general');
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
        print(sprintf('%s/0_t_esat.png', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
        pew=plot(tas_line, esatw_line./esati_line, '-k', 'linewidth', 1.2);
        xlabel('T (K)');
        ylabel('$e_{\mathrm{sat, \,w}} / e_{\mathrm{sat,\,i}}$ (unitless)');
        % make_title_type_lat(type, lat_bound, lat_pole, par);
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'xminortick', 'on', 'yminortick', 'on');
        folder=sprintf('/project2/tas1/miyawaki/projects/002/figures/general');
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
        print(sprintf('%s/0_t_esat_ratio.png', folder), '-dpng', '-r300');
        close;

        rh_ylim_lo = 50;
        rh_ylim_up = 120;

        % compute RH w.r.t. ice
        if any(strcmp(type, {'era5', 'erai', 'era5c'}))
            % e = calc_esat(srfc.d2m, par.frz);
            e = 1e2*(1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(srfc.d2m-273.15)./(247.15+(srfc.d2m-273.15)));
            u_eff = 1e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'rea', 'gcm'}))
            e = esatw .* srfc.hurs/100;
            u_eff = 1e-3; % effective velocity
            ylim_lo = -1;
            ylim_up = 15;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'merra2', 'merra2c'}))
            e = calc_e(srfc.PS, srfc.QV2M, par);
            u_eff = 1e-3; % effective velocity
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
            u_eff = 1e-3; % effective velocity
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
            e = calc_esat(srfc.dew2, par.frz);
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
        lh_est = par.Ls * u_eff*esati./(par.Rv*tas);
        lh_estw = par.L * u_eff*esatw./(par.Rv*tas);

        tas = squeeze(nanmean(tas, 1));
        ts = squeeze(nanmean(ts, 1));
        e = squeeze(nanmean(e));
        esatw = squeeze(nanmean(esatw, 1));
        esati = squeeze(nanmean(esati, 1));
        rhw = squeeze(nanmean(rhw, 1));
        rhi = squeeze(nanmean(rhi, 1));
        lh_est = squeeze(nanmean(lh_est, 1));
        lh_estw = squeeze(nanmean(lh_estw, 1));

        % lat_bound_list = [-85 -80 -70 70 80 85];
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
                    rhi_lat = interp1(grid.dim2.lat, rhi, lat);
                    rhi_lat = nansum(rhi_lat.*clat_mon)/nansum(clat);
                    esati_lat = interp1(grid.dim2.lat, esati, lat);
                    esati_lat = nansum(esati_lat.*clat_mon)/nansum(clat);
                    lh_est_lat = interp1(grid.dim2.lat, lh_est, lat);
                    lh_est_lat = nansum(lh_est_lat.*clat_mon)/nansum(clat);
                    lh_estw_lat = interp1(grid.dim2.lat, lh_estw, lat);
                    lh_estw_lat = nansum(lh_estw_lat.*clat_mon)/nansum(clat);

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    pts=plot([1:12], circshift(ts_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{s}$ (K)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ts_ylim_lo ts_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_ts', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    ptas=plot([1:12], circshift(tas_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{2\,\mathrm{m}}$ (K)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ts_ylim_lo ts_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_tas', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    prhw=plot([1:12], 100 * circshift(rhw_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$\mathrm{RH_{w,\,2\,m}}$ (\%)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [rh_ylim_lo rh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_rhw', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    prhi=plot([1:12], 100 * circshift(rhi_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$\mathrm{RH_{i,\,2\,m}}$ (\%)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [rh_ylim_lo rh_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_rhi', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    ptas=plot([1:12], circshift(tas_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{2\,\mathrm{m}}$ (K)');
                    yyaxis left
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = par.blue;
                    print(sprintf('%s/0_mon_tas_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    pts=plot([1:12], circshift(ts_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{s}$ (K)');
                    yyaxis left
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = par.blue;
                    print(sprintf('%s/0_mon_ts_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    ei=plot([1:12], circshift(esati_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$e_{\mathrm{esat,\, i}}$ (hPa)');
                    yyaxis left
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = par.blue;
                    print(sprintf('%s/0_mon_esati_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh', folder), '-dpng', '-r300');
                    if par.make_tikz
                        matlab2tikz(sprintf('%s/0_mon_lh.tex', folder));
                    end
                    close;

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
                    print(sprintf('%s/0_mon_lh_est_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhe=plot([1:12], circshift(lh_estw_lat,shiftby,2), '-.', 'color', par.blue, 'linewidth', 1.2);
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_estw_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    axis off;
                    axis([10,11,10,11])
                    legend([lhf lhe], 'LH', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, i}}T^{-1}$', 'location', 'northwest', 'orientation', 'horizontal');
                    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 3.3 0.5])
                    print(sprintf('%s/0_mon_lh_est_lh_legonly', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    lhew=plot([1:12], circshift(lh_estw_lat,shiftby,2), '-.', 'color', par.blue, 'linewidth', 1.2);
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    axis off;
                    axis([10,11,10,11])
                    legend([lhf lhe lhew], 'LH', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, i}}T^{-1}$', '$\mathrm{LH_{pred}}\propto e_{\mathrm{esat,\, w}}T^{-1}$', 'location', 'northwest', 'orientation', 'horizontal');
                    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 5.5 0.5])
                    print(sprintf('%s/0_mon_lh_estw_lh_legonly', folder), '-dpng', '-r300');
                    close;

                end

            end % for mse dse
        end % for land

end % for function
