function plot_dlh_polar_asymline_basic(type, par)
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
        esatw = (1.0007 + 3.46e-6*ps/100).*6.1121.*exp(17.966*(tas-273.15)./(247.15+(tas-273.15)));

        % compute saturation vapor pressure w.r.t. ice
        esati = (1.0003 + 4.18e-6*ps/100).*6.1115.*exp(22.452*(tas-273.15)./(272.55+(tas-273.15)));

        % compute RH w.r.t. ice
        if any(strcmp(type, {'era5', 'erai', 'era5c'}))
            e = calc_esat(srfc.d2m, par.frz);
            u_eff = 1e-3; % effective velocity
            ylim_lo = -1;
            ylim_up = 15;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif any(strcmp(type, {'rea', 'gcm'}))
            esat = calc_esat(tas, par.frz);
            e = esat .* srfc.hurs;
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -1;
            ylim_up = 15;
            ts_ylim_lo = 210;
            ts_ylim_up = 280;
        elseif strcmp(type, 'hahn')
            e = calc_e(srfc.PS, srfc.huss, par);
            u_eff = 0.8e-3; % effective velocity
            ylim_lo = -5;
            ylim_up = 25;
            ts_ylim_lo = 210;
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
        lh_est = par.Ls * u_eff*100*esati./(par.Rv*tas);
        lh_estw = par.L * u_eff*100*esatw./(par.Rv*tas);

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

                    if lat_bound>0
                        r1_nh_lat = interp1(grid.dim2.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                        r1_nh_lat = nansum(r1_nh_lat.*clat_mon)/nansum(clat);
                        ra_nh_lat = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat);
                        ra_nh_lat = nansum(ra_nh_lat.*clat_mon)/nansum(clat);
                        lh_nh_lat = interp1(grid.dim3.lat, lh, lat);
                        lh_nh_lat = nansum(lh_nh_lat.*clat_mon)/nansum(clat);
                        tas_nh_lat = interp1(grid.dim2.lat, tas, lat);
                        tas_nh_lat = nansum(tas_nh_lat.*clat_mon)/nansum(clat);
                        ts_nh_lat = interp1(grid.dim2.lat, ts, lat);
                        ts_nh_lat = nansum(ts_nh_lat.*clat_mon)/nansum(clat);
                        esati_nh_lat = interp1(grid.dim2.lat, esati, lat);
                        esati_nh_lat = nansum(esati_nh_lat.*clat_mon)/nansum(clat);
                        lh_est_nh_lat = interp1(grid.dim2.lat, lh_est, lat);
                        lh_est_nh_lat = nansum(lh_est_nh_lat.*clat_mon)/nansum(clat);
                        lh_estw_nh_lat = interp1(grid.dim2.lat, lh_estw, lat);
                        lh_estw_nh_lat = nansum(lh_estw_nh_lat.*clat_mon)/nansum(clat);
                    else
                        r1_sh_lat = interp1(grid.dim2.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                        r1_sh_lat = nansum(r1_sh_lat.*clat_mon)/nansum(clat);
                        ra_sh_lat = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat);
                        ra_sh_lat = nansum(ra_sh_lat.*clat_mon)/nansum(clat);
                        lh_sh_lat = interp1(grid.dim3.lat, lh, lat);
                        lh_sh_lat = nansum(lh_sh_lat.*clat_mon)/nansum(clat);
                        tas_sh_lat = interp1(grid.dim2.lat, tas, lat);
                        tas_sh_lat = nansum(tas_sh_lat.*clat_mon)/nansum(clat);
                        ts_sh_lat = interp1(grid.dim2.lat, ts, lat);
                        ts_sh_lat = nansum(ts_sh_lat.*clat_mon)/nansum(clat);
                        esati_sh_lat = interp1(grid.dim2.lat, esati, lat);
                        esati_sh_lat = nansum(esati_sh_lat.*clat_mon)/nansum(clat);
                        lh_est_sh_lat = interp1(grid.dim2.lat, lh_est, lat);
                        lh_est_sh_lat = nansum(lh_est_sh_lat.*clat_mon)/nansum(clat);
                        lh_estw_sh_lat = interp1(grid.dim2.lat, lh_estw, lat);
                        lh_estw_sh_lat = nansum(lh_estw_sh_lat.*clat_mon)/nansum(clat);
                    end
                end

                    var_text = '$\Delta R_1$';
                    figure(); clf; box on; hold all;
                    ylim_lo = 0.7; 
                    ylim_up = 1.7;
                    rcemax = par.ep;
                    if rcemax > ylim_lo
                        vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    end
                    raemin = par.ga;
                    if raemin < ylim_up
                        vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    end
                    % line([1 12], r1_ann_var*[1 1], 'linewidth', 0.5, 'color', 'k');
                    plot([1:12], circshift(nanmean(r1_nh_lat)*ones(1,12),0,2), 'k', 'linewidth', 0.5);
                    plot([1:12], circshift(nanmean(r1_sh_lat)*ones(1,12),0,2), ':k', 'linewidth', 0.5);
                    pnh=plot([1:12], circshift(r1_nh_lat,0,2), 'k');
                    psh=plot([1:12], circshift(r1_sh_lat,6,2), ':k');
                    ylabel(sprintf('$R_1$ (unitless)'));
                    % xlabel('Month (relative to winter solstice)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_r1_asym', folder), '-dpng', '-r300');
                    if par.make_tikz
                        matlab2tikz(sprintf('%s/0_mon_r1_asym.tex', folder));
                    end
                    close;

                    var_text = '$\Delta R_a$';
                    figure(); clf; box on; hold all;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot([1:12], circshift(nanmean(ra_nh_lat)*ones(1,12),0,2), 'color', par.gray, 'linewidth', 0.5);
                    plot([1:12], circshift(nanmean(ra_sh_lat)*ones(1,12),0,2), ':', 'color', par.gray, 'linewidth', 0.5);
                    pnh=plot([1:12], circshift(ra_nh_lat,0,2), 'color', par.gray);
                    psh=plot([1:12], circshift(ra_sh_lat,6,2), ':', 'color', par.gray);
                    ylabel(sprintf('$R_a$ (W m$^{-2}$)'));
                    xlabel('Month (relative to winter solstice)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-150 0], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_ra_asym', folder), '-dpng', '-r300');
                    close;

                    var_text = '$\mathrm{LH}$';
                    figure(); clf; box on; hold all;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % line([1 12], 5*[1 1], 'linewidth', 0.5, 'color', par.blue);
                    % line([1 12], r1_ann_var*[1 1], 'linewidth', 0.5, 'color', 'k');
                    % plot([1:12], circshift(nanmean(lh_nh_lat)*ones(1,12),0,2), 'color', par.blue, 'linewidth', 0.5);
                    % plot([1:12], circshift(nanmean(lh_sh_lat)*ones(1,12),0,2), ':', 'color', par.blue, 'linewidth', 0.5);
                    pnh=plot([1:12], circshift(lh_nh_lat,0,2), 'color', par.blue);
                    psh=plot([1:12], circshift(lh_sh_lat,6,2), ':', 'color', par.blue);
                    ylabel(sprintf('$\\mathrm{LH}$ (W m$^{-2}$)'));
                    % xlabel('Month (relative to winter solstice)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-5 25], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_asym', folder), '-dpng', '-r300');
                    if par.make_tikz
                        matlab2tikz(sprintf('%s/0_mon_lh_asym.tex', folder));
                    end
                    close;

                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    pts_nh=plot([1:12], circshift(ts_nh_lat,0,2), '-k', 'linewidth', 1.2);
                    pts_sh=plot([1:12], circshift(ts_sh_lat,6,2), ':k', 'linewidth', 1.2);
                    axis off;
                    axis([100 101 100 101])
                    legend([pts_nh pts_sh], 'Northern Hemisphere', 'Southern Hemisphere', 'location', 'northwest', 'orientation', 'horizontal');
                    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 4.5 0.5])
                    print(sprintf('%s/0_asym_leg', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    pts_nh=plot([1:12], circshift(ts_nh_lat,0,2), '-k', 'linewidth', 1.2);
                    pts_sh=plot([1:12], circshift(ts_sh_lat,6,2), ':k', 'linewidth', 1.2);
                    ylabel('$T_{s}$ (K)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    xlabel('Month (relative to winter solstice)');
                    % legend([pts_nh, pts_sh], 'NH', 'SH');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ts_ylim_lo ts_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_ts_asym', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % line([1 12], 260*[1 1], 'linewidth', 0.5, 'color', 'k');
                    ptas_nh=plot([1:12], circshift(tas_nh_lat,0,2), '-k', 'linewidth', 1.2);
                    ptas_sh=plot([1:12], circshift(tas_sh_lat,6,2), ':k', 'linewidth', 1.2);
                    ylabel('$T_{2\,\mathrm{m}}$ (K)');
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    xlabel('Month (relative to winter solstice)');
                    % legend([ptas_nh, ptas_sh], 'NH', 'SH');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'ylim', [ts_ylim_lo ts_ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_tas_asym', folder), '-dpng', '-r300');
                    close;

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

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % lhe=plot([1:12], circshift(lh_est_lat,shiftby,2), '--', 'color', par.blue, 'linewidth', 1.2);
                    % lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    % set(gca, 'ylim', [ylim_lo ylim_up], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % print(sprintf('%s/0_mon_lh_est_lh', folder), '-dpng', '-r300');
                    % close;

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

            end % for mse dse
        end % for land

end % for function
