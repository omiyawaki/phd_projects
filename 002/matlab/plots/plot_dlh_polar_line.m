function plot_dlh_polar_line(type, par)
    if strcmp(type, 'echam')
        make_dirs(type, par)

        prefix = make_prefix(type, par);
        prefix_proc = make_prefix_proc(type, par);
        plotdir = make_plotdir(type, par);

        load(sprintf('%s/grid.mat', prefix)); % read grid data
        % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
        load(sprintf('%s/rad.mat', prefix)); % read melting of ice
        load(sprintf('%s/srfc.mat', prefix)); % read melting of ice
        load(sprintf('%s/friac.mat', prefix)); friac = echamvar; % read melting of ice
        load(sprintf('%s/ahfres.mat', prefix)); ahfres = echamvar; % read melting of ice
        load(sprintf('%s/ahfliac.mat', prefix)); ahfliac = echamvar; % read LH over ice
        load(sprintf('%s/ahfllac.mat', prefix)); ahfllac = echamvar; % read LH over land
        load(sprintf('%s/ahflwac.mat', prefix)); ahflwac = echamvar; % read LH over water
        % load(sprintf('%s/ameltdepth.mat', prefix)); %ameltdepth = echamvar; % read LH over water
        % load(sprintf('%s/ameltfrac.mat', prefix)); %ameltfrac = echamvar; % read LH over water
        load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
        % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

        % sftlf = nanmean(sftlf, 1); % zonal average
        % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

        % take zonal mean of LH components
        tas = squeeze(nanmean(srfc.temp2, 1));
        ts = squeeze(nanmean(srfc.tsurf, 1));
        swabs = squeeze(nanmean(rad.srad0, 1));
        friac = squeeze(nanmean(friac, 1));
        ahfres = squeeze(nanmean(ahfres, 1));
        ahfliac = -squeeze(nanmean(ahfliac, 1));
        ahfllac = -squeeze(nanmean(ahfllac, 1));
        ahflwac = -squeeze(nanmean(ahflwac, 1));
        % ameltdepth = squeeze(nanmean(ameltdepth, 1));
        % ameltfrac = squeeze(nanmean(ameltfrac, 1));

        % compute saturation vapor pressure w.r.t. water
        esatw = (1.0007 + 3.46e-6*srfc.aps/100).*6.1121.*exp(17.966*(srfc.temp2-273.15)./(247.15+(srfc.temp2-273.15)));

        % compute saturation vapor pressure w.r.t. ice
        esati = (1.0003 + 4.18e-6*srfc.aps/100).*6.1115.*exp(22.452*(srfc.temp2-273.15)./(272.55+(srfc.temp2-273.15)));

        % compute RH w.r.t. ice
        e = calc_esat(srfc.dew2, par.frz);
        rhw = e./esatw;
        rhi = e./esati;

        e = squeeze(nanmean(e));
        esatw = squeeze(nanmean(esatw, 1));
        esati = squeeze(nanmean(esati, 1));
        rhw = squeeze(nanmean(rhw, 1));
        rhi = squeeze(nanmean(rhi, 1));

        rhw(1,:)
        rhi(1,:)
        return

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
                    swabs_lat = interp1(grid.dim2.lat, swabs, lat);
                    swabs_lat = nansum(swabs_lat.*clat_mon)/nansum(clat);
                    tas_lat = interp1(grid.dim2.lat, tas, lat);
                    tas_lat = nansum(tas_lat.*clat_mon)/nansum(clat);
                    ts_lat = interp1(grid.dim2.lat, ts, lat);
                    ts_lat = nansum(ts_lat.*clat_mon)/nansum(clat);
                    friac_lat = interp1(grid.dim2.lat, friac, lat);
                    friac_lat = nansum(friac_lat.*clat_mon)/nansum(clat);
                    ahfres_lat = interp1(grid.dim2.lat, ahfres, lat);
                    ahfres_lat = nansum(ahfres_lat.*clat_mon)/nansum(clat);
                    ahfliac_lat = interp1(grid.dim2.lat, ahfliac, lat);
                    ahfliac_lat = nansum(ahfliac_lat.*clat_mon)/nansum(clat);
                    ahfllac_lat = interp1(grid.dim2.lat, ahfllac, lat);
                    ahfllac_lat = nansum(ahfllac_lat.*clat_mon)/nansum(clat);
                    ahflwac_lat = interp1(grid.dim2.lat, ahflwac, lat);
                    ahflwac_lat = nansum(ahflwac_lat.*clat_mon)/nansum(clat);
                    esati_lat = interp1(grid.dim2.lat, esati, lat);
                    esati_lat = nansum(esati_lat.*clat_mon)/nansum(clat);
                    % ameltdepth_lat = interp1(grid.dim2.lat, ameltdepth, lat);
                    % ameltdepth_lat = nansum(ameltdepth_lat.*clat_mon)/nansum(clat);
                    % ameltfrac_lat = interp1(grid.dim2.lat, ameltfrac, lat);
                    % ameltfrac_lat = nansum(ameltfrac_lat.*clat_mon)/nansum(clat);

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '--c', 'linewidth', 1.2);
                    lhfl=plot([1:12], circshift(ahfllac_lat,shiftby,2), ':c', 'linewidth', 1.2);
                    lhfw=plot([1:12], circshift(ahflwac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([lhf, lhfi, lhfl, lhfw], '$\mathrm{LH_{tot}}$', '$\mathrm{LH_{ice}}$', '$\mathrm{LH_{land}}$', '$\mathrm{LH_{water}}$',  'location', 'eastoutside');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                    set(gca, 'ylim', [-1 20], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh', folder), '-dpng', '-r300');
                    close;

                    % NOLEG ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '--c', 'linewidth', 1.2);
                    lhfl=plot([1:12], circshift(ahfllac_lat,shiftby,2), ':c', 'linewidth', 1.2);
                    lhfw=plot([1:12], circshift(ahflwac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'ylim', [-1 20], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_noleg', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    icef=plot([1:12], circshift(friac_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('Sea ice fraction (unitless)');
                    set(gca, 'ydir', 'reverse');
                    yyaxis left
                    lhfw=plot([1:12], circshift(ahflwac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH over water (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = 'c';
                    print(sprintf('%s/0_mon_icef_lhw', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    ei=plot([1:12], circshift(esati_lat,shiftby,2), '-c', 'linewidth', 1.2);
                    ylabel('$e_{\mathrm{sat,\, i}}$ (hPa)');
                    yyaxis left
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'c';
                    ax.YAxis(1).Color = 'k';
                    print(sprintf('%s/0_mon_esati_lh', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    melt=plot([1:12], circshift(ahfres_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '--c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([melt, lhfi], '$\mathrm{F_{melt}}$', '$\mathrm{LH_{ice}}$', 'location', 'eastoutside');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                    set(gca, 'ylim', [-1 50], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_melting_lhi', folder), '-dpng', '-r300');
                    close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % yyaxis right
                    % ptas=plot([1:12], circshift(ameltfrac_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('Meltpond fraction (unitless)');
                    % yyaxis left
                    % lhfi=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % ax = gca;
                    % ax.YAxis(2).Color = 'k';
                    % ax.YAxis(1).Color = par.blue;
                    % print(sprintf('%s/0_mon_ameltfrac_lhi', folder), '-dpng', '-r300');
                    % close;

                    % % ALL lat x mon dependence of RCE and RAE
                    % figure(); clf; hold all; box on;
                    % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % yyaxis right
                    % ptas=plot([1:12], circshift(ameltdepth_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    % ylabel('Meltpond depth (m)');
                    % yyaxis left
                    % lhfi=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    % make_title_type_lat(type, lat_bound, lat_pole, par);
                    % % xlabel('Month');
                    % ylabel(sprintf('LH (Wm$^{-2}$)'));
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    % ax = gca;
                    % ax.YAxis(2).Color = 'k';
                    % ax.YAxis(1).Color = par.blue;
                    % print(sprintf('%s/0_mon_ameltdepth_lhi', folder), '-dpng', '-r300');
                    % close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    ptas=plot([1:12], circshift(tas_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{2\,\mathrm{m}}$ (K)');
                    yyaxis left
                    lhfi=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = par.blue;
                    print(sprintf('%s/0_mon_tas_lhi', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    pts=plot([1:12], circshift(ts_lat,shiftby,2), '-k', 'linewidth', 1.2);
                    ylabel('$T_{s}$ (K)');
                    yyaxis left
                    lhfi=plot([1:12], circshift(lh_lat,shiftby,2), '-', 'color', par.blue, 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = 'k';
                    ax.YAxis(1).Color = par.blue;
                    print(sprintf('%s/0_mon_ts_lhi', folder), '-dpng', '-r300');
                    close;

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    yyaxis right
                    psw=plot([1:12], circshift(swabs_lat,shiftby,2), 'color', par.yellow, 'linewidth', 1.2);
                    ylabel('$\mathrm{SW_{sfc}}$ (K)');
                    yyaxis left
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('LH over water (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    ax = gca;
                    ax.YAxis(2).Color = par.yellow;
                    ax.YAxis(1).Color = 'c';
                    print(sprintf('%s/0_mon_swabs_lhi', folder), '-dpng', '-r300');
                    close;


                end

            end % for mse dse
        end % for land

    else
        error('This function only works for ECHAM.');
    end

end % for function
