function plot_dlh_polar_line(type, par)
    if strcmp(type, 'echam')
        make_dirs(type, par)

        prefix = make_prefix(type, par);
        prefix_proc = make_prefix_proc(type, par);
        plotdir = make_plotdir(type, par);

        load(sprintf('%s/grid.mat', prefix)); % read grid data
        % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
        load(sprintf('%s/ahfliac.mat', prefix)); % read LH over ice
        load(sprintf('%s/ahfllac.mat', prefix)); % read LH over land
        load(sprintf('%s/ahflwac.mat', prefix)); % read LH over water
        load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
        % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

        % sftlf = nanmean(sftlf, 1); % zonal average
        % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

        % take zonal mean of LH components
        ahfliac = -squeeze(nanmean(ahfliac, 1));
        ahfllac = -squeeze(nanmean(ahfllac, 1));
        ahflwac = -squeeze(nanmean(ahflwac, 1));

        % lat_bound_list = [-85 -80 -70 70 80 85];
        lat_bound_list = [-80 80];

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
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
                    if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby=0; monlabel=par.monlabel;
                    else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby=6; monlabel=par.monlabelsh; end;
                    clat = cosd(lat); % cosine of latitude for cosine weighting
                    clat_mon = repmat(clat', [1 12]);

                    folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g', plotdir, fw, land, lat_bound);
                    if ~exist(folder, 'dir'); mkdir(folder); end;

                    [lh, sh] = rename_stf(type, flux_z, land);
                    lh_lat = interp1(grid.dim3.lat, lh, lat);
                    lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                    ahfliac_lat = interp1(grid.dim2.lat, ahfliac, lat);
                    ahfliac_lat = nansum(ahfliac_lat.*clat_mon)/nansum(clat);
                    ahfllac_lat = interp1(grid.dim2.lat, ahfllac, lat);
                    ahfllac_lat = nansum(ahfllac_lat.*clat_mon)/nansum(clat);
                    ahflwac_lat = interp1(grid.dim2.lat, ahflwac, lat);
                    ahflwac_lat = nansum(ahflwac_lat.*clat_mon)/nansum(clat);

                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '--c', 'linewidth', 1.2);
                    lhfl=plot([1:12], circshift(ahfllac_lat,shiftby,2), ':c', 'linewidth', 1.2);
                    lhfw=plot([1:12], circshift(ahflwac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([lhf, lhfi, lhfl, lhfw], '$\mathrm{LH_{tot}}$', '$\mathrm{LH_{ice}}$', '$\mathrm{LH_{land}}$', '$\mathrm{LH_{water}}$',  'location', 'eastoutside');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                    set(gca, 'ylim', [-1 12], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh', folder), '-dpng', '-r300');
                    close;

                    % NOLEG ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                    lhfi=plot([1:12], circshift(ahfliac_lat,shiftby,2), '--c', 'linewidth', 1.2);
                    lhfl=plot([1:12], circshift(ahfllac_lat,shiftby,2), ':c', 'linewidth', 1.2);
                    lhfw=plot([1:12], circshift(ahflwac_lat,shiftby,2), '-.c', 'linewidth', 1.2);
                    make_title_type_lat(type, lat_bound, lat_pole);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'ylim', [-1 12], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lh_noleg', folder), '-dpng', '-r300');
                    close;


                end

            end % for mse dse
        end % for land

    else
        error('This function only works for ECHAM.');
    end

end % for function
