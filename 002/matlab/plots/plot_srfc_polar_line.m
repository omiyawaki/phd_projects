function plot_dlh_polar_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/srfc.mat', prefix)); % read surface data
    load(sprintf('%s/rad.mat', prefix)); % read surface data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % read fluxes
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [-85 -80 -70 70 80 85];
    lat_bound_list = [-80 80];

    if any(strcmp(type, {'gcm', 'jra55'}))
        rlus = squeeze(nanmean(rad.rlus));
        rlds = squeeze(nanmean(rad.rlds));
    elseif strcmp(type, 'echam')
        rlus = squeeze(nanmean(-rad.tradsu));
        rlds = squeeze(nanmean(-rad.tradsu+rad.trads));
    end

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

        [lh, sh] = rename_stf(type, flux_z, land);

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

                res = lh + sh + flux_z.(land).swsfc + flux_z.(land).lwsfc;
                res_lat = interp1(grid.dim2.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                lh_lat = interp1(grid.dim2.lat, lh, lat);
                lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                sh_lat = interp1(grid.dim2.lat, sh, lat);
                sh_lat = nansum(sh_lat.*clat_mon)/nansum(clat);
                swsfc_lat  = interp1(grid.dim2.lat, flux_z.(land).swsfc, lat);
                swsfc_lat = nansum(swsfc_lat.*clat_mon)/nansum(clat);
                lwsfc_lat  = interp1(grid.dim2.lat, flux_z.(land).lwsfc, lat);
                lwsfc_lat = nansum(lwsfc_lat.*clat_mon)/nansum(clat);
                rlus_lat  = interp1(grid.dim2.lat, rlus, lat);
                rlus_lat = nansum(rlus_lat.*clat_mon)/nansum(clat);
                rlds_lat  = interp1(grid.dim2.lat, rlds, lat);
                rlds_lat = nansum(rlds_lat.*clat_mon)/nansum(clat);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                sw=plot([1:12], circshift(swsfc_lat,shiftby,2), 'color', par.yellow);
                lw=plot([1:12], circshift(lwsfc_lat,shiftby,2), 'color', par.green);
                resf=plot([1:12], circshift(res_lat,shiftby,2), 'k');
                make_title_type_lat(type, lat_bound, lat_pole, par);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([sw, lw, lhf, shf, resf], '$\mathrm{SW_{SFC}}$', '$\mathrm{LW_{SFC}}$', 'LH', 'SH', 'Resiudal', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-80 80], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_srfc', folder), '-dpng', '-r300');
                close;

                if any(strcmp(type, {'gcm', 'jra55', 'echam'}))
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    lw=plot([1:12], circshift(lwsfc_lat,shiftby,2), 'color', par.green);
                    lwup=plot([1:12], circshift(rlus_lat,shiftby,2), '--', 'color', par.green);
                    lwdn=plot([1:12], circshift(rlds_lat,shiftby,2), ':', 'color', par.green);
                    make_title_type_lat(type, lat_bound, lat_pole, par);
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([lw, lwup, lwdn], '$\mathrm{LW_{SFC}}$', '$\mathrm{LW_{SFC, UP}}$', '$\mathrm{LW_{SFC, DN}}$', 'location', 'eastoutside');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                    set(gca, 'ylim', [0 400], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_lwsfc', folder), '-dpng', '-r300');
                    close;
                end

            end

        end % for mse dse
    end % for land

end % for function
