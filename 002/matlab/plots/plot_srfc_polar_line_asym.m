function plot_dlh_polar_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/srfc.mat', prefix)); % read surface data
    load(sprintf('%s/hydro.mat', prefix)); % read evap data
    load(sprintf('%s/snevap.mat', prefix)); % read evap over snow data
    load(sprintf('%s/snmelt.mat', prefix)); % read snowmelt data
    load(sprintf('%s/rad.mat', prefix)); % read surface data
    load(sprintf('%s/radcs.mat', prefix)); % read surface data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % read fluxes
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    rhsfc = comp_rh(type, srfc, par); % compute surface relative humidity
    rhsfc = squeeze(nanmean(rhsfc, 1));

    evap = rename_evap(type, hydro); % use common name for evaporation
    evap = squeeze(nanmean(evap, 1));

    snevap = squeeze(nanmean(snevap, 1));
    snmelt = squeeze(nanmean(snmelt, 1));

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [-85 -80 -70 70 80 85];
    lat_bound_list = 80*[-1 1];

    if any(strcmp(type, {'gcm', 'jra55'}))
        rlus = squeeze(nanmean(rad.rlus,1));
        rlds = squeeze(nanmean(rad.rlds,1));
        rldscs = squeeze(nanmean(radcs.rldscs,1));
        rsdscs = squeeze(nanmean(radcs.rsdscs,1));
        rsuscs = squeeze(nanmean(radcs.rsuscs,1));
        lwsfc_cs = rlus-rldscs;
        swsfc_cs = rsuscs-rsdscs;
    elseif any(strcmp(type, {'era5', 'era5c', 'erai'}))
        lwsfc_cs = squeeze(nanmean(-radcs.strc,1));
        swsfc_cs = squeeze(nanmean(-radcs.ssrc,1));
    elseif strcmp(type, 'merra2')
        lwsfc_cs = squeeze(nanmean(-radcs.LWGNTCLR,1));
        swsfc_cs = squeeze(nanmean(-radcs.SWGNTCLR,1));
    elseif strcmp(type, 'echam')
        rlus = squeeze(nanmean(-rad.tradsu,1));
        rlds = squeeze(nanmean(-rad.tradsu+rad.trads,1));
        lwsfc_cs = squeeze(nanmean(-radcs.trafs,1));
        swsfc_cs = squeeze(nanmean(-radcs.srafs,1));
    end

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

        [lh, sh] = rename_stf(type, flux_z, land);

        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby_nh=0; monlabel=par.monlabelnh;
                else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby_sh=6; monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                res = lh + sh + flux_z.(land).swsfc + flux_z.(land).lwsfc;
                if lat_bound>0;
                    res_lat_nh = interp1(grid.dim2.lat, res, lat);
                    res_lat_nh = nansum(res_lat_nh.*clat_mon)/nansum(clat);
                    lh_lat_nh = interp1(grid.dim2.lat, lh, lat);
                    lh_lat_nh = nansum(lh_lat_nh.*clat_mon)/nansum(clat);
                    sh_lat_nh = interp1(grid.dim2.lat, sh, lat);
                    sh_lat_nh = nansum(sh_lat_nh.*clat_mon)/nansum(clat);
                    swsfc_lat_nh  = interp1(grid.dim2.lat, flux_z.(land).swsfc, lat);
                    swsfc_lat_nh = nansum(swsfc_lat_nh.*clat_mon)/nansum(clat);
                    lwsfc_lat_nh  = interp1(grid.dim2.lat, flux_z.(land).lwsfc, lat);
                    lwsfc_lat_nh = nansum(lwsfc_lat_nh.*clat_mon)/nansum(clat);
                    lwsfc_cs_lat_nh  = interp1(grid.dim2.lat, lwsfc_cs, lat);
                    lwsfc_cs_lat_nh = nansum(lwsfc_cs_lat_nh.*clat_mon)/nansum(clat);
                    swsfc_cs_lat_nh  = interp1(grid.dim2.lat, swsfc_cs, lat);
                    swsfc_cs_lat_nh = nansum(swsfc_cs_lat_nh.*clat_mon)/nansum(clat);
                    rhsfc_lat_nh = interp1(grid.dim2.lat, rhsfc, lat);
                    rhsfc_lat_nh = nansum(rhsfc_lat_nh.*clat_mon)/nansum(clat);
                    evap_lat_nh = interp1(grid.dim2.lat, evap, lat);
                    evap_lat_nh = nansum(evap_lat_nh.*clat_mon)/nansum(clat);
                    snevap_lat_nh = interp1(grid.dim2.lat, snevap, lat);
                    snevap_lat_nh = nansum(snevap_lat_nh.*clat_mon)/nansum(clat);
                    snmelt_lat_nh = interp1(grid.dim2.lat, snmelt, lat);
                    snmelt_lat_nh = nansum(snmelt_lat_nh.*clat_mon)/nansum(clat);
                else
                    res_lat_sh = interp1(grid.dim2.lat, res, lat);
                    res_lat_sh = nansum(res_lat_sh.*clat_mon)/nansum(clat);
                    lh_lat_sh = interp1(grid.dim2.lat, lh, lat);
                    lh_lat_sh = nansum(lh_lat_sh.*clat_mon)/nansum(clat);
                    sh_lat_sh = interp1(grid.dim2.lat, sh, lat);
                    sh_lat_sh = nansum(sh_lat_sh.*clat_mon)/nansum(clat);
                    swsfc_lat_sh  = interp1(grid.dim2.lat, flux_z.(land).swsfc, lat);
                    swsfc_lat_sh = nansum(swsfc_lat_sh.*clat_mon)/nansum(clat);
                    lwsfc_lat_sh  = interp1(grid.dim2.lat, flux_z.(land).lwsfc, lat);
                    lwsfc_lat_sh = nansum(lwsfc_lat_sh.*clat_mon)/nansum(clat);
                    lwsfc_cs_lat_sh  = interp1(grid.dim2.lat, lwsfc_cs, lat);
                    lwsfc_cs_lat_sh = nansum(lwsfc_cs_lat_sh.*clat_mon)/nansum(clat);
                    swsfc_cs_lat_sh  = interp1(grid.dim2.lat, swsfc_cs, lat);
                    swsfc_cs_lat_sh = nansum(swsfc_cs_lat_sh.*clat_mon)/nansum(clat);
                    rhsfc_lat_sh = interp1(grid.dim2.lat, rhsfc, lat);
                    rhsfc_lat_sh = nansum(rhsfc_lat_sh.*clat_mon)/nansum(clat);
                    evap_lat_sh = interp1(grid.dim2.lat, evap, lat);
                    evap_lat_sh = nansum(evap_lat_sh.*clat_mon)/nansum(clat);
                    snevap_lat_sh = interp1(grid.dim2.lat, snevap, lat);
                    snevap_lat_sh = nansum(snevap_lat_sh.*clat_mon)/nansum(clat);
                    snmelt_lat_sh = interp1(grid.dim2.lat, snmelt, lat);
                    snmelt_lat_sh = nansum(snmelt_lat_sh.*clat_mon)/nansum(clat);
                end
            end

            folder = sprintf('%s/dmse/%s/%s/0_polar_asym_lat_%g', plotdir, fw, land, abs(lat_bound_list(1)));
            if ~exist(folder, 'dir'); mkdir(folder); end;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            sw_nh=plot([1:12], -circshift(swsfc_lat_nh,shiftby_nh,2), 'color', par.yellow);
            sw_sh=plot([1:12], -circshift(swsfc_lat_sh,shiftby_sh,2), '--', 'color', par.yellow);
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{SW_{SFC}}$ (Wm$^{-2}$)'));
            legend([sw_nh, sw_sh], 'NH', 'SH', 'location', 'northeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', [0 120], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_sw', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            lw_nh=plot([1:12], -circshift(lwsfc_lat_nh,shiftby_nh,2), 'color', par.green);
            lw_sh=plot([1:12], -circshift(lwsfc_lat_sh,shiftby_sh,2), '--', 'color', par.green);
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{LW_{SFC}}$ (Wm$^{-2}$)'));
            legend([lw_nh, lw_sh], 'NH', 'SH', 'location', 'southeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', [-120 0], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_lw', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            lh_nh=plot([1:12], -circshift(lh_lat_nh,shiftby_nh,2), 'color', par.blue);
            lh_sh=plot([1:12], -circshift(lh_lat_sh,shiftby_sh,2), '--', 'color', par.blue);
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{LH}$ (Wm$^{-2}$)'));
            legend([lh_nh, lh_sh], 'NH', 'SH', 'location', 'southeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', [-20 10], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_lh', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            sh_nh=plot([1:12], -circshift(sh_lat_nh,shiftby_nh,2), 'color', par.orange);
            sh_sh=plot([1:12], -circshift(sh_lat_sh,shiftby_sh,2), '--', 'color', par.orange);
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{SH}$ (Wm$^{-2}$)'));
            legend([sh_nh, sh_sh], 'NH', 'SH', 'location', 'southeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', [-50 50], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_sh', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            rhsfc_nh=plot([1:12], circshift(rhsfc_lat_nh,shiftby_nh,2), 'k');
            rhsfc_sh=plot([1:12], circshift(rhsfc_lat_sh,shiftby_sh,2), '--k');
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{RH_{SFC}}$ (1)'));
            legend([rhsfc_nh, rhsfc_sh], 'NH', 'SH', 'location', 'southeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', [0 1], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_rhsfc', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            evap_nh=plot([1:12], circshift(evap_lat_nh,shiftby_nh,2), 'k');
            evap_sh=plot([1:12], circshift(evap_lat_sh,shiftby_sh,2), '--k');
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{E}$ (kg m$^{-2}$ s$^{-1}$)'));
            legend([evap_nh, evap_sh], 'NH', 'SH', 'location', 'northeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'ylim', 1e-6*[-1 6], 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_evap', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            snevap_nh=plot([1:12], circshift(snevap_lat_nh,shiftby_nh,2), 'k');
            snevap_sh=plot([1:12], circshift(snevap_lat_sh,shiftby_sh,2), '--k');
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{E_{snow}}$ (kg m$^{-2}$ s$^{-1}$)'));
            legend([snevap_nh, snevap_sh], 'NH', 'SH', 'location', 'northeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_snevap', folder), '-dpng', '-r300');
            close;

            % ALL lat x mon dependence of RCE and RAE
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            snmelt_nh=plot([1:12], circshift(snmelt_lat_nh,shiftby_nh,2), 'k');
            snmelt_sh=plot([1:12], circshift(snmelt_lat_sh,shiftby_sh,2), '--k');
            make_title_type(type);
            xlabel('Months (relative to winter)');
            ylabel(sprintf('$\\mathrm{Snowmelt}$ (kg m$^{-2}$ s$^{-1}$)'));
            legend([snmelt_nh, snmelt_sh], 'NH', 'SH', 'location', 'northeast');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_snmelt', folder), '-dpng', '-r300');
            close;

        end % for mse dse
    end % for land

end % for function
