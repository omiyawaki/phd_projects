function plot_dmse_polar_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    % load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

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

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat);
                ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                stf = flux_z.(land).stf.(fw);
                stf_lat = interp1(grid.dim3.lat, stf, lat);
                stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
                [lh, sh] = rename_stf(type, flux_z, land);
                lh_lat = interp1(grid.dim3.lat, lh, lat);
                lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                sh_lat = interp1(grid.dim3.lat, sh, lat);
                sh_lat = nansum(sh_lat.*clat_mon)/nansum(clat);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                make_title_type_lat(type, lat_bound, lat_pole);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res lhf shf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                % NOLEG ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                make_title_type_lat(type, lat_bound, lat_pole);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                % legend([ra res stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse_noleg', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
                ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat);
                dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
                res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat);
                dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
                stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
                dstf = flux_z.(land).stf.(fw) - stf_ann;
                dstf_lat = interp1(grid.dim3.lat, dstf, lat);
                dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12],  circshift(dra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(dres_lat,shiftby,2), 'color', par.maroon);
                dstf=plot([1:12], circshift(dstf_lat,shiftby,2), 'color', par.blue);
                make_title_type_lat(type, lat_bound, lat_pole);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse', folder), '-dpng', '-r300');
                close;

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat);
                ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
                sw = flux_z.(land).sw;
                sw_lat = interp1(grid.dim3.lat, sw, lat);
                sw_lat = nansum(sw_lat.*clat_mon)/nansum(clat);
                lw = flux_z.(land).lw;
                lw_lat = interp1(grid.dim3.lat, lw, lat);
                lw_lat = nansum(lw_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
                sw=plot([1:12], circshift(sw_lat,shiftby,2), 'color', par.blue);
                lw=plot([1:12], circshift(lw_lat,shiftby,2), 'color', par.maroon);
                make_title_type_lat(type, lat_bound, lat_pole);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra sw lw], '$R_a$', '$\mathrm{Net SW}$', '$\mathrm{Net LW}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_ra', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
                ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat);
                dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
                sw_ann = repmat(nanmean(flux_z.(land).sw, 2), [1 12]);
                sw_ann_lat = interp1(grid.dim3.lat, sw_ann, lat);
                sw_ann_lat = nansum(sw_ann_lat.*clat_mon)/nansum(clat);
                dsw = flux_z.(land).sw - sw_ann;
                dsw_lat = interp1(grid.dim3.lat, dsw, lat);
                dsw_lat = nansum(dsw_lat.*clat_mon)/nansum(clat);
                lw_ann = repmat(nanmean(flux_z.(land).lw, 2), [1 12]);
                lw_ann_lat = interp1(grid.dim3.lat, lw_ann, lat);
                lw_ann_lat = nansum(lw_ann_lat.*clat_mon)/nansum(clat);
                dlw = flux_z.(land).lw - lw_ann;
                dlw_lat = interp1(grid.dim3.lat, dlw, lat);
                dlw_lat = nansum(dlw_lat.*clat_mon)/nansum(clat);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12], circshift(dra_lat,shiftby,2), 'color', 0.5*[1 1 1]);
                dsw=plot([1:12], circshift(dsw_lat,shiftby,2), 'color', par.blue);
                dlw=plot([1:12], circshift(dlw_lat,shiftby,2), 'color', par.maroon);
                make_title_type_lat(type, lat_bound, lat_pole);
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dsw dlw], '$\Delta R_a$', '$\Delta(\mathrm{Net SW})$', '$\Delta (\mathrm{Net LW})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dra', folder), '-dpng', '-r300');
                close;

                if strcmp(type, 'gcm')
                    stf = flux_z.(land).stf.(fw);
                    stf_lat = interp1(grid.dim3.lat, stf, lat);
                    stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
                    hfls = flux_z.(land).hfls;
                    hfls_lat = interp1(grid.dim3.lat, hfls, lat);
                    hfls_lat = nansum(hfls_lat.*clat_mon)/nansum(clat);
                    hfss = flux_z.(land).hfss;
                    hfss_lat = interp1(grid.dim3.lat, hfss, lat);
                    hfss_lat = nansum(hfss_lat.*clat_mon)/nansum(clat);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], circshift(stf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(hfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(hfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\mathrm{LH+SH}$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_stf', folder), '-dpng', '-r300');
                    close;

                    stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                    stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
                    stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
                    dstf = flux_z.(land).stf.(fw) - stf_ann;
                    dstf_lat = interp1(grid.dim3.lat, dstf, lat);
                    dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
                    hfls_ann = repmat(nanmean(flux_z.(land).hfls, 2), [1 12]);
                    hfls_ann_lat = interp1(grid.dim3.lat, hfls_ann, lat);
                    hfls_ann_lat = nansum(hfls_ann_lat.*clat_mon)/nansum(clat);
                    dhfls = flux_z.(land).hfls - hfls_ann;
                    dhfls_lat = interp1(grid.dim3.lat, dhfls, lat);
                    dhfls_lat = nansum(dhfls_lat.*clat_mon)/nansum(clat);
                    hfss_ann = repmat(nanmean(flux_z.(land).hfss, 2), [1 12]);
                    hfss_ann_lat = interp1(grid.dim3.lat, hfss_ann, lat);
                    hfss_ann_lat = nansum(hfss_ann_lat.*clat_mon)/nansum(clat);
                    dhfss = flux_z.(land).hfss - hfss_ann;
                    dhfss_lat = interp1(grid.dim3.lat, dhfss, lat);
                    dhfss_lat = nansum(dhfss_lat.*clat_mon)/nansum(clat);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], circshift(dstf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(dhfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(dhfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole));
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\Delta(\mathrm{LH+SH})$', '$\Delta\mathrm{LH}$', '$\Delta\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dstf', folder), '-dpng', '-r300');
                    close;

                end

            end

        end % for mse dse
    end % for land

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    % elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
    %         dlat = 0.25; % step size for standard lat grid
    %         if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby=0; monlabel=par.monlabel;
    %         else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby=6; monlabel=par.monlabelsh; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dmse/%s/0_poleward_of_lat_%g', plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         ra = flux_z.lo.ra.(fw);
    %         ra_lat = interp1(grid.dim3.lat, ra, lat);
    %         ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.ra.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.ra.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RA lat x mon dependence of RCE and RAE
    %         var_text = '$R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$R_a$', '$R_{a,\mathrm{\,L+O}}$', '$R_{a,\mathrm{\,L}}$', '$R_{a,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_ra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         ra_ann = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
    %         ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
    %         dra = flux_z.lo.ra.(fw) - ra_ann;
    %         dra_lat = interp1(grid.dim3.lat, dra, lat);
    %         dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
    %         ra_ann_l = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.ra.(fw) - ra_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         ra_ann_o = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.ra.(fw) - ra_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RA lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(R_a)$', '$\Delta(R_{a,\mathrm{\,L+O}})$', '$\Delta(R_{a,\mathrm{\,L}})$', '$\Delta(R_{a,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res = flux_z.lo.res.(fw);
    %         res_lat = interp1(grid.dim3.lat, res, lat);
    %         res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.res.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.res.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RES lat x mon dependence of RCE and RAE
    %         var_text = '$\nabla\cdot F_m$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(res_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\nabla\cdot F_m$', '$\nabla\cdot F_{m,\mathrm{\,L+O}}$', '$\nabla\cdot F_{m,\mathrm{\,L}}$', '$\nabla\cdot F_{m,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\nabla\\cdot F_m$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_res_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res_ann = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
    %         res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
    %         dres = flux_z.lo.res.(fw) - res_ann;
    %         dres_lat = interp1(grid.dim3.lat, dres, lat);
    %         dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
    %         res_ann_l = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.res.(fw) - res_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         res_ann_o = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.res.(fw) - res_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RES lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta(\nabla\cdot F_m)$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dres_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\nabla\cdot F_m)$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L+O}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta(\\nabla\\cdot F_m)$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dres_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf = flux_z.lo.stf.(fw);
    %         stf_lat = interp1(grid.dim3.lat, stf, lat);
    %         stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.stf.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.stf.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % STF lat x mon dependence of RCE and STFE
    %         var_text = '$\mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(stf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\mathrm{LH+SH}$', '$(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$(\mathrm{LH+SH})_{\mathrm{L}}$', '$(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_stf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf_ann = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
    %         stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
    %         dstf = flux_z.lo.stf.(fw) - stf_ann;
    %         dstf_lat = interp1(grid.dim3.lat, dstf, lat);
    %         dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
    %         stf_ann_l = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.stf.(fw) - stf_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         stf_ann_o = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.stf.(fw) - stf_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA STF lat x mon dependence of RCE and STFE
    %         var_text = '$\Delta \mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dstf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\mathrm{LH+SH})$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\delta \\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dstf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %     end
    % end

end % for function
