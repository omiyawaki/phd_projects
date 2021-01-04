function plot_dmse_polar_line_topocomp(type, par)
    make_dirs(type, par)

    if strcmp(type, 'echam') & strcmp(par.echam.clim, 'echr0001')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        prefix_comp=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, 'echr0023');
        prefix_proc_comp=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, 'echr0023');
    else
        error('This function is only relevant for ECHAM AGCM simulation testing topography.');
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    tmp = load(sprintf('%s/%s/flux_z.mat', prefix_proc_comp, par.lat_interp)); flux_z_comp = tmp.flux_z; clear tmp; % load lat x mon RCAE data

    lat_bound_list = [80];

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
                latn_pole = 90; latn = lat_bound:dlat:latn_pole; monlabel=par.monlabel;
                clatn = cosd(latn); % cosine of latnitude for cosine weighting
                clatn_mon = repmat(clatn', [1 12]);
                lats_pole = -90; lats = -lat_bound:-dlat:lats_pole; monlabel=par.monlabelsh;
                clats = cosd(lats); % cosine of latsitude for cosine weighting
                clats_mon = repmat(clats', [1 12]);

                folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g_topocomp', par.plotdir, fw, land, -lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                ra = flux_z.(land).ra.(fw);
                ra_latn = interp1(grid.dim3.lat, ra, latn);
                ra_latn = nansum(ra_latn.*clatn_mon)/nansum(clatn);
                ra_lats = interp1(grid.dim3.lat, ra, lats);
                ra_lats = nansum(ra_lats.*clats_mon)/nansum(clats);
                res = flux_z.(land).res.(fw);
                res_latn = interp1(grid.dim3.lat, res, latn);
                res_latn = nansum(res_latn.*clatn_mon)/nansum(clatn);
                res_lats = interp1(grid.dim3.lat, res, lats);
                res_lats = nansum(res_lats.*clats_mon)/nansum(clats);
                stf = flux_z.(land).stf.(fw);
                stf_latn = interp1(grid.dim3.lat, stf, latn);
                stf_latn = nansum(stf_latn.*clatn_mon)/nansum(clatn);
                stf_lats = interp1(grid.dim3.lat, stf, lats);
                stf_lats = nansum(stf_lats.*clats_mon)/nansum(clats);

                comp_ra = flux_z_comp.(land).ra.(fw);
                comp_ra_latn = interp1(grid.dim3.lat, comp_ra, latn);
                comp_ra_latn = nansum(comp_ra_latn.*clatn_mon)/nansum(clatn);
                comp_ra_lats = interp1(grid.dim3.lat, comp_ra, lats);
                comp_ra_lats = nansum(comp_ra_lats.*clats_mon)/nansum(clats);
                comp_res = flux_z_comp.(land).res.(fw);
                comp_res_latn = interp1(grid.dim3.lat, comp_res, latn);
                comp_res_latn = nansum(comp_res_latn.*clatn_mon)/nansum(clatn);
                comp_res_lats = interp1(grid.dim3.lat, comp_res, lats);
                comp_res_lats = nansum(comp_res_lats.*clats_mon)/nansum(clats);
                comp_stf = flux_z_comp.(land).stf.(fw);
                comp_stf_latn = interp1(grid.dim3.lat, comp_stf, latn);
                comp_stf_latn = nansum(comp_stf_latn.*clatn_mon)/nansum(clatn);
                comp_stf_lats = interp1(grid.dim3.lat, comp_stf, lats);
                comp_stf_lats = nansum(comp_stf_lats.*clats_mon)/nansum(clats);

                [lh, sh] = rename_stf(type, flux_z, land);
                lh_latn = interp1(grid.dim3.lat, lh, latn);
                lh_latn = nansum(lh_latn.*clatn_mon)/nansum(clatn);
                lh_lats = interp1(grid.dim3.lat, lh, lats);
                lh_lats = nansum(lh_lats.*clats_mon)/nansum(clats);
                sh_latn = interp1(grid.dim3.lat, sh, latn);
                sh_latn = nansum(sh_latn.*clatn_mon)/nansum(clatn);
                sh_lats = interp1(grid.dim3.lat, sh, lats);
                sh_lats = nansum(sh_lats.*clats_mon)/nansum(clats);

                [comp_lh, comp_sh] = rename_stf(type, flux_z_comp, land);
                comp_lh_latn = interp1(grid.dim3.lat, comp_lh, latn);
                comp_lh_latn = nansum(comp_lh_latn.*clatn_mon)/nansum(clatn);
                comp_lh_lats = interp1(grid.dim3.lat, comp_lh, lats);
                comp_lh_lats = nansum(comp_lh_lats.*clats_mon)/nansum(clats);
                comp_sh_latn = interp1(grid.dim3.lat, comp_sh, latn);
                comp_sh_latn = nansum(comp_sh_latn.*clatn_mon)/nansum(clatn);
                comp_sh_lats = interp1(grid.dim3.lat, comp_sh, lats);
                comp_sh_lats = nansum(comp_sh_lats.*clats_mon)/nansum(clats);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lats ,6,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lats,6,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lats,6,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lats,6,2), 'color', par.orange);
                ra_f = plot([1:12], circshift(comp_ra_lats ,6,2),':', 'color', 0.5*[1 1 1]);
                res_f = plot([1:12], circshift(comp_res_lats,6,2), ':', 'color', par.maroon);
                lhf_f = plot([1:12], circshift(comp_lh_lats ,6,2),  ':', 'color', par.blue);
                shf_f = plot([1:12], circshift(comp_sh_lats ,6,2),  ':', 'color', par.orange);
                title(sprintf('%s AGCM, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound, lats_pole));
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res lhf shf ra_f res_f lhf_f shf_f], '$R_a^{C,S}$', '$\nabla\cdot F_m^{C,S}$', '$\mathrm{LH}^{C,S}$', '$\mathrm{SH}^{C,S}$','$R_a^{F,S}$', '$\nabla\cdot F_m^{F,S}$', '$\mathrm{LH}^{F,S}$', '$\mathrm{SH}^{F,S}$', 'location', 'eastoutside', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse_sh_comp', folder), '-dpng', '-r300');
                close;

                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12],  circshift(ra_lats ,6,2) - circshift(comp_ra_lats ,6,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(res_lats,6,2) - circshift(comp_res_lats,6,2), 'color', par.maroon);
                dlh=plot([1:12],  circshift(lh_lats,6,2) -  circshift(comp_lh_lats ,6,2), 'color', par.blue);
                dsh=plot([1:12],  circshift(sh_lats,6,2) -  circshift(comp_sh_lats ,6,2), 'color', par.orange);
                title(sprintf('%s AGCM, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound, lats_pole));
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dlh dsh], '$R_a^C - R_a^F$', '$\nabla\cdot F_m^C - \nabla\cdot F_m^F$', '$\mathrm{LH}^C-\mathrm{LH}^F$', '$\mathrm{SH}^C-\mathrm{SH}^F$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse_sh_comp', folder), '-dpng', '-r300');
                close;

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_latn ,0,2), '--', 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_latn,0,2), '--', 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_latn,0,2),  '--', 'color', par.blue);
                shf=plot([1:12], circshift(sh_latn,0,2),  '--', 'color', par.orange);
                ra_f=plot([1:12], circshift(comp_ra_lats ,6,2),':', 'color', 0.5*[1 1 1]);
                res_f=plot([1:12], circshift(comp_res_lats,6,2), ':', 'color', par.maroon);
                lhf_f=plot([1:12], circshift(comp_lh_lats ,6,2),  ':', 'color', par.blue);
                shf_f=plot([1:12], circshift(comp_sh_lats ,6,2),  ':', 'color', par.orange);
                title(sprintf('%s AGCM', upper(type)));
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res lhf shf ra_f res_f lhf_f shf_f], '$R_a^{C,N}$', '$\nabla\cdot F_m^{C,N}$', '$\mathrm{LH}^{C,N}$', '$\mathrm{SH}^{C,N}$','$R_a^{F,S}$', '$\nabla\cdot F_m^{F,S}$', '$\mathrm{LH}^{F,S}$', '$\mathrm{SH}^{F,S}$', 'location', 'eastoutside', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse_flatsh_toponh', folder), '-dpng', '-r300');
                close;

                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12],  circshift(comp_ra_lats ,6,2) - circshift(ra_latn ,0,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(comp_res_lats,6,2) - circshift(res_latn,0,2), 'color', par.maroon);
                dlh=plot([1:12],  circshift(comp_lh_lats,6,2) -  circshift(lh_latn ,0,2), 'color', par.blue);
                dsh=plot([1:12],  circshift(comp_sh_lats,6,2) -  circshift(sh_latn ,0,2), 'color', par.orange);
                title(sprintf('%s AGCM', upper(type)));
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dlh dsh], '$R_a^{F,S} - R_a^{C,N}$', '$\nabla\cdot F_m^{F,S} - \nabla\cdot F_m^{C,N}$', '$\mathrm{LH}^{F,S}-\mathrm{LH}^{C,N}$', '$\mathrm{SH}^{F,S}-\mathrm{SH}^{C,N}$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse_flatsh_toponh', folder), '-dpng', '-r300');
                close;

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

    %         folder = sprintf('%s/dmse/%s/0_poleward_of_lat_%g', par.plotdir, fw, lat_bound);
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
