function plot_dmse_toasfc_midlatitude_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [5 10 15 20 -5 -10 -15 -20];
    lat_bound_list = [5 -5];

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center;
                else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dmse_toasfc/%s/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                rtoa = flux_z.(land).rtoa;
                rtoa_lat = interp1(grid.dim3.lat, rtoa, lat);
                rtoa_lat = nansum(rtoa_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                sfc = flux_z.(land).sfc.(fw);
                sfc_lat = interp1(grid.dim3.lat, sfc, lat);
                sfc_lat = nansum(sfc_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                sfc=plot([1:12], sfc_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([rtoa res sfc], '$F_{\mathrm{TOA}}$', '$\nabla\cdot F_m$', '$F_\mathrm{SFC}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_toa_sfc', folder), '-dpng', '-r300');
                close;

                % NO LEGEND ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                sfc=plot([1:12], sfc_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_toa_sfc_noleg', folder), '-dpng', '-r300');
                close;

                sw = flux_z.(land).sw;
                sw_lat = interp1(grid.dim3.lat, sw, lat);
                sw_lat = nansum(sw_lat.*clat_mon)/nansum(clat);
                olr = flux_z.(land).olr;
                olr_lat = interp1(grid.dim3.lat, olr, lat);
                olr_lat = nansum(olr_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                shf = flux_z.(land).shf.(fw);
                shf_lat = interp1(grid.dim3.lat, shf, lat);
                shf_lat = nansum(shf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                sw=plot([1:12], sw_lat, 'color', par.yellow);
                olr=plot([1:12], olr_lat, 'color', par.green);
                res=plot([1:12], res_lat, 'color', par.maroon);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([sw olr res shf], '$\mathrm{SWABS}$', '$\mathrm{OLR}$', '$\nabla\cdot F_m$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swabs_olr_shf', folder), '-dpng', '-r300');
                close;

                sw_ann = repmat(nanmean(flux_z.(land).sw, 2), [1 12]);
                olr_ann = repmat(nanmean(flux_z.(land).olr, 2), [1 12]);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                shf_ann = repmat(nanmean(flux_z.(land).shf.(fw), 2), [1 12]);

                dsw = flux_z.(land).sw - sw_ann;
                dsw_lat = interp1(grid.dim3.lat, dsw, lat);
                dsw_lat = nansum(dsw_lat.*clat_mon)/nansum(clat);
                dolr = flux_z.(land).olr - olr_ann;
                dolr_lat = interp1(grid.dim3.lat, dolr, lat);
                dolr_lat = nansum(dolr_lat.*clat_mon)/nansum(clat);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat);
                dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
                dshf = flux_z.(land).shf.(fw) - shf_ann;
                dshf_lat = interp1(grid.dim3.lat, dshf, lat);
                dshf_lat = nansum(dshf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dsw=plot([1:12],  dsw_lat, 'color', par.yellow);
                dolr=plot([1:12], dolr_lat, 'color', par.green);
                dres=plot([1:12], dres_lat, 'color', par.maroon);
                dshf=plot([1:12], dshf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dsw dolr dres dshf], '$\mathrm{\Delta SWABS}$', '$\mathrm{\Delta OLR}$', '$\Delta(\nabla\cdot F_m)$', '$\Delta F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dswabs_dolr_dshf', folder), '-dpng', '-r300');
                close;

                rtoa = flux_z.(land).rtoa;
                rtoa_lat = interp1(grid.dim3.lat, rtoa, lat);
                rtoa_lat = nansum(rtoa_lat.*clat_mon)/nansum(clat);
                swsfc = flux_z.(land).swsfc;
                swsfc_lat = interp1(grid.dim3.lat, swsfc, lat);
                swsfc_lat = nansum(swsfc_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                shf = flux_z.(land).shf.(fw);
                shf_lat = interp1(grid.dim3.lat, shf, lat);
                shf_lat = nansum(shf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'color', 'k');
                swsfc=plot([1:12], swsfc_lat, 'color', par.yellow);
                res=plot([1:12], res_lat, 'color', par.maroon);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([rtoa swsfc res shf], '$F_{\mathrm{TOA}}$', '$F_\mathrm{SW,\,SFC}$', '$\nabla\cdot F_m$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_rtoa_swsfc_shf', folder), '-dpng', '-r300');
                close;

                % just SWSFC and SHF lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', par.yellow);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([swsfc shf], '$F_\mathrm{SW,\,SFC}$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_shf', folder), '-dpng', '-r300');
                close;

                % NOLEGEND just SWSFC and SHF lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', par.yellow);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_shf_noleg', folder), '-dpng', '-r300');
                close;

                lwsfc = flux_z.(land).lwsfc;
                lwsfc_lat = interp1(grid.dim3.lat, lwsfc, lat);
                lwsfc_lat = nansum(lwsfc_lat.*clat_mon)/nansum(clat);

                % SWSFC and LWSFC lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', 'b');
                lwsfc=plot([1:12], lwsfc_lat, 'color', 'r');
                max(lwsfc_lat)
                min(lwsfc_lat)
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([swsfc lwsfc], '$F_\mathrm{SW,\,SFC}$', '$F_\mathrm{LW,\,SFC}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_lwsfc', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
