function plot_dr1_polar_line_repl(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2')
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
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [80];

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'L+O';
        elseif strcmp(land, 'l'); land_text = 'L';
        elseif strcmp(land, 'o'); land_text = 'O';
        end

        if strcmp(type, 'echam')
            if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                latn_pole = 90; latn = lat_bound:dlat:latn_pole; monlabel=par.monlabel;
                clatn = cosd(latn); % cosine of latnitude for cosine weighting
                clatn_mon = repmat(clatn', [1 12]);
                lats_pole = -90; lats = -lat_bound:-dlat:lats_pole; % monlabel=par.monlabelsh;
                clats = cosd(lats); % cosine of latsitude for cosine weighting
                clats_mon = repmat(clats', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_poleward_of_lat_-%g_repl_n', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                % fm and ra components
                fm_latn = interp1(grid.dim3.lat, flux_z.(land).res.(fw), latn);
                fm_latn = nansum(fm_latn.*clatn_mon)/nansum(clatn);
                fm_lats = interp1(grid.dim3.lat, flux_z.(land).res.(fw), lats);
                fm_lats = nansum(fm_lats.*clats_mon)/nansum(clats);

                ra_latn = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), latn);
                ra_latn = nansum(ra_latn.*clatn_mon)/nansum(clatn);
                ra_lats = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), lats);
                ra_lats = nansum(ra_lats.*clats_mon)/nansum(clats);

                stf_latn = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), latn);
                stf_latn = nansum(stf_latn.*clatn_mon)/nansum(clatn);
                stf_lats = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), lats);
                stf_lats = nansum(stf_lats.*clats_mon)/nansum(clats);

                [lh, sh] = rename_stf(type, flux_z, land);

                lh_latn = interp1(grid.dim3.lat, lh, latn);
                lh_latn = nansum(lh_latn.*clatn_mon)/nansum(clatn);
                lh_lats = interp1(grid.dim3.lat, lh, lats);
                lh_lats = nansum(lh_lats.*clats_mon)/nansum(clats);

                sh_latn = interp1(grid.dim3.lat, sh, latn);
                sh_latn = nansum(sh_latn.*clatn_mon)/nansum(clatn);
                sh_lats = interp1(grid.dim3.lat, sh, lats);
                sh_lats = nansum(sh_lats.*clats_mon)/nansum(clats);

                % REPLACE SH Ra with NH
                r1z_latn = fm_latn./ra_latn;
                r1z_lats = fm_lats./ra_lats;

                r1z_lats_alt = (ra_lats + stf_lats)./ra_lats;

                r1z_lats_ra_n = (circshift(ra_latn,6,2) + lh_lats + sh_lats)./circshift(ra_latn,6,2);
                r1z_lats_stf_n = (ra_lats + circshift(stf_latn,6,2))./ra_lats;
                r1z_lats_lh_n = (ra_lats + circshift(lh_latn,6,2) + sh_lats)./ra_lats;
                r1z_lats_sh_n = (ra_lats + circshift(sh_latn,6,2) + lh_lats)./ra_lats;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                ylim_lo = 0.2;
                ylim_up = 2;
                raemin = par.ga;
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                all_s=plot([1:12], circshift(r1z_lats,6,2), 'k');
                all_n=plot([1:12], r1z_latn, '--k');
                plot([1:12], circshift(r1z_lats_alt,6,2), ':k');
                % north=plot([1:12], r1z_latn, '--k');
                ra_n=plot([1:12], circshift(r1z_lats_ra_n,6,2), 'color', 0.5*[1 1 1]);
                % stf_n=plot([1:12], circshift(r1z_lats_stf_n,6,2), '-.k');
                lh_n=plot([1:12], circshift(r1z_lats_lh_n,6,2), 'color', par.blue);
                sh_n=plot([1:12], circshift(r1z_lats_sh_n,6,2), 'color', par.orange);
                xlabel('Month (relative to winter solstice)');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s', upper(type), echamtext));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                % xlabel('Month');
                legend([all_s all_n ra_n lh_n sh_n], '$R_1^S$', '$R_1^N$', '$\frac{R_a^N + \mathrm{LH}^S+\mathrm{SH}^S}{R_a^N}$', '$\frac{R_a^S+\mathrm{LH}^N+\mathrm{SH}^S}{R_a^S}$', '$\frac{R_a^S+\mathrm{LH}^S+\mathrm{SH}^N}{R_a^S}$', 'location', 'eastoutside', 'numcolumns', 1);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                print(sprintf('%s/0_mon_dr1z_repl', folder), '-dpng', '-r300');
                close;

                r1z_ann_latn = nanmean(r1z_latn,2);
                r1z_ann_lats = nanmean(r1z_lats,2);
                lh_ann_latn = nanmean(lh_latn,2);
                sh_ann_latn = nanmean(sh_latn,2);
                ra_ann_latn = nanmean(ra_latn,2);
                lh_ann_lats = nanmean(lh_lats,2);
                sh_ann_lats = nanmean(sh_lats,2);
                ra_ann_lats = nanmean(ra_lats,2);
                comp_ann_ra = -(lh_ann_latn+sh_ann_latn)./ra_ann_latn.^2.*(ra_ann_lats-ra_ann_latn);
                comp_ann_lh = (lh_ann_lats-lh_ann_latn)./ra_ann_latn;
                comp_ann_sh = (sh_ann_lats-sh_ann_latn)./ra_ann_latn;
                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                % line([0.5,4.5],[0,0], 'linewidth', 0.5, 'color', 'k');
                % ylim_lo = 0.2;
                % ylim_up = 2;
                % raemin = par.ga;
                % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                dr1=bar(1, r1z_ann_lats - r1z_ann_latn, 'k');
                dra=bar(2, comp_ann_ra, 'facecolor', 0.5*[1 1 1]);
                dlh=bar(3, comp_ann_lh, 'facecolor', par.blue);
                dsh=bar(4, comp_ann_sh, 'facecolor', par.orange);
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xtick', [1:4], 'xticklabel', {'$R_1$','$R_a$','LH','SH'}, 'yminortick', 'on', 'tickdir', 'out');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s', upper(type), echamtext));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                % xlabel('Month');
                % legend([all_s all_n ra_n lh_n sh_n], '$R_1^S$', '$R_1^N$', '$\frac{R_a^N + \mathrm{NH}^S+\mathrm{SH}^S}{R_a^N}$', '$\frac{R_a^S+\mathrm{LH}^N+\mathrm{SH}^S}{R_a^S}$', '$\frac{R_a^S+\mathrm{LH}^S+\mathrm{SH}^N}{R_a^S}$', 'location', 'eastoutside', 'numcolumns', 1);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                print(sprintf('%s/0_mon_dr1z_ann_decomp', folder), '-dpng', '-r300');
                close;

                comp_ra = -(lh_latn+sh_latn)./ra_latn.^2.*(circshift(ra_lats,6,2)-ra_latn);
                comp_lh = (circshift(lh_lats,6,2)-lh_latn)./ra_latn;
                comp_sh = (circshift(sh_lats,6,2)-sh_latn)./ra_latn;
                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                line([1,12],[0,0], 'linewidth', 0.5, 'color', 'k');
                % ylim_lo = 0.2;
                % ylim_up = 2;
                % raemin = par.ga;
                % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                dr1=plot([1:12], circshift(r1z_lats,6,2) - r1z_latn, 'k');
                res=plot([1:12], (circshift(r1z_lats,6,2)-r1z_latn)-(comp_ra+comp_lh+comp_sh), '-.k');
                % plot([1:12], circshift(r1z_lats_alt,6,2), ':k');
                % % north=plot([1:12], r1z_latn, '--k');
                dra_n=plot([1:12], comp_ra, 'color', 0.5*[1 1 1]);
                % % stf_n=plot([1:12], circshift(r1z_lats_stf_n,6,2), '-.k');
                dlh_n=plot([1:12], comp_lh, 'color', par.blue);
                dsh_n=plot([1:12], comp_sh, 'color', par.orange);
                xlabel('Month (relative to winter solstice)');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s', upper(type)));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s', upper(type), echamtext));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                % xlabel('Month');
                legend([dr1 res dra_n dlh_n dsh_n], '$R_1^S-R_1^N$', 'Residual', '$-\frac{\mathrm{LH}^N+\mathrm{SH}^N}{(R_a^N)^2}(R_a^S-R_a^N)$', '$\frac{1}{R_a^N}(\mathrm{LH}^S-\mathrm{LH}^N)$', '$\frac{1}{R_a^N}(\mathrm{SH}^S-\mathrm{SH}^N)$', 'location', 'eastoutside', 'numcolumns', 1);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                print(sprintf('%s/0_mon_dr1z_decomp', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
