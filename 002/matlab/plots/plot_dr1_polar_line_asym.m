function plot_dr1_polar_line_asym(type, par)
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
                lats_pole = -90; lats = -lat_bound:-dlat:lats_pole; monlabel=par.monlabelsh;
                clats = cosd(lats); % cosine of latsitude for cosine weighting
                clats_mon = repmat(clats', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_poleward_of_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                % R1 computed after zonal averaging
                r1z_latn = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), latn);
                r1z_latn = nansum(r1z_latn.*clatn_mon)/nansum(clatn);

                r1z_lats = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lats);
                r1z_lats = nansum(r1z_lats.*clats_mon)/nansum(clats);

                % hemispheric difference
                dr1z = circshift(r1z_lats,6,2) - r1z_latn;

                % fm and ra components
                fm_latn = interp1(grid.dim3.lat, flux_z.(land).res.(fw), latn);
                fm_latn = nansum(fm_latn.*clatn_mon)/nansum(clatn);
                fm_lats = interp1(grid.dim3.lat, flux_z.(land).res.(fw), lats);
                fm_lats = nansum(fm_lats.*clats_mon)/nansum(clats);
                dfm = circshift(fm_lats,6,2) - fm_latn;

                ra_latn = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), latn);
                ra_latn = nansum(ra_latn.*clatn_mon)/nansum(clatn);
                ra_lats = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), lats);
                ra_lats = nansum(ra_lats.*clats_mon)/nansum(clats);
                dra = circshift(ra_lats,6,2) - ra_latn;

                stf_latn = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), latn);
                stf_latn = nansum(stf_latn.*clatn_mon)/nansum(clatn);
                stf_lats = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), lats);
                stf_lats = nansum(stf_lats.*clats_mon)/nansum(clats);
                dstf = circshift(stf_lats,6,2) - stf_latn;

                lh_latn = interp1(grid.dim3.lat, flux_z.(land).hfls, latn);
                lh_latn = nansum(lh_latn.*clatn_mon)/nansum(clatn);
                lh_lats = interp1(grid.dim3.lat, flux_z.(land).hfls, lats);
                lh_lats = nansum(lh_lats.*clats_mon)/nansum(clats);
                dlh = circshift(lh_lats,6,2) - lh_latn;

                sh_latn = interp1(grid.dim3.lat, flux_z.(land).hfss, latn);
                sh_latn = nansum(sh_latn.*clatn_mon)/nansum(clatn);
                sh_lats = interp1(grid.dim3.lat, flux_z.(land).hfss, lats);
                sh_lats = nansum(sh_lats.*clats_mon)/nansum(clats);
                dsh = circshift(sh_lats,6,2) - sh_latn;

                comp1 = dfm./(ra_latn);
                pc2_latn = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw).^2, latn);
                pc2_latn = nansum(pc2_latn.*clatn_mon)/nansum(clatn);
                % comp2 = -fm_latn./(ra_latn).^2.*dra;
                comp2 = -pc2_latn.*dra;

                comp1a = dstf./(ra_latn);
                pc2a_latn = interp1(grid.dim3.lat, flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw).^2, latn);
                pc2a_latn = nansum(pc2a_latn.*clatn_mon)/nansum(clatn);
                % comp2a = -stf_latn./(ra_latn).^2.*dra;
                comp2a = -pc2a_latn.*dra;

                comp1b = dlh./(ra_latn);
                comp2b = dsh./(ra_latn);
                % comp3b = -stf_latn./(ra_latn).^2.*dra;
                comp3b = -pc2a_latn.*dra;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1z, 'k');
                res=plot([1:12], comp1+comp2-dr1z, '-.k');
                c1=plot([1:12], comp1, '-', 'color', par.maroon);
                c2=plot([1:12], comp2, '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                ylabel(sprintf('$R_1^{\\mathrm{S-N}}$ (unitless)'));
                legend([tot res c1 c2], '$R_{1,\,\mathrm{S-N}}$', 'Residual', '$\frac{\nabla\cdot F_m^{\mathrm{N}}}{R_a^{\mathrm{N}}} (\nabla\cdot F_m)^{\mathrm{S-N}}$', '$\Delta R_a$', 'location', 'eastoutside', 'numcolumns', 2);
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_asym', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z STF DECOMP lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1z, 'k');
                res=plot([1:12], comp1a+comp2a-dr1z, '-.k');
                c1=plot([1:12], comp1a, '-', 'color', par.blue);
                c2=plot([1:12], comp2a, '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                ylabel(sprintf('$R_1^{\\mathrm{S-N}}$ (unitless)'));
                legend([tot res c1 c2], '$R_{1,\,\mathrm{S-N}}$', 'Residual', '$\frac{(\mathrm{LH+SH})^{\mathrm{S-N}}}{R_a^{\mathrm{N}}}$', '$\Delta R_a$', 'location', 'eastoutside', 'numcolumns', 2);
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_asym_stf', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z lh+sh DECOMP lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1z, 'k');
                res=plot([1:12], comp1b+comp2b+comp3b-dr1z, '-.k');
                c1=plot([1:12], comp1b, '-', 'color', par.blue);
                c2=plot([1:12], comp2b, '-', 'color', par.orange);
                c3=plot([1:12], comp3b, '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s', par.gcm.clim));
                    else
                        title(sprintf('%s', par.model));
                    end
                end
                ylabel(sprintf('$R_1^{\\mathrm{S-N}}$ (unitless)'));
                legend([tot res c1 c2], '$R_{1,\,\mathrm{S-N}}$', 'Residual', '$\frac{(\mathrm{LH+SH})^{\mathrm{S-N}}}{R_a^{\mathrm{N}}}$', '$\Delta R_a$', 'location', 'eastoutside', 'numcolumns', 2);
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_asym_lh_sh', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
