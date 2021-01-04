function plot_dr1_line(type, par)
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
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_eval_list = [-85 -45 0 45 85]; % Latitude to evaluate R1 seasonality

    f_vec = assign_fw(type, par);
    for f = f_vec; fw = f{1};
        for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
            folder = sprintf('%s/dr1/%s/0_lat_%g', par.plotdir, fw, lat_eval);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat_eval);

            dr1 = flux_z.lo.r1.(fw) - r1_ann;
            dr1_lat = interp1(grid.dim3.lat, dr1, lat_eval);

            r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);

            r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);

            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]);
            ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]);
            figure(); clf; hold all; box on;
            if abs(lat_eval)<60
                rcemax = par.ep-r1_ann_lat(1);
                vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
            else
                raemin = par.ga-r1_ann_lat(1);
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            end
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dr1_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
            legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\Delta R_1$ (unitless)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
            close;
        end
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
                folder = sprintf('%s/dr1/%s/%s/0_lat_%g', par.plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat_eval);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat_eval);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1 = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);

                comp2_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2 = fm_ann./(ra_ann).^2.*delta_stf;
                comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo); ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up); ylim_up = 1; end;
                figure(); clf; hold all; box on;
                if abs(lat_eval)<60
                    % rce=plot([1:12], par.ep-r1_ann_lat, 'color', par.orange);
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                else
                    % rae=plot([1:12], par.ga-r1_ann_lat, 'color', par.blue);
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper(type), var_text, land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', par.model, var_text, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                if abs(lat_eval)<60
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                else
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                c12=plot([1:12], comp1_lat+comp2_lat, '-.k');
                c1=plot([1:12], comp1_lat, '--k');
                c2=plot([1:12], comp2_lat, ':k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper(type), var_text, land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', par.model, var_text, land_text, lat_eval)); end;
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta (LH + SH)$', 'location', 'eastoutside');
                xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land
end % for function
