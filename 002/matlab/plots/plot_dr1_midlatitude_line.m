function plot_dr1_midlatitude_line(type, par)
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

    lat_bound_list = [-5 5];

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
                if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
                else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                % R1 computed before zonal averaging
                r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                r1_lat = nansum(r1_lat.*clat_mon)/nansum(clat);
                comp1r_lat = interp1(grid.dim3.lat, flux_z.(land).comp1.(fw), lat);
                comp1r_lat = nansum(comp1r_lat.*clat_mon)/nansum(clat);
                comp2r_lat = interp1(grid.dim3.lat, flux_z.(land).comp2.(fw), lat);
                comp2r_lat = nansum(comp2r_lat.*clat_mon)/nansum(clat);

                % R1 computed after zonal averaging
                r1z_lat = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                r1z_lat = nansum(r1z_lat.*clat_mon)/nansum(clat);

                % R1 lat x mon dependence of RCE and RAE
                var_text = '$R_1$';
                figure(); clf; hold all; box on;
                ylim_lo = -0.5;
                ylim_up = 1.5;
                    rcemax = par.ep;
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    raemin = par.ga;
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([1 12], [1 1], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(r1_lat,shiftby, 2), 'k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_r1', folder), '-dpng', '-r300');
                close;

                % R1 computed at each lat x lon RES and RA
                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                r1z_ann = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat = interp1(grid.dim3.lat, r1z_ann, lat);
                r1z_ann_lat = nansum(r1z_ann_lat.*clat_mon)/nansum(clat);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                dr1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - r1z_ann;
                dr1z_lat = interp1(grid.dim3.lat, dr1z, lat);
                dr1z_lat = nansum(dr1z_lat.*clat_mon)/nansum(clat);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1a = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1a_lat = interp1(grid.dim3.lat, comp1a, lat);
                comp1a_lat = nansum(comp1a_lat.*clat_mon)/nansum(clat);

                comp2a_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2a = fm_ann./(ra_ann).^2.*delta_stf;
                comp2a_lat = interp1(grid.dim3.lat, comp2a, lat);
                comp2a_lat = nansum(comp2a_lat.*clat_mon)/nansum(clat);

                % % DELTA R1 lat x mon dependence of RCE and RAE
                % var_text = '$\Delta R_1$';
                % figure(); clf; hold all; box on;
                % colororder({'k', 'k'});
                % yyaxis left
                % ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1_ann_lat(1)+max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                % ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % yyaxis right
                % ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % rcemax = par.ep-r1_ann_lat(1);
                % if rcemax > ylim_lo
                %     vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % raemin = par.ga-r1_ann_lat(1);
                % if raemin < ylim_up
                %     vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                % tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                % elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                % ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                % close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                c12=plot([1:12], circshift(comp1a_lat+comp2a_lat, shiftby, 2), '-.k');
                c1=plot([1:12],  circshift(comp1a_lat, shiftby, 2), '--k');
                c2=plot([1:12],  circshift(comp2a_lat, shiftby, 2), ':k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta (LH + SH)$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_alt', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                comp1s = delta_fm./ra_ann;
                comp1s_lat = interp1(grid.dim3.lat, comp1s, lat);
                comp1s_lat = nansum(comp1s_lat.*clat_mon)/nansum(clat);

                delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                comp2s = -fm_ann./(ra_ann).^2.*delta_ra;
                comp2s_lat = interp1(grid.dim3.lat, comp2s, lat);
                comp2s_lat = nansum(comp2s_lat.*clat_mon)/nansum(clat);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'southoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp', folder), '-dpng', '-r300');
                close;

                % NOLEGEND DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_noleg', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -0.2;
                ylim_up = 0.5;
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -r1z_ann_lat(1)+ylim_lo;
                ylim_up = -r1z_ann_lat(1)+ylim_up;
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2) - circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext , -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside', 'NumColumns', 2);
                legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'southoutside', 'orientation', 'horizontal');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = min([r1z_lat,-0.4]);
                ylim_up = max([r1z_lat,0.6]);
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -r1z_ann_lat(1)+ylim_lo;
                ylim_up = -r1z_ann_lat(1)+ylim_up;
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2) - circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside', 'NumColumns', 2);
                legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'southoutside', 'orientation', 'horizontal', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_2x2_legend', folder), '-dpng', '-r300');
                close;

                % NO LEGEND DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = min([r1z_lat,-0.4]);
                ylim_up = max([r1z_lat,0.6]);
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -r1z_ann_lat(1)+ylim_lo;
                ylim_up = -r1z_ann_lat(1)+ylim_up;
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2) - circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), echamtext, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm');
                    if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.gcm.clim, -lat_bound+lat_center, lat_bound+lat_center));
                    else
                        title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center));
                    end
                end
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside');
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_noleg', folder), '-dpng', '-r300');
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
    %         if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
    %         else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dr1/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
    %         r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

    %         dr1 = flux_z.lo.r1.(fw) - r1_ann;
    %         dr1_lat = interp1(grid.dim3.lat, dr1, lat);
    %         dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

    %         % DELTA R1 lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_1$';
    %         ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
    %         ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
    %         figure(); clf; hold all; box on;
    %             rcemax = par.ep-r1_ann_lat(1);
    %             vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
    %             patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    %             raemin = par.ga-r1_ann_lat(1);
    %             vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    %             patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat, shiftby, 2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat, shiftby, 2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat, shiftby, 2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound+lat_center, lat_bound+lat_center));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound+lat_center, lat_bound+lat_center));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
    %         legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_1$ (unitless)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
    %         close;
    %     end
    % end

end % for function
