function plot_dr1_midlatitude_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
    load(sprintf('%s/flux_t.mat', prefix_proc)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    par.lat_bound_list = [-10 10];
    center = 50;

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
        
            for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);
            
                dlat = 0.25; % step size for standard lat grid
                if par.lat_bound>0; par.lat_center=center; lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center; par.shiftby=0; par.monlabel=par.monlabelnh;
                else; par.lat_center=-center; lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center; par.shiftby=6; par.monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                par.folder = sprintf('%s/dr1/%s/%s/0_midlatitude_lat_%g_to_%g', plotdir, fw, land, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
                if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

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

                % R1 computed at each lat x lon RES and RA
                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                r1z_ann = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat = interp1(grid.dim3.lat, r1z_ann, lat);
                r1z_ann_lat = nansum(r1z_ann_lat.*clat_mon)/nansum(clat);

                % annual mean values
                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                % compute deviation from annual mean
                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                dr1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - r1z_ann;
                dr1z_lat = interp1(grid.dim3.lat, dr1z, lat);
                dr1z_lat = nansum(dr1z_lat.*clat_mon)/nansum(clat);

                % DIV FM and STF DECOMP
                % compute mse flux divergence component
                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1a = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1a_lat = interp1(grid.dim3.lat, comp1a, lat);
                comp1a_lat = nansum(comp1a_lat.*clat_mon)/nansum(clat);

                comp2a_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2a = fm_ann./(ra_ann).^2.*delta_stf;
                comp2a_lat = interp1(grid.dim3.lat, comp2a, lat);
                comp2a_lat = nansum(comp2a_lat.*clat_mon)/nansum(clat);

                % DIVFM and RA DECOMP
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                comp1s = delta_fm./ra_ann;
                comp1s_lat = interp1(grid.dim3.lat, comp1s, lat);
                comp1s_lat = nansum(comp1s_lat.*clat_mon)/nansum(clat);

                delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                comp2s = -fm_ann./(ra_ann).^2.*delta_ra;
                comp2s_lat = interp1(grid.dim3.lat, comp2s, lat);
                comp2s_lat = nansum(comp2s_lat.*clat_mon)/nansum(clat);

                % Set y axis limits of plots
                if strcmp(type, 'echam')
                    ymin = -0.8;
                    ymax = 0.8;
                else
                    ymin = -0.6;
                    ymax = 0.8;
                end

                % MAKE PLOTS
                plot_dr1(r1z_lat, r1z_ann_lat, dr1z_lat, comp1s_lat, comp2s_lat, '', ymin, ymax, type, fw, par)
                plot_dr1(r1z_lat, r1z_ann_lat, dr1z_lat, comp1s_lat, comp2s_lat, '_noleg', ymin, ymax, type, fw, par)

            end

        end % for mse dse
    end % for land

end % for function

% plot dr1 and its decomposition
function plot_dr1(r1_var, r1_ann_var, dr1_var, comp1, comp2, leg, ymin, ymax, type, fw, par)

        var_text = '$\Delta R_1$';
        figure(); clf; hold all; box on;
        
        colororder({'k', 'k'});
        
        yyaxis left
        ylim_lo = min([r1_var,ymin]);
        ylim_up = max([r1_var,ymax]);
        tot=plot([1:12], circshift(r1_var,par.shiftby,2), 'k');
        ylabel(sprintf('$R_1$ (unitless)'));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
        
        yyaxis right
        ylim_lo = -r1_ann_var(1)+ylim_lo;
        ylim_up = -r1_ann_var(1)+ylim_up;
        rcemax = par.ep-r1_ann_var(1);
        if rcemax > ylim_lo
            vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
            patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
        end
        raemin = par.ga-r1_ann_var(1);
        if raemin < ylim_up
            vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
            patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
        end
        
        line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
        tot=plot([1:12], circshift(dr1_var,par.shiftby,2), 'k');
        res=plot([1:12], circshift(dr1_var,par.shiftby,2) - circshift(comp1+comp2,par.shiftby,2), '-.k');
        c1=plot([1:12],  circshift(comp1,par.shiftby,2), '-', 'color', par.maroon);
        c2=plot([1:12],  circshift(comp2,par.shiftby,2), '-', 'color', 0.5*[1 1 1]);
        
        make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
        ylabel(sprintf('$\\Delta R_1$ (unitless)'));
        
        if leg == ""
            if ~strcmp(fw, 'mse_old')
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside', 'orientation', 'vertical');
            else
                ylabel(sprintf('$\\Delta R_1^*$ (unitless)'));
                legend([tot res c1 c2], '$\Delta R_1^*$', 'Residual', '$\frac{\Delta (\partial_t h + \nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\partial_t h + \nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside', 'orientation', 'vertical');
            end
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
        elseif leg == "_noleg"
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        end
        
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/0_mon_dr1z_decomp%s', par.folder, leg), '-dpng', '-r300');
        close;

end
