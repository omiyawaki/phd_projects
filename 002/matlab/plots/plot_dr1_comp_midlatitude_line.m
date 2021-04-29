function plot_dr1_comp_midlatitude_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % read grid data

    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    par.lat_bound_list = [-10 10];
    center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_midlatitude_lat(center, par);
        
        load(sprintf('%s/dr1_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound));
        if strcmp(type, 'rea') %| (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
            load(sprintf('%s/dr2_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound));
        end
        
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam')
                if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
            end

            for f = 1:length(par.(type).fw); fw = par.(type).fw{f};
                r1z.(fw) = flux_z.(land).res.(fw) ./ flux_z.(land).ra.(fw);
                r1z.(fw) = interp1(grid.dim2.lat, r1z.(fw), lat);
                r1z_lat.(fw) = nansum(clat_mon.*r1z.(fw))/nansum(clat);
            end

            f_vec = assign_fw(type, par);
            
            par.folder = sprintf('%s/dr1/mse/%s/0_midlatitude_lat_%g_to_%g', plotdir, land, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
            if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

            var_text = '$\Delta R_1$';
            figure(); clf; hold all; box on;
            
            ylim_lo = -0.6;
            ylim_up = 1.6;
            rcemax = par.ep;
            if rcemax > ylim_lo
                vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
            end
            raemin = par.ga;
            if raemin < ylim_up
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            end
            mse_old=plot([1:12], circshift(r1z_lat.mse_old,par.shiftby,2), 'k');
            % mse_lat=plot([1:12], circshift(r1z_lat.mse_lat,par.shiftby,2), 'color', par.blue);
            mse=plot([1:12], circshift(r1z_lat.mse,par.shiftby,2), '--k');
            if isfield(par,'lat_center')
                make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
            else
                make_title_type_lat(type, par.lat_bound, par.lat_pole, par);
            end
            ylabel(sprintf('$R_1$ (unitless)'));
            % legend([mse_old mse_lat mse], '$\frac{\partial_t h + \nabla\cdot F_m}{R_a}$', '$\frac{\partial_t (c_p T) + \nabla\cdot F_m}{R_a}$', '$\frac{\nabla\cdot F_m}{R_a}$', 'location', 'eastoutside', 'orientation', 'vertical');
            legend([mse_old mse], '$\frac{\partial_t h + \nabla\cdot F_m}{R_a}$', '$\frac{\nabla\cdot F_m}{R_a}$', 'location', 'eastoutside', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_r1z_comp', par.folder);
            print(plotname, '-dpng', '-r300');
            close;

        end % for land
    end % lat bound

end % for function
