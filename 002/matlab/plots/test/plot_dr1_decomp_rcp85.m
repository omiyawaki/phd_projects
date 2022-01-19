function plot_dr1_decomp_rcp85(type, par)
    
    if ~strcmp(par.gcm.clim, 'historical')
        error('This analysis must be run with the GCM historical setting.')
    end

    make_dirs(type, par)

    % LOAD HISTORICAL DATA
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    % LOAD RCP85 DATA
    par2 = par;
    par2.gcm.clim = 'rcp85';
    par2.gcm.yr_span = '207001-209912';
    prefix2 = make_prefix(type, par2);
    prefix_proc2 = make_prefix_proc(type, par2);

    tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux_z=tmp.flux_z; clear tmp; % read grid data

    tmp=load(sprintf('%s/grid.mat', prefix2)); grid2=tmp.grid; clear tmp; % read grid data
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc2)); flux_z2=tmp.flux_z; clear tmp; % read grid data

    lat_bound_list = [-80 80];

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] =  make_polar_lat(par);
        
        load(sprintf('%s/dr1_poleward_of_lat_%g.mat', prefix_proc, par.lat_bound));

        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % R1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

            for f = 1:length(par.(type).fw); fw = par.(type).fw{f};
                r1z.(fw) = flux_z.(land).res.(fw) ./ flux_z.(land).ra.(fw);
                r1z.(fw) = interp1(grid.dim2.lat, r1z.(fw), lat);
                r1z_lat.(fw) = nansum(clat_mon.*r1z.(fw))/nansum(clat);

                r1z2.(fw) = flux_z2.(land).res.(fw) ./ flux_z2.(land).ra.(fw);
                r1z2.(fw) = interp1(grid2.dim2.lat, r1z2.(fw), lat);
                r1z_lat2.(fw) = nansum(clat_mon.*r1z2.(fw))/nansum(clat);
            end

            f_vec = assign_fw(type, par);

            par.folder = sprintf('%s/dr1/mse_old/%s/0_poleward_of_lat_%g', plotdir, land, par.lat_bound);
            if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

            figure(); clf; hold all; box on;
            
            ylim_lo = 0.7;
            ylim_up = 1.7;
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
            control=plot([1:12], circshift(r1z_lat.mse_old,par.shiftby,2), 'k');
            flat=plot([1:12], circshift(r1z_lat2.mse_old,par.shiftby,2), '--k');
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', 'CMIP5 RCP8.5$-$historical', par.lat_bound, par.lat_pole));
            ylabel(sprintf('$R_1$ (unitless)'));
            legend([control flat], 'Control AA', 'Flat AA', 'location', 'northeast', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_r1z_rcp85', par.folder);
            print(plotname, '-dpng', '-r300');
            close;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % Ra
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

            ra = interp1(grid.dim2.lat, flux_z.(land).ra.mse_old, lat);
            ra_lat = nansum(clat_mon.*ra)/nansum(clat);

            ra2 = interp1(grid.dim2.lat, flux_z2.(land).ra.mse_old, lat);
            ra_lat2 = nansum(clat_mon.*ra2)/nansum(clat);

            figure(); clf; hold all; box on;
            
            ylim_lo = -180;
            ylim_up = 50;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            control=plot([1:12], circshift(ra_lat,par.shiftby,2), 'color', par.gray);
            flat=plot([1:12], circshift(ra_lat2,par.shiftby,2), '--', 'color', par.gray);
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', 'CMIP5 RCP8.5$-$historical', par.lat_bound, par.lat_pole));
            ylabel(sprintf('$R_a$ (W m$^{-2}$)'));
            legend([control flat], 'Control AA', 'Flat AA', 'location', 'northeast', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_ra_rcp85', par.folder);
            print(plotname, '-dpng', '-r300');
            close;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % MSE flux div + storage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

            res = interp1(grid.dim2.lat, flux_z.(land).res.mse_old, lat);
            res_lat = nansum(clat_mon.*res)/nansum(clat);

            res2 = interp1(grid.dim2.lat, flux_z2.(land).res.mse_old, lat);
            res_lat2 = nansum(clat_mon.*res2)/nansum(clat);

            figure(); clf; hold all; box on;
            
            ylim_lo = -180;
            ylim_up = 50;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            control=plot([1:12], circshift(res_lat,par.shiftby,2), 'color', par.maroon);
            flat=plot([1:12], circshift(res_lat2,par.shiftby,2), '--', 'color', par.maroon);
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', 'CMIP5 RCP8.5$-$historical', par.lat_bound, par.lat_pole));
            ylabel(sprintf('$\\partial_t m + \\partial_y (vm)$ (W m$^{-2}$)'));
            legend([control flat], 'Control AA', 'Flat AA', 'location', 'northeast', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_res_rcp85', par.folder);
            print(plotname, '-dpng', '-r300');
            close;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % DELTA ALL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            
            lh = interp1(grid.dim2.lat, flux_z.(land).hfls, lat);
            lh_lat = nansum(clat_mon.*lh)/nansum(clat);

            lh2 = interp1(grid.dim2.lat, flux_z2.(land).hfls, lat);
            lh_lat2 = nansum(clat_mon.*lh2)/nansum(clat);

            sh = interp1(grid.dim2.lat, flux_z.(land).hfss, lat);
            sh_lat = nansum(clat_mon.*sh)/nansum(clat);

            sh2 = interp1(grid.dim2.lat, flux_z2.(land).hfss, lat);
            sh_lat2 = nansum(clat_mon.*sh2)/nansum(clat);

            figure(); clf; hold all; box on;
            
            ylim_lo = -50;
            ylim_up = 50;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            dra=plot([1:12], circshift(ra_lat2-ra_lat,par.shiftby,2), '-', 'color', par.gray);
            dres=plot([1:12], circshift(res_lat2 - res_lat,par.shiftby,2), '-', 'color', par.maroon);
            dlh=plot([1:12], circshift(lh_lat2 - lh_lat,par.shiftby,2), '-', 'color', par.blue);
            dsh=plot([1:12], circshift(sh_lat2 - sh_lat,par.shiftby,2), '-', 'color', par.orange);
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', 'CMIP5 RCP8.5$-$historical', par.lat_bound, par.lat_pole));
            ylabel(sprintf('Energy flux$_{\\mathrm{Flat-Control}}$ (W m$^{-2}$)'));
            legend([dra dres dlh dsh], '$R_a$', '$\partial_t m + \partial_y (vm)$', 'LH', 'SH', 'location', 'eastoutside', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_dflux_rcp85', par.folder);
            print(plotname, '-dpng', '-r300');
            close;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % DECOMPOSITION
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

            dr1 = r1z_lat2.mse_old - r1z_lat.mse_old;
            dra = ra_lat2 - ra_lat;
            dres = res_lat2 - res_lat;
            rad = -r1z_lat.mse_old .* (dra./ra_lat);
            dyn = r1z_lat.mse_old .* (dres./res_lat);
            resi = dr1 - rad - dyn;

            figure(); clf; hold all; box on;
            
            ylim_lo = -0.5;
            ylim_up = 0.5;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            total=plot([1:12], circshift(dr1,par.shiftby,2), '-k');
            radcomp=plot([1:12], circshift(rad,par.shiftby,2), '-', 'color', par.gray);
            dyncomp=plot([1:12], circshift(dyn,par.shiftby,2), '-', 'color', par.maroon);
            residual=plot([1:12], circshift(resi,par.shiftby,2), '-.k');
            title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', 'CMIP5 RCP8.5$-$historical', par.lat_bound, par.lat_pole));
            ylabel(sprintf('$R_{1,\\mathrm{ flat}}-R_{1,\\mathrm{ control}}$ (unitless)'));
            legend([total radcomp dyncomp residual], '$\Delta R_1$', '$-R_{1,\mathrm{ control}}\frac{\Delta R_a}{R_{a,\mathrm{ control}}}$', '$R_{1,\mathrm{ control}}\frac{\Delta (\partial_t m + \partial_y(vm))}{(\partial_t m + \partial_y(vm))_{\mathrm{ control}}}$', 'Residual', 'location', 'eastoutside', 'orientation', 'vertical');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            plotname = sprintf('%s/0_mon_decompr1z_rcp85', par.folder);
            print(plotname, '-dpng', '-r300');
            close;
            

        end % for land
    end % lat bound

end % for function
