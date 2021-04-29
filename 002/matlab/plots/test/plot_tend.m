function plot_tend(type, par) 

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    path_don = '/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology';
    
    % load grid
    tmp = load(sprintf('%s/grid.mat', prefix)); grid = tmp.grid; clear tmp;
    tmp = load(path_don); grid_don.lat = tmp.lat; grid_don.lon = tmp.lon; clear tmp;
    
    % load mse tendency data
    tmp = load(sprintf('%s/tend.mat', prefix));
    tend = tmp.tend.tend;
    tendmon = tmp.tend.tendmon;
    % cT = tmp.tend.int;
    % qL = tmp.tend.lat;
    % gz = tmp.tend.pot;
    clear tmp;
    % load donohoe tendency
    tmp = load(path_don); tend_don = tmp.TETEN; clear tmp;
    % reshape to (lon x lat x lev)
    tend_don = permute(tend_don, [3 2 1]);
    
    % take zonal mean
    tend = squeeze(nanmean(tend,1));
    tendmon = squeeze(nanmean(tendmon,1));
    % cT = squeeze(nanmean(cT,1));
    % qL = squeeze(nanmean(qL,1));
    % gz = squeeze(nanmean(gz,1));
    tend_don = squeeze(nanmean(tend_don,1));
    
    % plot mon x lat structure of tendency
    [mesh_lat, mesh_mon] = meshgrid(1:12, grid.dim3.lat);
    [mesh_don_lat, mesh_don_mon] = meshgrid(1:12, grid_don.lat);
    
    plot_tend_mon_lat(tend, mesh_lat, mesh_mon, '', plotdir, type, par);
    plot_tend_mon_lat(tend_don, mesh_don_lat, mesh_don_mon, '_don', plotdir, type, par);
    
    % plot the difference
    tend_doni = interp1(grid_don.lat, tend_don, grid.dim3.lat);
    tend_diff = tend - tend_doni;
    
    % plot_tend_mon_lat(tend_diff, mesh_lat, mesh_mon, '_diff', plotdir, type, par);
    % plot_tend_lat(type, tend, cT, qL, gz, grid, plotdir, par);
    plot_tend_lat(type, tend, tendmon, grid, plotdir, par);

end

function plot_tend_lat(type, tend, tendmon, grid, plotdir, par)

    %% TEMPORAL COMPARISON
    lat_list = [10 -10 80 -80];
    lat_center = 50;
    
    for lb = 1:length(lat_list); par.lat_bound = lat_list(lb);

        if par.lat_bound > 0
            hemi = 'nh';
            shiftby = 0;
            monlabel = par.monlabelnh;
        else
            hemi = 'sh';
            shiftby = 6;
            monlabel = par.monlabelsh;
        end

        if abs(par.lat_bound) > 45  % polar
            lat_print = 'polar';
            [lat, clat, clat_mon, par] = make_polar_lat(par);
            lat1 = par.lat_bound;
            lat2 = par.lat_pole;
        else % midlats
            lat_print = 'mid';
            [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);
            lat1 = par.lat_center-par.lat_bound;
            lat2 = par.lat_center+par.lat_bound;
        end

        % interpolate at select latitude
        tend_lat = interp1(grid.dim2.lat, tend, lat); % multiply by 100 to express in percentage
        tend_lat = nansum(clat_mon .* tend_lat)/nansum(clat); % area averaged weighting
        tendmon_lat = interp1(grid.dim2.lat, tendmon, lat); % multiply by 100 to express in percentage
        tendmon_lat = nansum(clat_mon .* tendmon_lat)/nansum(clat); % area averaged weighting
        % cT_lat = interp1(grid.dim2.lat, cT, lat); % multiply by 100 to express in percentage
        % cT_lat = nansum(clat_mon .* cT_lat)/nansum(clat); % area averaged weighting
        % qL_lat = interp1(grid.dim2.lat, qL, lat); % multiply by 100 to express in percentage
        % qL_lat = nansum(clat_mon .* qL_lat)/nansum(clat); % area averaged weighting
        % gz_lat = interp1(grid.dim2.lat, gz, lat); % multiply by 100 to express in percentage
        % gz_lat = nansum(clat_mon .* gz_lat)/nansum(clat); % area averaged weighting

        figure(); clf; hold all; box on;

        plot(1:12, circshift(tend_lat, shiftby), '-', 'color', 'k');
        plot(1:12, circshift(tendmon_lat, shiftby), '--', 'color', 'k');
        % plot(1:12, circshift(cT_lat, shiftby), '-', 'color', par.orange);
        % plot(1:12, circshift(qL_lat, shiftby), '-', 'color', par.blue);
        % plot(1:12, circshift(gz_lat, shiftby), '-', 'color', par.green);
        ylabel(sprintf('Energy Flux (W m$^{-2}$)'));
        % legend('$\partial_t h$', '$\partial_t (c_p T)$', '$\partial_t (qL)$', '$\partial_t (gz)$');
        legend('Daily', 'Monthly');
        make_title_type_lat(type, lat1, lat2, par);
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'yminortick', 'on', 'ylim', [-50 50]) 
        print(sprintf('%s/tend/tend_decomp_%s_%s.png', plotdir, hemi, lat_print), '-dpng', '-r300');
        close;

    end % lat bounds

end

function plot_tend_mon_lat(tend, mesh_lat, mesh_mon, label, plotdir, type, par)

    % lat x mon dependence of mse tendency
    figure(); clf; hold all; box on;
    cmp = colCog(20);
    colormap(cmp);
    
    if label == "_diff"
        contourf(mesh_lat, mesh_mon, tend, [-20 -10:1:10 20], 'linecolor', 'w');
    else
        contourf(mesh_lat, mesh_mon, tend, [-100 -50:5:50 100], 'linecolor', 'w');
    end
    %contour(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
    
    if label == ""
        make_title_type(type, par);
    elseif label == "_don"
        title('ERA-I (Donohoe)');
    elseif label == "_diff"
        if strcmp(type, 'era5c');
            title(sprintf('%s -- ERA-I (Donohoe)', 'ERA5'));
        else
            title(sprintf('%s -- ERA-I (Donohoe)', upper(type)));
        end
    end
    
    if label == "_diff"
        var_text = '$\Delta \partial_t h$';
        caxis([-20 20]);
        cb = colorbar('limits', [-20 20], 'ytick', [-20:10:20], 'location', 'eastoutside');
    else
        var_text = '$\partial_t h$';
        caxis([-50 50]);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
    end
    
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
    xlabel('Month'); ylabel('Latitude (deg)');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    foldername = sprintf('%s/tend', plotdir);
    if ~exist(foldername, 'dir'); mkdir(foldername); end;
    print(sprintf('%s/0_tend_mon_lat%s', foldername, label), '-dpng', '-r300');
    close;

end
