function plot_tend(type, par) 

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    path_don = '/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology';
    
    % load grid
    tmp = load(sprintf('%s/grid.mat', prefix)); grid = tmp.grid; clear tmp;
    tmp = load(path_don); grid_don.lat = tmp.lat; grid_don.lon = tmp.lon; clear tmp;
    
    % load main mse tendency data
    tmp = load(sprintf('%s/tend.mat', prefix)); tend = tmp.tend.tend; clear tmp;
    % load donohoe tendency
    tmp = load(path_don); tend_don = tmp.TETEN; clear tmp;
    % reshape to (lon x lat x lev)
    tend_don = permute(tend_don, [3 2 1]);
    
    % take zonal mean
    tend = squeeze(nanmean(tend,1));
    tend_don = squeeze(nanmean(tend_don,1));
    
    % plot mon x lat structure of tendency
    [mesh_lat, mesh_mon] = meshgrid(1:12, grid.dim3.lat);
    [mesh_don_lat, mesh_don_mon] = meshgrid(1:12, grid_don.lat);
    
    plot_tend_mon_lat(tend, mesh_lat, mesh_mon, '', plotdir, type, par);
    plot_tend_mon_lat(tend_don, mesh_don_lat, mesh_don_mon, '_don', plotdir, type, par);
    
    % plot the difference
    tend_doni = interp1(grid_don.lat, tend_don, grid.dim3.lat);
    tend_diff = tend - tend_doni;
    
    plot_tend_mon_lat(tend_diff, mesh_lat, mesh_mon, '_diff', plotdir, type, par);

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
