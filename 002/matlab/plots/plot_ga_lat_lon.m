function plot_ga_lat_lon(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    % load GCM data
    prefix=make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    load(sprintf('%s/grid.mat', prefix));

    tmp=load(sprintf('%s/si_bl_%g/ga_malr_diff_t_%g.mat', prefix_proc, par.si_bl, par.si_up)); ga=tmp.ga_malr_diff_t.lo; clear tmp; % read clear sky albedoedo data

    if strcmp(type, 'rea')
        parera = par;
        parera.lat_interp = 'native';
        prefixera=make_prefix('era5c', parera);
        temp=load(sprintf('%s/grid.mat', prefixera)); grid_era = temp.grid; clear temp;
        temp=load(sprintf('%s/sftlf.mat', prefixera)); sftlf_era = temp.sftlf; clear temp;

        sftlf_era = nanmean(sftlf_era, 3);

    end

    [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim3.lon, grid.dim3.lat);
    [mesh_ll_lat_era, mesh_ll_lon_era] = meshgrid(grid_era.dim3.lon, grid_era.dim3.lat);

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % lon x lat of ga
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % lat x lon of RCE and RAE
        var_text = '$[(\Gamma - \Gamma_m)/\Gamma_m]_{0.7}^{0.3}$ (\%)';
        figure(); clf; hold all; box on;
        cmp = colCog(20);
        colormap(cmp);
        contourf(mesh_ll_lat, mesh_ll_lon, ga.(time)', [-50:5:50], 'linecolor', 'none');
        contour(mesh_ll_lat, mesh_ll_lon, ga.(time)', [0 0], 'linecolor', 0.5*[1 1 1]);
        caxis([-50 50]);
        cb = colorbar('limits', [-50 50], 'ytick', [-50:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        contour(mesh_ll_lat_era, mesh_ll_lon_era, sftlf_era', 0.5*[1 1], 'k');
        make_title_type_time(type, time, par);
        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_lo/ga_lon_lat_%s', plotdir, time), '-dpng', '-r300');
        close;
    end


end
