function plot_sice_rea(type, par)
% plot sicedow depth
    make_dirs(type, par)

    % load data
    % [~, ~, ~, lat, par] = load_flux(type, par);

    prefix = make_prefix(type, par);
    plotdir = make_plotdir(type, par);
    load(sprintf('%s/grid.mat', prefix)); % read grid data

    rea_list = {'era5c', 'merra2c', 'jra55'};

    sice_mat = nan([length(grid.dim2.lon) length(grid.dim2.lat) 12 length(rea_list)]);
    lh_mat = nan([length(grid.dim2.lon) length(grid.dim2.lat) 12 length(rea_list)]);
    sftlf_mat = nan([length(grid.dim2.lon) length(grid.dim2.lat) length(rea_list)]);

    for tt = 1:length(rea_list); ttype = rea_list{tt};
        prefix = make_prefix(ttype, par);
        tmp = load(sprintf('%s/grid.mat', prefix)); % read grid data
        grid_tmp = tmp.grid; clear tmp;
        load(sprintf('%s/sice.mat', prefix)); % read sea ice
        load(sprintf('%s/sftlf.mat', prefix)); % read land fraction to plot coastline
        load(sprintf('%s/stf.mat', prefix)); % read lh
        lh = rename_lh(ttype, stf);

        if length(size(sftlf)) > 2
            sftlf = squeeze(nanmean(sftlf,3));
        end

        % interpolate to common lon grid
        sice_i = interp1(grid_tmp.dim2.lon, sice, grid.dim2.lon);
        lh_i = interp1(grid_tmp.dim2.lon, lh, grid.dim2.lon);
        sftlf_i = interp1(grid_tmp.dim2.lon, sftlf, grid.dim2.lon);

        sice_i = permute(sice_i, [2 1 3]);
        lh_i = permute(lh_i, [2 1 3]);
        sftlf_i = permute(sftlf_i, [2 1]);

        % interpolate to common lat grid
        sice_i = interp1(grid_tmp.dim2.lat, sice_i, grid.dim2.lat);
        lh_i = interp1(grid_tmp.dim2.lat, lh_i, grid.dim2.lat);
        sftlf_i = interp1(grid_tmp.dim2.lat, sftlf_i, grid.dim2.lat);

        sice_i = permute(sice_i, [2 1 3]);
        lh_i = permute(lh_i, [2 1 3]);
        sftlf_i = permute(sftlf_i, [2 1]);

        % save
        sice_mat(:,:,:,tt) = sice_i;
        lh_mat(:,:,:,tt) = lh_i;
        sftlf_mat(:,:,tt) = sftlf_i;

        clear sice sice_i lh lh_i sftlf sftlf_i

    end

    % take reanalysis mean
    sice = nanmean(sice_mat, 4);
    lh = nanmean(lh_mat, 4);
    sftlf = nanmean(sftlf_mat, 3);

    % annual mean
    sice_ann = nanmean(sice, 3);
    lh_ann = nanmean(lh, 3);

    folder_prefix = sprintf('%s/sice/', plotdir);
    folder_list = {'nh_hl', 'sh_hl'};
    for fs = 1:length(folder_list); folder_suffix = folder_list{fs};
        foldername = sprintf('%s%s', folder_prefix, folder_suffix);
        if ~exist(foldername, 'dir')
            mkdir(foldername);
        end
    end

    % NH
    idx_nh = grid.dim2.lat >= 45;
    idx_sh = grid.dim2.lat <= -45;

    [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim2.lon, grid.dim2.lat);
    [mesh_ll_lat_nh, mesh_ll_lon_nh] = meshgrid(grid.dim2.lon, grid.dim2.lat(idx_nh));
    [mesh_ll_lat_sh, mesh_ll_lon_sh] = meshgrid(grid.dim2.lon, grid.dim2.lat(idx_sh));

    % SEA ICE
    var_text = 'Sea ice fraction (1)';
    figure(); clf; hold all; box on;
    colormap('parula')
    contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice_ann(:,idx_nh))', 0:0.1:1, 'linecolor', 'none', 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
    % contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice_ann(:,idx_nh))', 0.5*[1 1], 'r', 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 70*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 60*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 50*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    caxis([0 1]);
    cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, var_text);
    xlabel('eastings (km)');
    ylabel('northings (km)');
    title('Reanalysis mean, ANN')
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
    set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
    print(sprintf('%s/sice/nh_hl/stereo_sice_ll_ann', plotdir), '-dpng', '-r300');

    % LH
    var_text = 'LH (W m$^{-2}$)';
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp)
    contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh_ann(:,idx_nh))', 0:2.5:50, 'linecolor', 'none', 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
    % contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 70*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 60*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 50*[1 1], ':', 'linecolor', 0.5*[1 1 1], 'km');
    % contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh(:,idx_nh,mon))', 273.15*[1 1], 'r', 'km');
    caxis([-50 50]);
    cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, var_text);
    xlabel('eastings (km)');
    ylabel('northings (km)');
    title('Reanalysis mean, ANN')
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
    set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
    print(sprintf('%s/sice/nh_hl/stereo_lh_ll_ann', plotdir), '-dpng', '-r300');

    return


    for mon = 1:12
        mon_str = make_mon_str(mon);
        [xmin, xmax] = make_lon_bounds(type);

        %%%%%%%% STEREOGRAPHIC PROJECTION %%%%%%%%%%%

        %%%%%% NORTH %%%%%%%%%

        % SEA ICE
        var_text = 'Sea ice fraction (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice(:,idx_nh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice(:,idx_nh,mon))', 0.5*[1 1], 'r', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_sice_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % LH
        var_text = 'LH (W m$^{-2}$)';
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp)
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh(:,idx_nh,mon))', 0:2.5:50, 'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        % contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh(:,idx_nh,mon))', 273.15*[1 1], 'r', 'km');
        caxis([-50 50]);
        cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_lh_ll_%02d', plotdir, mon), '-dpng', '-r300');


        %%%%%% SOUTH %%%%%%%%%

        % SEA ICE
        var_text = 'Sea ice fraction (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(sice(:,idx_sh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.5*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(sice(:,idx_sh,mon))', 0.5*[1 1], 'r', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_sice_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % LH
        var_text = 'LH (W m$^{-2}$)';
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp)
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(lh(:,idx_sh,mon))', 0:2.5:50, 'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.5*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        % contourpsn(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(lh(:,idx_sh,mon))', 273.15*[1 1], 'r', 'km');
        caxis([-50 50]);
        cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_lh_ll_%02d', plotdir, mon), '-dpng', '-r300');


    end

end
