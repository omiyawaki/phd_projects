function plot_alb_comp(type, par)
% plot albedoow depth
    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    if ~strcmp(type, 'echam')
        error('This analysis only works for ECHAM')
    end

    % load data
    [~, ~, ~, lat, ~] = load_flux(type, par);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    prefix_proc = make_prefix_proc(type, par);
    tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    tmp=load(sprintf('%s/albedo.mat', prefix)); albedo=tmp.albedo; clear tmp; % read clear sky albedoedo data
    tmp=load(sprintf('%s/srfc.mat', prefix)); tas=tmp.srfc.temp2; clear tmp; % read clear sky albedoedo data
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data

    % load additional data
    type2='era5c';
    prefix2=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type2, par.(type2).yr_span);
    prefix_proc2 = make_prefix_proc(type2, par);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);

    tmp=load(sprintf('%s/grid.mat', prefix2)); grid2=tmp.grid; clear tmp; % read grid data
    tmp=load(sprintf('%s/albedo.mat', prefix2)); albedo2=tmp.albedo; clear tmp; % read clear sky albedoedo data
    tmp=load(sprintf('%s/srfc.mat', prefix2)); tas2=tmp.srfc.t2m; clear tmp; % read clear sky albedoedo data
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc2)); flux2=tmp.flux_z; clear tmp; % read clear sky albedoedo data
    
    lh = -flux.lo.ahfl;
    lh2 = -flux2.lo.slhf;

    r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
    r12 = flux2.lo.res.mse_old./flux2.lo.ra.mse_old;
    
    albedo = squeeze(nanmean(albedo,1)); % zonal mean
    albedo2 = squeeze(nanmean(albedo2,1)); % zonal mean
    tas = squeeze(nanmean(tas,1)); % zonal mean
    tas2 = squeeze(nanmean(tas2,1)); % zonal mean

    lat_aa = linspace(-80,-90,50);
    clat = cosd(lat_aa);
    
    % area average over AA
    albedo_i = interp1(grid.dim2.lat, albedo, lat_aa);
    albedo_aa = nansum(clat'.*albedo_i) / nansum(clat);
    r1_i = interp1(grid.dim2.lat, r1, lat_aa);
    r1_aa = nansum(clat'.*r1_i) / nansum(clat);
    lh_i = interp1(grid.dim2.lat, lh, lat_aa);
    lh_aa = nansum(clat'.*lh_i) / nansum(clat);
    tas_i = interp1(grid.dim2.lat, tas, lat_aa);
    tas_aa = nansum(clat'.*tas_i) / nansum(clat);

    albedo2_i = interp1(grid2.dim2.lat, albedo2, lat_aa);
    albedo2_aa = nansum(clat'.*albedo2_i) / nansum(clat);
    r12_i = interp1(grid2.dim2.lat, r12, lat_aa);
    r12_aa = nansum(clat'.*r12_i) / nansum(clat);
    lh2_i = interp1(grid2.dim2.lat, lh2, lat_aa);
    lh2_aa = nansum(clat'.*lh2_i) / nansum(clat);
    tas2_i = interp1(grid2.dim2.lat, tas2, lat_aa);
    tas2_aa = nansum(clat'.*tas2_i) / nansum(clat);

    % mon x lat of albedoow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    l1=plot(1:12, circshift(albedo_aa,6), '-k');
    l2=plot(1:12, circshift(albedo2_aa,6), '--k');
    ylabel('Surface albedo (unitless)');
    make_title_type(type,par)
    legend([l1, l2], 'ECHAM', 'ERA5', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelsh, 'ylim', [0 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/albedo_mon_%s', plotdir, type2), '-dpng', '-r300');
    close;

    % mon x lat of albedoow depth
    figure(); clf; hold all; box on;
    raemin=0.9;
    ylim_up =2.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    l1=plot(1:12, circshift(r1_aa,6), '-k');
    l2=plot(1:12, circshift(r12_aa,6), '--k');
    ylabel('$R_1$ (unitless)');
    make_title_type(type,par)
    legend([l1, l2], 'ECHAM', 'ERA5', 'location', 'northeast')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelsh, 'ylim', [0.5 ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/r1_mon_%s', plotdir, type2), '-dpng', '-r300');
    close;

    % mon x lat of lhow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    line([1 12], [0 0], 'linewidth', 0.5)
    l1=plot(1:12, circshift(lh_aa,6), '-', 'color', par.blue);
    l2=plot(1:12, circshift(lh2_aa,6), '--', 'color', par.blue);
    ylabel('LH (W m$^{-2}$)');
    make_title_type(type,par)
    legend([l1, l2], 'ECHAM', 'ERA5', 'location', 'northwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelsh, 'ylim', [-1 10], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/lh_mon_%s', plotdir, type2), '-dpng', '-r300');
    close;

    % mon x lat of tasow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    l1=plot(1:12, circshift(tas_aa,6), '-k');
    l2=plot(1:12, circshift(tas2_aa,6), '--k');
    ylabel('2 m temperature (K)');
    make_title_type(type,par)
    legend([l1, l2], 'ECHAM', 'ERA5', 'location', 'southeast')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelsh, 'ylim', [180 260], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/tas_mon_%s', plotdir, type2), '-dpng', '-r300');
    close;

end
