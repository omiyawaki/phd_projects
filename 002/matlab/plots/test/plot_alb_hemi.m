function plot_alb_hemi(type, par)
% plot albedoow depth
    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    % load data
    [~, ~, ~, lat, ~] = load_flux(type, par);
    prefix=make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    tmp=load(sprintf('%s/albedo.mat', prefix)); albedo=tmp.albedo; clear tmp; % read clear sky albedoedo data
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data

    r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
    
    albedo = squeeze(nanmean(albedo,1)); % zonal mean

    lat_nh = linspace(80,90,50);
    clat_nh = cosd(lat_nh);
    lat_sh = linspace(-80,-90,50);
    clat_sh = cosd(lat_sh);
    
    % area average over NH
    albedo_i = interp1(grid.dim2.lat, albedo, lat_nh);
    albedo_nh = nansum(clat_nh'.*albedo_i) / nansum(clat_nh);
    r1_i = interp1(grid.dim2.lat, r1, lat_nh);
    r1_nh = nansum(clat_nh'.*r1_i) / nansum(clat_nh);

    % area average over SH
    albedo_i = interp1(grid.dim2.lat, albedo, lat_sh);
    albedo_sh = nansum(clat_sh'.*albedo_i) / nansum(clat_sh);
    r1_i = interp1(grid.dim2.lat, r1, lat_sh);
    r1_sh = nansum(clat_sh'.*r1_i) / nansum(clat_sh);

    % mon x lat of albedoow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    lnh=plot(1:12, circshift(albedo_nh,0), '-k');
    lsh=plot(1:12, circshift(albedo_sh,6), '--k');
    if strcmp(type, 'hahn');
        ylabel('$\mathrm{SW_{sfc, up}/SW_{sfc, down}}$ (unitless)');
    else
        ylabel('Surface albedo (unitless)');
    end
    make_title_type(type,par)
    legend([lnh, lsh], 'NH', 'SH (shifted 6 mo.)', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [0.5 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/albedo_mon_hemi', plotdir), '-dpng', '-r300');
    close;
    
    % mon x lat of albedoow depth
    figure(); clf; hold all; box on;
    raemin=0.9;
    ylim_up =2.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    lnh=plot(1:12, circshift(r1_nh,0), '-k');
    lsh=plot(1:12, circshift(r1_sh,6), ':k');
    ylabel('$R_1$ (unitless)');
    make_title_type(type,par)
    legend([lnh, lsh], 'NH', 'SH (shifted 6 mo.)', 'location', 'northeast')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [0.5 ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/r1_mon_hemi', plotdir), '-dpng', '-r300');
    close;

end
