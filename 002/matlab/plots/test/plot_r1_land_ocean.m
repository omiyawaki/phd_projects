function plot_r1_land_ocean(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    [~, ~, ~, lat, ~] = load_flux(type, par);

    lat_mid = linspace(40,60,50);
    if sh; lat_mid = -lat_mid; end;

    clat = cosd(lat_mid);

    % load data
    prefix=make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/sftlf.mat', prefix));
    load(sprintf('%s/flux_z.mat', prefix_proc));
    % load(sprintf('%s/si_bl_0.9/ga_malr_bl_diff_poleward_of_lat_80.mat', prefix_proc));

    % zonal mean land fraction
    sftlf = squeeze(nanmean(sftlf, 1));

    % r1 land and ocean
    r1 = flux_z.lo.res.mse_old./flux_z.lo.ra.mse_old;
    r1_land = flux_z.l.res.mse_old./flux_z.l.ra.mse_old;
    r1_ocean = flux_z.o.res.mse_old./flux_z.o.ra.mse_old;

    % area average over AA
    sftlf_i = interp1(grid.dim2.lat, sftlf, lat_mid);
    sftlf_mid = nansum(clat'.*sftlf_i) / nansum(clat);
    r1_i = interp1(grid.dim2.lat, r1, lat_mid);
    r1_mid = nansum(clat'.*r1_i) / nansum(clat);
    r1_land_i = interp1(grid.dim2.lat, r1_land, lat_mid);
    r1_land_mid = nansum(clat'.*r1_land_i) / nansum(clat);
    r1_ocean_i = interp1(grid.dim2.lat, r1_ocean, lat_mid);
    r1_ocean_mid = nansum(clat'.*r1_ocean_i) / nansum(clat);

    % scale land and ocean r1 by area fraction
    r1_land_mid = sftlf_mid.*r1_land_mid;
    r1_ocean_mid = (1-sftlf_mid).*r1_ocean_mid;

    % compute residual
    r1_res_mid = r1_mid - r1_land_mid - r1_ocean_mid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of r1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_res=plot(1:12, circshift(r1_res_mid,shiftby), '-.k');
    h=plot(1:12, circshift(r1_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(r1_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(r1_ocean_mid,shiftby), '-', 'color', par.blue);
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_mon_mid', plotdir), '-dpng', '-r300');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % leg only
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % axis off;
    % axis([100 101 100 101])
    % title('');
    % legend(flip(leg), '5 m', '15 m', '25 m', '40 m', '50 m', 'location', 'northwest', 'numcolumns', 5, 'orientation', 'horizontal');
    % set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 6.2 0.5])
    % print(sprintf('%s/alb/mld_legend', plotdir), '-dpng', '-r300');

    close;

end
