function plot_alb_mld(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

% plot albedoow depth
    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    if ~strcmp(type, 'echam')
        error('This analysis only works for ECHAM')
    end

    % configurations to loop over
    par.echam_clims = {            "rp000126",... % 50 m
                                     "rp000148",... % 45 m
                                     "rp000134",... % 40 m
                                     "rp000146",... % 35 m
                                     "rp000130",... % 30 m
                                     "rp000144"}; % 25 m
                             
    [~, ~, ~, lat, ~] = load_flux(type, par);

    lat_aa = linspace(80,90,50);
    if sh; lat_aa = -lat_aa; end;

    clat = cosd(lat_aa);
        
    albedo_aa={}; r1_aa={};
    n_clims = length(par.echam_clims);
    for i = 1:n_clims; clim=par.echam_clims{i};
        par.echam.clim = clim;
        % load data
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, clim);
        prefix_proc = make_prefix_proc(type, par);
        tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
        tmp=load(sprintf('%s/albedo.mat', prefix)); albedo=tmp.albedo; clear tmp; % read clear sky albedoedo data
        tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data

        r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
        albedo = squeeze(nanmean(albedo,1)); % zonal mean
    
        % area average over AA
        albedo_i = interp1(grid.dim2.lat, albedo, lat_aa);
        albedo_aa.(clim) = nansum(clat'.*albedo_i) / nansum(clat);
        r1_i = interp1(grid.dim2.lat, r1, lat_aa);
        r1_aa.(clim) = nansum(clat'.*r1_i) / nansum(clat);
        
        clear grid albedo flux r1 albedo_i r1_i

    end                                    

    % mon x lat of albedoow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims; clim=par.echam_clims{i};
        leg(i)=plot(1:12, circshift(albedo_aa.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
    end
    ylabel('Surface albedo (unitless)');
    title('ECHAM w/ ice')
    legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/albedo_mon_icemld', plotdir), '-dpng', '-r300');
    close;
    
    % mon x lat of r1
    figure(); clf; hold all; box on;
    raemin=0.9;
    ylim_up =2.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    for i = 1:n_clims; clim=par.echam_clims{i};
        leg(i)=plot(1:12, circshift(r1_aa.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
    end
    ylabel('$R_1$ (unitless)');
    title('ECHAM w/ ice')
    legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'northeast')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0.5 ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/r1_mon_icemld', plotdir), '-dpng', '-r300');
    close;

end
