function plot_r1_mid_mld(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

% plot albedoow depth
    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    if ~strcmp(type, 'echam')
        error('This analysis only works for ECHAM')
    end

    par.echam_clims = {"rp000046",... % 50 m
                       "rp000135",... % 40 m
                        "rp000145",... % 25 m
                       "rp000141",... % 15 m
                        "rp000086",... % 5 m
                        "rp000172"}; % 3 m

    [~, ~, ~, lat, ~] = load_flux(type, par);

    lat_mid = linspace(40,60,50);
    if sh; lat_mid = -lat_mid; end;

    clat = cosd(lat_mid);
        
    r1_mid={};
    n_clims = length(par.echam_clims);
    for i = 1:n_clims; clim=par.echam_clims{i};
        par.echam.clim = clim;
        % load data
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, clim);
        prefix_proc = make_prefix_proc(type, par);
        tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
        tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data

        r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
    
        % area average over AA
        r1_i = interp1(grid.dim2.lat, r1, lat_mid);
        r1_mid.(clim) = nansum(clat'.*r1_i) / nansum(clat);
        
        clear grid flux r1 r1_i

    end                                    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r1 seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-1.1;
    ylim_up =0.8;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    for i = 1:n_clims; clim=par.echam_clims{i};
        % leg(i)=plot(1:12, circshift(r1_mid.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
        leg(i)=plot(1:12, circshift(r1_mid.(clim),shiftby), '-');
    end
    ylabel('$R_1$ (unitless)');
    title('AQUA w/o ice, Midlatitudes')
    % legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'northeast')
    leg=legend(flip(leg), '3 m', '5 m', '15 m', '25 m', '40 m', '50 m', 'location', 'southeast')
    leg.ItemTokenSize = [10,10];
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/mid_mld/r1_mon_mld', plotdir), '-dpng', '-r300');
    if par.make_tikz
        matlab2tikz(sprintf('%s/mid_mld/r1_mon_mld.tex', plotdir));
    end
    close;

end
