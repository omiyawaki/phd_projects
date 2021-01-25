function plot_olr_ts(type, par)
    make_dirs(type, par);

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/radcs.mat', prefix)); % read clear sky radiation data
    load(sprintf('%s/srfc.mat', prefix)); % read surface data
    load(sprintf('%s/flux.mat', prefix_proc)); % read processed flux data
    
    % load tas
    tas = load_tas(srfc, type, par);
    clear srfc;
    
    % load clear sky OLR
    olrcs = load_olrcs(radcs, type, par);
    clear radcs;
    
    % OLRcs vs Ts regression
    X = [ones(length(tas(:)),1) tas(:)];
    beta = X\(-olrcs(:));
    % OLRas vs Ts regression
    beta_as = X\(-flux.olr(:));
    
    % regression line
    tvec = linspace(200,320,100);
    olrreg = beta(1) + beta(2)*tvec;
    Rsq = 1 - sum((-olrcs(:) - (beta(1)+beta(2)*tas(:))).^2)/sum((-olrcs(:) - mean(-olrcs(:))).^2);
    olrreg_as = beta_as(1) + beta_as(2)*tvec;
    Rsq_as = 1 - sum((-flux.olr(:) - (beta_as(1)+beta_as(2)*tas(:))).^2)/sum((-flux.olr(:) - mean(-flux.olr(:))).^2);
    
    % repeat for zonal mean data
    tas_z = squeeze(nanmean(tas,1));
    olrcs_z = squeeze(nanmean(olrcs,1));
    X_z = [ones(length(tas_z(:)),1) tas_z(:)];
    beta_z = X_z\(-olrcs_z(:));
    olrreg_z = beta_z(1) + beta_z(2)*tvec;
    Rsq_z = 1 - sum((-olrcs_z(:) - (beta_z(1)+beta_z(2)*tas_z(:))).^2)/sum((-olrcs_z(:) - mean(-olrcs_z(:))).^2);
    olras_z = squeeze(nanmean(flux.olr,1));
    beta_as_z = X_z\(-olras_z(:));
    olrreg_as_z = beta_as_z(1) + beta_as_z(2)*tvec;
    Rsq_as_z = 1 - sum((-olras_z(:) - (beta_as_z(1)+beta_as_z(2)*tas_z(:))).^2)/sum((-olras_z(:) - mean(-olras_z(:))).^2);
    
    figure(); clf; hold all; box on;
    plot(tas(:), -flux.olr(:), '.k', 'markersize', 1);
    plot(tvec, olrreg_as, '-r');
    text(210,270,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta_as(1), beta_as(2)));
    text(280,125,sprintf('$R^2 = %.2f$', Rsq_as));
    xlabel('$T_\mathrm{2\, m}$');
    ylabel('$\mathrm{OLR_{all}}$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    %set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
    print(sprintf('%s/olr_ts/all_sky_tas.png', plotdir), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    plot(tas(:), -olrcs(:), '.k', 'markersize', 1);
    plot(tvec, olrreg, '-r');
    text(210,300,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta(1), beta(2)));
    text(280,125,sprintf('$R^2 = %.2f$', Rsq));
    xlabel('$T_\mathrm{2\, m}$');
    ylabel('$\mathrm{OLR_{clear}}$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    %set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
    print(sprintf('%s/olr_ts/clear_sky_tas.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(tas_z(:), -olrcs_z(:), '.k', 'markersize', 1);
    plot(tvec, olrreg_z, '-r');
    text(210,300,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta_z(1), beta_z(2)));
    text(280,125,sprintf('$R^2 = %.2f$', Rsq_z));
    xlabel('$T_\mathrm{2\, m}$');
    ylabel('$\mathrm{OLR_{clear}}$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    %set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
    print(sprintf('%s/olr_ts/clear_sky_tas_zonmean.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(tas_z(:), -olras_z(:), '.k', 'markersize', 1);
    plot(tvec, olrreg_as_z, '-r');
    text(210,270,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta_as_z(1), beta_as_z(2)));
    text(280,125,sprintf('$R^2 = %.2f$', Rsq_as_z));
    xlabel('$T_\mathrm{2\, m}$');
    ylabel('$\mathrm{OLR_{all}}$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    %set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
    print(sprintf('%s/olr_ts/all_sky_tas_zonmean.png', plotdir), '-dpng', '-r300');
    close;
    
    % extratropics only
    tas_z_et = tas_z(grid.dim2.lat<=-25 | grid.dim2.lat>=25, :);
    olrcs_z_et = olrcs_z(grid.dim2.lat<=-25 | grid.dim2.lat>=25, :);
    X_z_et = [ones(length(tas_z_et(:)),1) tas_z_et(:)];
    beta_z_et = X_z_et\(-olrcs_z_et(:));
    olrreg_z_et = beta_z_et(1) + beta_z_et(2)*tvec;
    Rsq_z_et = 1 - sum((-olrcs_z_et(:) - (beta_z_et(1)+beta_z_et(2)*tas_z_et(:))).^2)/sum((-olrcs_z_et(:) - mean(-olrcs_z_et(:))).^2);

    figure(); clf; hold all; box on;
    plot(tas_z_et(:), -olrcs_z_et(:), '.k', 'markersize', 1);
    plot(tvec, olrreg_z_et, '-r');
    text(210,300,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta_z_et(1), beta_z_et(2)));
    text(280,125,sprintf('$R^2 = %.2f$', Rsq_z_et));
    xlabel('$T_\mathrm{2\, m}$');
    ylabel('$\mathrm{OLR_{clear}}$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    %set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
    print(sprintf('%s/olr_ts/clear_sky_tas_zonmean_et.png', plotdir), '-dpng', '-r300');
    close;
    
end
