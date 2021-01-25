function plot_divfm_lapts(type, par)
    make_dirs(type, par);

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % read surface data
    load(sprintf('%s/flux.mat', prefix_proc)); % read processed flux data
    
    % load tas
    tas = load_tas(srfc, type, par);
    % take zonal average
    tas = squeeze(nanmean(tas));
    clear srfc;
    
    % load mse flux divergence
    divfm = flux.res.mse;
    divfm = squeeze(nanmean(divfm));
    clear flux;
    
    % compute 1/cos(phi) * d/dphi ( cos(phi) * dTs/dphi )
    clat = repmat(cosd(grid.dim2.lat), [1 12]);
    dphi = deg2rad(grid.dim2.lat(2)-grid.dim2.lat(1)); % uniform lat grid
    dtdlat = nan(size(tas));
    dtdlat(2:end-1,:) = (tas(3:end,:) - tas(1:end-2,:))/(2*dphi); % central difference
    dtdlat(1,:) = (tas(2,:) - tas(1,:))/dphi; % forward diff
    dtdlat(end,:) = (tas(end,:) - tas(end-1,:))/dphi; % backward diff
    dtdlat = clat .* dtdlat;
    lapt = nan(size(tas));
    lapt(2:end-1,:) = (dtdlat(3:end,:) - dtdlat(1:end-2,:))/(2*dphi); % central difference
    lapt(1,:) = (dtdlat(2,:) - dtdlat(1,:))/dphi; % forward diff
    lapt(end,:) = (dtdlat(end,:) - dtdlat(end-1,:))/dphi; % backward diff
    lapt = lapt./clat;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(tas,2), '-k');
    xlabel('latitude (deg)');
    ylabel('$T_{\mathrm{2\, m}}$ (K)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/tas.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(clat,2), '-k');
    xlabel('latitude (deg)');
    ylabel('$\cos(\phi)$ (unitless)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/clat.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(dtdlat,2), '-k');
    xlabel('latitude (deg)');
    ylabel('$\cos(\phi)\frac{\mathrm{d}T}{\mathrm{d}\phi}$ (K)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/dtdlat.png', plotdir), '-dpng', '-r300');
    close;
    
    % divFm vs lapt regression
    % X = [ones(length(lapt(:)),1) lapt(:)];
    % beta = X\(-divfm(:));
    beta = lapt(:)\(-divfm(:));
    
    % regression line
    tvec = linspace(-300,300,100);
    % divfmreg = beta(1) + beta(2)*tvec;
    divfmreg = beta*tvec;
    Rsq = 1 - sum((divfm(:) - (beta*-lapt(:))).^2)/sum((divfm(:) - mean(divfm(:))).^2);
    
    % compute diffusivity
    int_divfm = cumtrapz(deg2rad(grid.dim2.lat), clat.*divfm, 1);
    D = int_divfm./dtdlat;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(int_divfm,2), '-k');
    xlabel('latitude (deg)');
    ylabel('$\int\cos(\phi)\nabla\cdot F_m$ (W/m)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/int_divfm.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(D,2), '-k');
    % plot(tvec, olrreg, '-r');
    % text(210,300,sprintf('$\\mathrm{OLR=%.2f+%.2f}T$', beta(1), beta(2)));
    xlabel('latitude (deg)');
    ylabel('$D (W/K)$');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-2 2], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/diffusivity.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(-lapt(:), divfm(:), '.k', 'markersize', 1);
    plot(tvec, divfmreg, '-r');
    text(-350,150,sprintf('$\\nabla\\cdot F_m=%.2f\\frac{1}{\\cos(\\phi)}\\frac{\\partial}{\\partial \\phi}\\left(-\\cos(\\phi) \\frac{\\partial T}{\\partial \\phi} \\right)$', beta));
    text(100,-170,sprintf('$R^2 = %.2f$', Rsq));
    xlabel('$\frac{1}{\cos(\phi)}\frac{\partial}{\partial \phi}\left(-\cos(\phi) \frac{\partial T}{\partial \phi} \right)$ (K)');
    ylabel('$\mathrm{\nabla\cdot F_m}$ (W m$^{-2}$)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-400 400], 'yscale', 'linear', 'ylim', [-200 200], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/divfm_lapt.png', plotdir), '-dpng', '-r300');
    close;
    
    % analyze with extratropics only
    lapt_et = lapt(grid.dim2.lat<=-25 | grid.dim2.lat>=25,:);
    divfm_et = divfm(grid.dim2.lat<=-25 | grid.dim2.lat>=25,:);

    beta_et = lapt_et(:)\(-divfm_et(:));
    % regression line
    divfmreg_et = beta_et*tvec;
    Rsq_et = 1 - sum((divfm_et(:) - (beta_et*-lapt_et(:))).^2)/sum((divfm_et(:) - mean(divfm_et(:))).^2);
    
    figure(); clf; hold all; box on;
    plot(-lapt_et(:), divfm_et(:), '.k', 'markersize', 1);
    plot(tvec, divfmreg_et, '-r');
    text(-350,150,sprintf('$\\nabla\\cdot F_m=%.2f\\frac{1}{\\cos(\\phi)}\\frac{\\partial}{\\partial \\phi}\\left(-\\cos(\\phi) \\frac{\\partial T}{\\partial \\phi} \\right)$', beta_et));
    text(100,-170,sprintf('$R^2 = %.2f$', Rsq_et));
    xlabel('$\frac{1}{\cos(\phi)}\frac{\partial}{\partial \phi}\left(-\cos(\phi) \frac{\partial T}{\partial \phi} \right)$ (K)');
    ylabel('$\mathrm{\nabla\cdot F_m}$ (W m$^{-2}$)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-400 400], 'yscale', 'linear', 'ylim', [-200 200], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/divfm_lapt_et.png', plotdir), '-dpng', '-r300');
    close;
    
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, nanmean(-beta_et*lapt,2), '-k');
    plot(grid.dim2.lat, nanmean(divfm,2), '-', 'color', par.maroon);
    xlabel('latitude (deg)');
    ylabel('Energy flux (W m$^{-2}$)');
    make_title_type(type, par);
    %set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on')
    print(sprintf('%s/divfm_lapt/lapt.png', plotdir), '-dpng', '-r300');
    close;
    
end
