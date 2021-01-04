function plot_mlev(type, par)
    make_dirs(type, par)

    if strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    else
        error('This code only works for data output in the model vertical grid.')
    end
    load(sprintf('%s/grid.mat', prefix));

    % compute sigma from a and b
    ps_vert = repmat(ps_orig, [1 1 1 length(grid.dim3.a)]); % dims (lon x lat x time x plev)
    ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
    b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
    pa = a + b.*ps_vert;

    % time and zonal mean
    pa = squeeze(nanmean(nanmean(pa,1),4));

    % standard pl grid
    pl = 1e2*[1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1 0.5 0.2 0.1];

    % pressure x lat
    figure(); clf; hold all; box on;
    for lev = 1:length(grid.dim3.a)
        h1=plot(grid.dim3.lat, 1e-2*pa(:,lev), '-k', 'linewidth', 0.5);
    end
    for lev = 1:length(pl); plev = pl(lev);
        h2=line([-90 90], 1e-2*plev*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
    end
    legend([h1 h2], 'Model levels', 'Standard pressure levels', 'location', 'southoutside')
    xlabel('latitude (deg)'); ylabel('p (hPa)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on')
    print(sprintf('%s/mlev/mlev', par.plotdir), '-dpng', '-r300');
    close;

end
