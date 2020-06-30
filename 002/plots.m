clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.ep_swp = 0.3; %[0.25 0.3 0.35]; % threshold value for determining RCE and RAE
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget (option for ERA-Interim)
% if == stf, use SH and LH data from ERA-Interim to close energy budget (option for ERA-Interim)
% if == era5, use ERA5 radiative cooling and surface turbulent fluxes to close energy budget (only option for ERA5)
par.closure = 'era5';
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
% set default figure parameters
if 1
    par.ppos = [0 0 10/3 7/3];
    par.ppos_sq = [0 0 10/3 10/3];
    par.ppos_wide = [0 0 13/3 7/3];
    par.ppos_verywide = [0 0 16/3 7/3];
    par.fs = 10;
    set(0, 'DefaultLineLineWidth', 1.1);
    set(0, 'DefaultFigureUnits', 'inches', 'DefaultFigurePosition', par.ppos_wide);
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    par.navy = 0.2*[0, 0.447, 0.741];
    par.darkblue = 0.5*[0, 0.447, 0.741];
    par.blue = [0, 0.447, 0.741];
    par.orange = [0.85, 0.325, 0.098];
    par.yellow = [0.929, 0.694, 0.125];
    par.purple = [0.494, 0.184, 0.556];
    par.green = [0.466, 0.674, 0.188];
    par.cyan = [0.301, 0.745, 0.933];
    par.maroon = [0.635, 0.078, 0.184];
    par.brown = 0.5*[0.635, 0.078, 0.184];
    par.darkbrown = 0.2*[0.635, 0.078, 0.184];
    par.purple = 0.5*[0.4940, 0.1840, 0.5560];
    par.gray = 0.5*[1 1 1];
end
gcm_info

%% call functions
% plot_rad_lat(par)
% plot_rad_lon_lat(par)
plot_tediv_lat(par)

type = 'era5';
% choose_plots(type, par);
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % choose_plots(type, par);
end

% sweep through various threshold values
for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
    type = 'era5';
    % choose_plots_ep(type, par)
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        % choose_plots_ep(type, par)
    end
end

%% define functions
function choose_plots(type, par)
    plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
end
function plot_energy_lat(type, par) % latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, par] = load_flux(type, par); % load data

    for f = {'mse', 'dse'}; fw = f{1};
        for l = {'lo', 'l', 'o'}; land = l{1};
            % northward M/DSE transport, mon x lat
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
            figure(); clf; hold all;
            cmp = colCog(30);
            colormap(cmp);
            [C, h] = contour(mesh_lat, mesh_mon, vh_mon.(land).(fw)*10^-15, -5:1:5);
            clabel(C, h, [-4:2:4], 'fontsize', 6, 'interpreter', 'latex');
            caxis([-5 5]);
            xlabel('Month'); ylabel('latitude (deg)');
            title(sprintf('Northward %s Transport (PW)', upper(fw)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'fontsize', par.fs, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on')
            print(sprintf('%s/transport/%s/all/%s', par.plotdir, land, fw), '-dpng', '-r300');
            close;
        end
    end

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        for f = {'mse', 'dse'}; fw = f{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                plot(lat, flux.(land).(time).ra, 'color', par.gray); text(0, 0.75*interp1(lat,flux.(land).(time).ra,0), '\boldmath{$R_a$}', 'color', par.gray);
                if strcmp(fw, 'mse'); plot(lat,flux.(land).(time).res.mse, 'color', par.maroon); text(-42, 2*interp1(lat,flux.(land).(time).res.mse,-42), '\boldmath{$\nabla\cdot F_m$}', 'color', par.maroon);
                elseif strcmp(fw, 'dse'); plot(lat,flux.(land).(time).res.dse, '--', 'color', par.maroon); text(-30, 2*interp1(lat,flux.(land).(time).res.dse,-30), '\boldmath{$\nabla\cdot F_s$}', 'color', par.maroon);
                end
                if strcmp(type, 'era5') | strcmp(type, 'erai')
                    if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                    elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue);
                    end
                    plot(lat, -flux.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                elseif strcmp(type, 'gcm')
                    if strcmp(fw, 'mse'); plot(lat, flux.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                    elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue);
                    end
                    plot(lat, flux.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                end
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                end
                xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf 150]); end
                print(sprintf('%s/energy-flux/%s/%s/%s-all', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;

                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                line([-90 90], [1 1]*(1-par.ep), 'linewidth', 0.5, 'color', par.blue);
                if strcmp(fw, 'mse'); plot(lat,flux.(land).(time).r1.mse, '-k');
                elseif strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                end
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                end
                xlabel('latitude (deg)'); ylabel('$R_1$ (unitless)');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                print(sprintf('%s/energy-flux/%s/%s/%s-r1', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;
            % northward M/DSE transport
                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([0 0], [min(vh.(land).(time).(fw)) max(vh.(land).(time).(fw))]*10^-15, 'linewidth', 0.5, 'color', 'k');
                if strcmp(fw, 'mse'); plot(lat, vh.(land).(time).(fw)*10^-15, 'color', par.maroon);
                elseif strcmp(fw, 'dse'); plot(lat, vh.(land).(time).(fw)*10^-15, '--', 'color', par.maroon);
                end
                xlabel('latitude (deg)'); ylabel('PW')
                title(sprintf('Northward %s Transport, %s', upper(fw), upper(time)));
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                if strcmp(fw, 'dse'); set(gca, 'ytick', [-5:5]); end;
                print(sprintf('%s/transport/%s/%s/%s', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;
            % MSE/DSE transport plotted together
                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([0 0], [min([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse]) max([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse])]*10^-15, 'linewidth', 0.5, 'color', 'k');
                h_mse = plot(lat, vh.(land).(time).mse*10^-15, 'color', par.maroon);
                h_dse = plot(lat, vh.(land).(time).dse*10^-15, '--', 'color', par.maroon);
                h_lh = plot(lat, (vh.(land).(time).mse-vh.(land).(time).dse)*10^-15, ':', 'color', par.maroon);
                legend([h_mse h_dse h_lh], '$F_m$', '$F_s$', '$F_{m}-F_{s}$', 'location', 'eastoutside');
                xlabel('latitude (deg)'); ylabel('PW')
                title(sprintf('Northward Energy Transport, %s', upper(time)));
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                print(sprintf('%s/transport/%s/%s/all', par.plotdir, land, time), '-dpng', '-r300');
                close;
            end % land
        end % end mse/dse loop

    end % time
end

function choose_plots_ep(type, par)
    plot_rcae(type, par) % plot RCE/RAE regimes
    % plot_temp(type, par) % plot temperature profiles
    % plot_ma_diff(type, par) % plot difference of temperature profile from moist adiabat
end
function plot_rcae(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g/rcae_z.mat', prefix_proc, par.lat_interp, par.ep)); % load lat x mon RCAE data
    load(sprintf('%s/%s/eps_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for f = {'mse', 'dse'}; fw = f{1};
            for fn = fieldnames(rcae_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'mse'); crit_text = '$F_m < \max(|F_m|)/2$ flag';
                    elseif strcmp(fw, 'mse'); crit_text = '$F_s < \max(|F_s|)/2$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'mse'); crit_text = '$F_m < \max(|F_m|)/4$ flag';
                    elseif strcmp(fw, 'mse'); crit_text = '$F_s < \max(|F_s|)/4$ flag'; end;
                end
                % lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_z.(land).(fw).(crit));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s'), upper(type), crit_text, land_text);
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, crit_text, land_text)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g/%s/%s/%s/0_rcae_mon_lat', par.plotdir, par.ep, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g/%s/%s/%s/%s/rcae_lat_lon', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_temp(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/ta.mat', type, par.lat_interp, par.ep));
        plev = grid.dim3.plev;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ta.mat', type, par.model, par.lat_interp, par.ep));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ma.mat', type, par.model, par.lat_interp, par.ep));
        plev = grid.dim3.plev/100;
    end
    for f = {'mse', 'dse'}; fw = f{1};
        for c = fieldnames(ta.rce.tp.(fw))'; crit = c{1};
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for l = {'lo', 'l', 'o'}; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

                % RCE and RAE separated into NH and SH
                    figure(); clf; hold all;
                    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
                    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
                    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
                    h_rae_nh = plot(ta.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
                    h_rae_sh = plot(ta.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('%s, %s', upper(time), land_text));
                    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE', 'NH ML RCE', 'SH ML RCE', 'NH RAE', 'SH RAE', 'location', 'northeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rcae_all', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time)-ma.rce.tp.(fw).(crit).(land).(time).ta, plev, 'color', par.maroon);
                    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time)-ma.rce.nh.(fw).(crit).(land).(time).ta, plev, 'color', par.orange);
                    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time)-ma.rce.sh.(fw).(crit).(land).(time).ta, plev, '--', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [100 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rce_diff', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % Tropical RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
                    h_rce_ma = plot(ma.rce.tp.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.maroon);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rce_tp', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % NH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, 'color', par.orange);
                    h_rce_ma = plot(ma.rce.nh.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('NH ML RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rce_nh', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % SH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, 'color', par.orange);
                    h_rce_ma = plot(ma.rce.sh.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('SH ML RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rce_sh', par.plotdir, par.ep, fw, crit, land, time), '-dpng', '-r300');
                    close;
                end % land
                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    h_rce_l = plot(ta.rce.nh.(fw).(crit).l.(time)-ma.rce.nh.(fw).(crit).l.(time).ta, plev, 'color', par.orange);
                    h_rce_o = plot(ta.rce.nh.(fw).(crit).o.(time)-ma.rce.nh.(fw).(crit).o.(time).ta, plev, ':', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    legend([h_rce_l, h_rce_o], 'Land', 'Ocean', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [100 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g/%s/%s/%s/%s/temp/rce_lo_diff', par.plotdir, par.ep, fw, crit, 'lo', time), '-dpng', '-r300');
                    close;
            end % time avg
        end % RCE/RAE definition
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
    figure(); clf; hold all;
    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
    h_rae_nh = plot(ta.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
    h_rae_sh = plot(ta.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
    xlabel('T (K)'); ylabel('p (hPa)');
    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE ($\pm 30^\circ$)', 'NH ML RCE ($>+30^\circ$)', 'SH ML RCE ($<-30^\circ$)', 'NH RAE', 'SH RAE', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp', par.plotdir), '-dpng', '-r300');
    close;
    % Legend for moist adiabat comparisons
    figure(); clf; hold all;
    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, '-k');
    h_rce_tp_ma = plot(ma.rce.tp.(fw).(crit).(land).(time).ta, plev, ':k');
    xlabel('T (K)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        legend([h_rce_tp h_rce_tp_ma], upper(type), 'Moist adiabat', 'location', 'eastoutside');
    elseif strcmp(type, 'gcm')
        legend([h_rce_tp h_rce_tp_ma], par.model, 'Moist adiabat', 'location', 'eastoutside');
    end
    title(upper(sprintf('%s', time)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp_ma', par.plotdir), '-dpng', '-r300');
    close;
end
function plot_ma_diff(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/ta_mon_lat.mat', type, par.lat_interp, par.ep));
        plev = grid.dim3.plev;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ta_mon_lat.mat', type, par.model, par.lat_interp, par.ep));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ma_mon_lat.mat', type, par.model, par.lat_interp, par.ep));
        plev = grid.dim3.plev/100;
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for plev_eval = [300:100:500] % evaluation plev in hPa
            diff = permute(ta.(land) - ma.(land).ta, [3 1 2]); % bring plev to front
            diff = squeeze(interp1(plev, diff, plev_eval)); % evaluate difference at plev_eval
            % lat x lon of RCE and RAE
            figure(); clf; hold all;
            cmp = colCog(40);
            colormap(cmp);
            imagesc([1:12], [lat(1) lat(end)], diff);
            caxis([-20 20]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-20 20], 'ytick', [-20:5:20], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$(T - T_m)_{%g \\,\\mathrm{hPa}}$ (K)', plev_eval));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ma_diff/plev_%g/%s/ma_diff_lat_lon', par.plotdir, plev_eval, land), '-dpng', '-r300');
            close;
        end
    end
end

function plot_rad_lat(par)
    load('/project2/tas1/miyawaki/projects/002/data/proc/comp/comp_zt.mat');

    for fn = {'tsr', 'ttr', 'net', 'ssr', 'swabs', 'str', 'ra', 'surface_turbulent_plus_LW'}; fname = fn{1};
        if strcmp(fname, 'tsr'); ytext = 'Net SW TOA';
        elseif strcmp(fname, 'ttr'); ytext = 'OLR';
        elseif strcmp(fname, 'net'); ytext = 'Net TOA';
        elseif strcmp(fname, 'ssr'); ytext = 'Net SW SFC';
        elseif strcmp(fname, 'str'); ytext = 'Net LW SFC';
        elseif strcmp(fname, 'swabs'); ytext = 'SWABS';
        elseif strcmp(fname, 'ra'); ytext = '$R_a$';
        elseif strcmp(fname, 'surface_turbulent_plus_LW'); ytext = 'LH + SH + Net LW SFC';
        end

        figure(); clf; hold all;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        h.erai = plot(lat, erai_zt.rad.(fname), 'color', par.yellow);
        h.era5 = plot(lat, era5_zt.rad.(fname), 'color', par.orange);
        if ~strcmp(fname, 'surface_turbulent_plus_LW'); h.ceres = plot(lat, ceres_zt.rad.(fname), 'color', par.green); end;
        if ~any(strcmp(fname, {'str', 'ra'})); h.don = plot(lat, don_zt.(fname), 'color', par.blue); end;
        xlabel('latitude (deg)'); ylabel(sprintf('%s (Wm$^{-2}$)', ytext));
        if any(strcmp(fname, {'str', 'ra'})); legend([h.erai h.era5 h.ceres], 'ERA-I', 'ERA5', 'CERES4.1', 'location', 'eastoutside');
        elseif strcmp(fname, 'surface_turbulent_plus_LW'); legend([h.erai h.era5 h.don], 'ERA-I', 'ERA5', 'DB13', 'location', 'eastoutside');
        else legend([h.erai h.era5 h.ceres h.don], 'ERA-I', 'ERA5', 'CERES4.1', 'DB13', 'location', 'eastoutside'); end;
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s', fname), '-dpng', '-r300');
        close;

        if ~any(strcmp(fname, {'str', 'ra'}));
            figure(); clf; hold all;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            h.erai = plot(lat, erai_zt.rad.(fname) - don_zt.(fname), 'color', par.yellow);
            h.era5 = plot(lat, era5_zt.rad.(fname) - don_zt.(fname), 'color', par.orange);
            if ~strcmp(fname, 'surface_turbulent_plus_LW'); h.ceres = plot(lat, ceres_zt.rad.(fname) - don_zt.(fname), 'color', par.green); end;
            xlabel('latitude (deg)'); ylabel(sprintf('$\\Delta$ (%s) (Wm$^{-2}$)', ytext));
            if ~strcmp(fname, 'surface_turbulent_plus_LW'); legend([h.erai h.era5 h.ceres], 'ERA-I $-$ DB13', 'ERA5 $-$ DB13', 'CERES4.1 $-$ DB13', 'location', 'eastoutside');
            else; legend([h.erai h.era5], 'ERA-I $-$ DB13', 'ERA5 $-$ DB13', 'location', 'eastoutside'); end;
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s-diff-db13', fname), '-dpng', '-r300');
            close;
        else
            figure(); clf; hold all;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            h.erai = plot(lat, erai_zt.rad.(fname) - ceres_zt.rad.(fname), 'color', par.yellow);
            h.era5 = plot(lat, era5_zt.rad.(fname) - ceres_zt.rad.(fname), 'color', par.orange);
            xlabel('latitude (deg)'); ylabel(sprintf('$\\Delta$ (%s) (Wm$^{-2}$)', ytext));
            legend([h.erai h.era5], 'ERA-I $-$ CERES4.1', 'ERA5 $-$ CERES4.1', 'location', 'eastoutside')
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s-diff-ceres', fname), '-dpng', '-r300');
            close;
        end

    end

end
function plot_rad_lon_lat(par)
    load('/project2/tas1/miyawaki/projects/002/data/proc/comp/ceres_2001_2009.mat');
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    [mesh_lat, mesh_lon] = meshgrid(grid.ceres.dim2.lon, lat);

    figure(); clf; hold all;
    cmp = colCog(48);
    colormap(cmp);
    contourf(mesh_lat, mesh_lon, ceres_t.rad.ra', [-360 -240:60:-120 -90:30:-30 -15], 'linecolor', 'none');
    contour(par.landlon+180, par.landlat, par.land, [1 1], 'linecolor', 0.3*[1 1 1], 'linewidth', 0.5);
    caxis([-360 360]);
    xlabel('longitude (deg)'); ylabel('latitude (deg)');
    title('$R_a$, CERES4.1 2001--2009');
    cb = colorbar('limits', [-360 -15], 'ytick', [-360 -240:60:-120 -90:30:-30 -15], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$R_a$ (Wm$^{-2}$)'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [0:60:360], 'xminortick', 'on', 'ytick', [-90:30:90], 'yminortick', 'on')
    print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lon_lat/ra_ceres'), '-dpng', '-r300');
    close;

end
function plot_tediv_lat(par)
    par.model = 'MPI-ESM-LR';
    gcm=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/gcm/%s/std/flux_zt.mat', par.model));
    don=load('/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology.mat');
    don_zt.TEDIV = squeeze(nanmean(nanmean(don.TEDIV, 1), 3));
    don_zt.TEDIV = interp1(don.lat, don_zt.TEDIV, gcm.lat);

    figure(); clf; hold all;
    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
    h_gcm = plot(gcm.lat, gcm.flux_zt.lo.ann.res.mse, 'color', par.maroon);
    h_don = plot(gcm.lat, don_zt.TEDIV, 'color', par.blue);
    xlabel('latitude (deg)'); ylabel(sprintf('$\\nabla \\cdot F_m$ (Wm$^{-2}$)'));
    legend([h_gcm h_don], par.model, 'DB13', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/tediv'), '-dpng', '-r300');
    close;

end

function [flux_zt, vh, vh_mon, lat, par] = load_flux(type, par)
    % load processed data/proc
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh_mon.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/flux_zt.mat', type, par.model, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh.mat', type, par.model, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh_mon.mat', type, par.model, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
    end

    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        for f = {'mse', 'dse'}; fw = f{1};
            for j = {'def', 'pe', 'cp', 'w500', 'vh2', 'vh4'}; crit = j{1};
                for l = {'lo', 'l', 'o'}; land = l{1};
                    for k = {'ann', 'djf', 'jja', 'mam', 'son'}; time = k{1};
                        if ~exist(sprintf('%s/eps_%g/%s/%s/%s/%s/temp', par.plotdir, par.ep, fw, crit, land, time), 'dir')
                            mkdir(sprintf('%s/eps_%g/%s/%s/%s/%s/temp', par.plotdir, par.ep, fw, crit, land, time));
                        end
                    end
                end
            end
        end
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        for plev_eval = [300:100:500]
            if ~exist(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land), 'dir')
                mkdir(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land));
            end
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time));
            end
            if ~exist(sprintf('%s/transport/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/transport/%s/%s', par.plotdir, land, time));
            end
        end
    end

    if ~exist(sprintf('%s/legends', par.plotdir), 'dir')
        mkdir(sprintf('%s/legends', par.plotdir));
    end
end
