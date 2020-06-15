clc; close all; clear variables;

% set default figure parameters
if 1
    par.ppos = [0 0 10/3 7/3];
    par.ppos_sq = [0 0 10/3 10/3];
    par.ppos_wide = [0 0 13/3 7/3];
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

% read CERES data
ceres.swdt = ncread('/project2/tas1/miyawaki/projects/002/data/raw/ceres/ceres_swdt_2001.nc', 'solar_mon');
ceres.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/ceres/ceres_swdt_2001.nc', 'lat');

% read ERA data
% era.swdt_raw = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/era_swdt_2001.nc', 'tisr');
% era.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/era_swdt_2001.nc', 'latitude');
era.swdt_raw = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/era_swdt_2001_method_2.nc', 'tisr');
era.lat = ncread('/project2/tas1/miyawaki/projects/002/data/raw/era-interim/rad/era_swdt_2001_method_2.nc', 'latitude');

% account for steps in ERA
for month=1:12
    era.swdt(:,:,month) = ( era.swdt_raw(:,:,month) )/86400;
end

% take zonal and annual averages
ceres.swdt_z = squeeze(nanmean(ceres.swdt, 1));
ceres.swdt_zt = squeeze(nanmean(ceres.swdt_z, 2));
era.swdt_z = squeeze(nanmean(era.swdt, 1));
era.swdt_zt = squeeze(nanmean(era.swdt_z, 2));
% interpolate era profile to ceres grid to compute difference
era.swdt_zti = interp1(era.lat, era.swdt_zt, ceres.lat);

% plot
figure(); clf; hold all;
h_ceres = plot(ceres.lat, ceres.swdt_zt, '-');
h_era = plot(era.lat, era.swdt_zt, '--');
legend([h_ceres h_era], 'CERES', 'ERA', 'location', 'south')
xlabel('latitude (deg)'); ylabel('Energy flux (Wm$^{-2}$)');
title('TOA Incoming SW Flux');
axis('tight');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
print('/project2/tas1/miyawaki/projects/002/figures/ceres/ceres_era_comp', '-dpng', '-r300');
close;

figure(); clf; hold all;
line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
plot(ceres.lat, era.swdt_zti - ceres.swdt_zt, '-k');
legend('ERA $-$ CERES', 'location', 'south')
xlabel('latitude (deg)'); ylabel('Energy flux (Wm$^{-2}$)');
title('TOA Incoming SW Flux');
axis('tight');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
print('/project2/tas1/miyawaki/projects/002/figures/ceres/ceres_era_diff', '-dpng', '-r300');
close;
