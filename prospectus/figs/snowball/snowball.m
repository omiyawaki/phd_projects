clc; close all; clear variables;

addpath(genpath('/home/miyawaki/doc/matlab'));

% global settings
ppos = [0 0 10/3 7/3];
set(0, 'DefaultLineLineWidth', 1.1);
set(0, 'DefaultFigureUnits', 'inches', 'DefaultFigurePosition', ppos);
set(0, 'DefaultTextInterpreter', 'tex');
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');
ppos_wide = [0 0 13/3 7/3];
fs = 9;
lv = 2.2e6;

navy = 0.2*[0, 0.447, 0.741];
darkblue = 0.5*[0, 0.447, 0.741];
blue = [0, 0.447, 0.741];
orange = [0.85, 0.325, 0.098];
yellow = [0.929, 0.694, 0.125];
purple = [0.494, 0.184, 0.556];
green = [0.466, 0.674, 0.188];
cyan = [0.301, 0.745, 0.933];
maroon = [0.635, 0.078, 0.184];
brown = 0.5*[0.635, 0.078, 0.184];
darkbrown = 0.2*[0.635, 0.078, 0.184];
purple = 0.5*[0.4940, 0.1840, 0.5560];

% load data
file = '/project2/tas1/arejaygraham/6HourSims/rjg_20170915_2/BOT_rjg_20170915_2_0128_0134_djf.nc';
file_ice='/project2/tas1/arejaygraham/6HourSims/rjg_20170908/BOT_rjg_20170908_0078_0083_djf.nc';

lat=ncread(file, 'lat');

swabs=squeeze(nanmean(nanmean(ncread(file,'srad0')-(ncread(file,'srads')),3),1));
swabs_ice=squeeze(nanmean(nanmean(ncread(file_ice,'srad0')-(ncread(file_ice,'srads')),3),1));
ghe=squeeze(nanmean(nanmean(ncread(file,'trad0')-(ncread(file,'trads')),3),1));
ghe_ice=squeeze(nanmean(nanmean(ncread(file_ice,'trad0')-(ncread(file_ice,'trads')),3),1));

ra = swabs + ghe;
ra_ice = swabs_ice + ghe_ice;

lhf=-squeeze(nanmean(nanmean(ncread(file,'ahfl'),3),1));
lhf_ice=-squeeze(nanmean(nanmean(ncread(file_ice,'ahfl'),3),1));

shf=-squeeze(nanmean(nanmean(ncread(file,'ahfs'),3),1));
shf_ice=-squeeze(nanmean(nanmean(ncread(file_ice,'ahfs'),3),1));

ei=squeeze(nanmean(nanmean(ncread(file,'trad0')+ncread(file,'srad0')-(ncread(file,'srads')+ncread(file,'trads')+ncread(file,'ahfs')+ncread(file,'ahfl')),3),1))';
ei_ice=squeeze(nanmean(nanmean(ncread(file_ice,'trad0')+ncread(file_ice,'srad0')-(ncread(file_ice,'srads')+ncread(file_ice,'trads')+ncread(file_ice,'ahfs')+ncread(file_ice,'ahfl')),3),1))';

% moist rce criteria
moist_conv = lhf+shf;
moist_rce_lat_s = interp1(( moist_conv(lat<0) + ra(lat<0) ), lat(lat<0), 0)
moist_rce_lat_n = interp1(( moist_conv(lat>0) + ra(lat>0) ), lat(lat>0), 0)

moist_conv_ice = lhf_ice+shf_ice;
moist_rce_ice_lat_n = interp1(( moist_conv_ice(lat>0 & lat<18) + ra_ice(lat>0 & lat<18) ), lat(lat>0 & lat<18), 0)
moist_rce_ice_lat_s = interp1(( moist_conv_ice(lat<-18 & lat>-60) + ra_ice(lat<-18 & lat>-60) ), lat(lat<-18 & lat>-60), 0)

dry_conv_ice = shf;
dry_rce_ice_lat_s = interp1(( dry_conv_ice(lat<0) + ra(lat<0) ), lat(lat<0), 0)
dry_rce_ice_lat_n = interp1(( dry_conv_ice(lat>0) + ra(lat>0) ), lat(lat>0), 0)

% rae criteria
rae_lat = interp1(( ei(lat>40) - ra(lat>40)' ), lat(lat>40), 0)
rae_ice_lat = interp1(( ei_ice(lat>30 & lat<60) - ra_ice(lat>30 & lat<60)' ), lat(lat>30 & lat<60), 0)

% MSE plots
figure(); clf; hold all;
plot(lat, ra, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lhf, '-', 'color', blue);
plot(lat, shf, '-', 'color', orange);
plot(lat, ei, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
line([moist_rce_lat_s moist_rce_lat_s], [-200 200], 'linestyle', '-', 'color', 'k');
line([moist_rce_lat_n moist_rce_lat_n], [-200 200], 'linestyle', '-', 'color', 'k');
line([rae_lat rae_lat], [-200 200], 'linestyle', '--', 'color', 'k');
hline(0, '-k');
legend('R_a', 'LH', 'SH', '\nabla\cdot(F_a)', 'location', 'eastoutside');
print('modern_energy_budget', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra+lhf+shf, '-k');
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
hline(50, '--k');
hline(-50, '--k');
print('modern_rce_balance', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra'-ei, '-k');
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
print('modern_rae_balance', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lhf, '-', 'color', blue);
plot(lat, shf, '-', 'color', orange);
plot(lat, ei, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
legend('R_a', 'LHF', 'SHF', '\nabla\cdot(F_a)', 'location', 'eastoutside');
print('modern_energy_budget_no_annotate', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra, '--', 'color', 0.3*[ 1 1 1]);
plot(lat, lhf, '--', 'color', blue);
plot(lat, shf, '--', 'color', orange);
plot(lat, ei, '--', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
legend('R_a', 'LHF', 'SHF', '\nabla\cdot(F_a)', 'location', 'eastoutside');
print('modern_energy_budget_dashed', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, lhf, '-', 'color', blue);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('modern_lhf', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, shf, '-', 'color', orange);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('modern_shf', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ei, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('modern_ei', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, abs(ei./(lhf')), '-', 'color', 0.3*[ 1 1 1]);
xlabel('latitude (deg N)'); ylabel('dimensionless');
title('| \nabla\cdot(F_a) / LHF |');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.01 100], 'yscale', 'log');
hline(0.1, '-k');
hline(0.5, '-k');
hline(10, '--k');
print('modern_lhf_fa_ratio', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, abs(ei./(shf')), '-', 'color', 0.3*[ 1 1 1]);
xlabel('latitude (deg N)'); ylabel('dimensionless');
title('| \nabla\cdot(F_a) / SHF |');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.01 100], 'yscale', 'log');
hline(0.1, '-k');
hline(10, '--k');
print('modern_shf_fa_ratio', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, abs(ei./(lhf'+shf')), '-', 'color', 0.3*[ 1 1 1]);
xlabel('latitude (deg N)'); ylabel('dimensionless');
title('| \nabla\cdot(F_a) / (LHF + SHF) |');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.01 100], 'yscale', 'log');
hline(0.1, '-k');
hline(0.5, '-k');
hline(10, '--k');
print('modern_lhf+shf_fa_ratio', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra_ice, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lhf_ice, '-', 'color', blue);
plot(lat, shf_ice, '-', 'color', orange);
plot(lat, ei_ice, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
line([moist_rce_ice_lat_s moist_rce_ice_lat_s], [-200 200], 'linestyle', '-', 'color', 'k');
line([moist_rce_ice_lat_n moist_rce_ice_lat_n], [-200 200], 'linestyle', '-', 'color', 'k');
line([rae_ice_lat rae_ice_lat], [-200 200], 'linestyle', '--', 'color', 'k');
hline(0, '-k');
legend('R_a', 'LH', 'SH', '\nabla\cdot(F_a)', 'location', 'eastoutside');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-60 60]);
print('snowball_energy_budget', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra_ice, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lhf_ice, '-', 'color', blue);
plot(lat, shf_ice, '-', 'color', orange);
plot(lat, ei_ice, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
hline(0, '-k');
legend('R_a', 'LH', 'SH', '\nabla\cdot(F_a)', 'location', 'eastoutside');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-60 60]);
print('snowball_energy_budget_no_annotate', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, lhf, '--', 'color', blue);
plot(lat, lhf_ice, '-', 'color', blue);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('snowball_lhf', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, shf, '--', 'color', orange);
plot(lat, shf_ice, '-', 'color', orange);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('snowball_shf', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ei, '--', 'color', maroon);
plot(lat, ei_ice, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-130 150]);
hline(0, '-k');
print('snowball_ei', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, lhf_ice+ra_ice, '-', 'color', 0.3*[ 1 1 1]);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
title('LH - R_a');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
hline(2, '--k');
print('snowball_rce_criteria', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, abs(ei_ice./(lhf_ice'+shf_ice')), '-', 'color', 0.3*[ 1 1 1]);
xlabel('latitude (deg N)'); ylabel('dimensionless');
title('| \nabla\cdot(F_a) / (LH + SH) |');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.01 100], 'yscale', 'log');
hline(0.1, '-k');
hline(10, '--k');
print('snowball_lhf_fa_ratio', '-dpng', '-r300');
close;

% DSE framework
precip=squeeze(nanmean(nanmean(ncread(file,'precip'),3),1));
precip_ice=squeeze(nanmean(nanmean(ncread(file_ice,'precip'),3),1));

lp = precip*lv;
lp_ice = precip_ice*lv;

dse=squeeze(nanmean(nanmean(ncread(file,'trad0')+ncread(file,'srad0')-(ncread(file,'srads')+ncread(file,'trads')+ncread(file,'ahfs')-lv*ncread(file,'precip')),3),1))';
dse_ice=squeeze(nanmean(nanmean(ncread(file_ice,'trad0')+ncread(file_ice,'srad0')-(ncread(file_ice,'srads')+ncread(file_ice,'trads')+ncread(file_ice,'ahfs')-lv*ncread(file_ice,'precip')),3),1))';

% DSE plots
figure(); clf; hold all;
plot(lat, ra, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lp, '-', 'color', blue);
plot(lat, shf, '-', 'color', orange);
plot(lat, dse, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
legend('R_a', 'LP', 'SH', '\nabla\cdot(vs)', 'location', 'eastoutside');
print('modern_dse', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra+lp+shf, '-k');
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
title('Net balance (DJF)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
hline(50, '--k');
hline(-50, '--k');
print('modern_dse_rce_balance', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra'-dse, '-k');
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
hline(0, '-k');
print('modern_dse_rae_balance', '-dpng', '-r300');
close;

figure(); clf; hold all;
plot(lat, ra_ice, '-', 'color', 0.3*[ 1 1 1]);
plot(lat, lp_ice, '-', 'color', blue);
plot(lat, shf_ice, '-', 'color', orange);
plot(lat, dse_ice, '-', 'color', maroon);
xlabel('latitude (deg N)'); ylabel('energy flux (W/m^2)');
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xlim', [-90 90], 'xtick', [-90:30:90]);
line([moist_rce_lat_s moist_rce_lat_s], [-200 200], 'linestyle', '-', 'color', 'k');
line([moist_rce_lat_n moist_rce_lat_n], [-200 200], 'linestyle', '-', 'color', 'k');
line([rae_lat rae_lat], [-200 200], 'linestyle', '--', 'color', 'k');
hline(0, '-k');
legend('R_a', 'LP', 'SH', '\nabla\cdot(vs)', 'location', 'eastoutside');
print('snowball_dse', '-dpng', '-r300');
close;

% compute tropical average of RCE distance
bal = ra+lp+shf;
bal_wgt = cosd(lat)'.*bal;
bal_glob = nansum(bal_wgt)/nansum(cosd(lat)')
[~, idx_south] = min(abs(lat + 30))
[~, idx_north] = min(abs(lat - 30))
lat_south = lat(idx_south)
lat_north = lat(idx_north)
bal_trop = nansum(bal_wgt(idx_north:idx_south))/nansum(cosd(lat(idx_north:idx_south))')

% compute tropical average of MSE distance
bal = ra+lhf+shf;
bal_wgt = cosd(lat)'.*bal;
bal_glob = nansum(bal_wgt)/nansum(cosd(lat)')
[~, idx_south] = min(abs(lat + 30))
[~, idx_north] = min(abs(lat - 30))
lat_south = lat(idx_south)
lat_north = lat(idx_north)
bal_trop = nansum(bal_wgt(idx_north:idx_south))/nansum(cosd(lat(idx_north:idx_south))')
