clc; close all; clear variables;

addpath(genpath('/home/miyawaki/doc/matlab'));

% global settings
ppos = [0 0 10/3 7/3];
set(0, 'DefaultLineLineWidth', 1.1);
set(0, 'DefaultFigureUnits', 'inches', 'DefaultFigurePosition', ppos);
set(0, 'DefaultTextInterpreter', 'none');
set(0, 'defaultAxesTickLabelInterpreter','none'); 
set(0, 'defaultLegendInterpreter','none');
ppos_wide = [0 0 13/3 7/3];
fs = 9;

% load data
temp = ncread('/project2/tas1/abacus/data1/tas/archive/Reanalysis/ERAinterim_alllevs/T/T1979_2012_jan_avg.nc', 'var130');
lat = ncread('/project2/tas1/abacus/data1/tas/archive/Reanalysis/ERAinterim_alllevs/T/T1979_2012_jan_avg.nc', 'lat');
lev = ncread('/project2/tas1/abacus/data1/tas/archive/Reanalysis/ERAinterim_alllevs/T/T1979_2012_jan_avg.nc', 'lev');
temp_mean = squeeze(nanmean(nanmean(temp,1),4));

figure(); clf; hold all;
cmax = 350; cintu = 350/100; cmp = colCog(2*cmax/cintu);
colormap(cmp);
cint = 2; cintpos = cint:cint:cint*100; cintneg = -fliplr(cintpos);
[C, h] = contourf(lat, lev/100, temp_mean', min(min(temp_mean)):5:max(max(temp_mean)), 'levellist', [200:5:300], 'showtext', 'on', 'textlist', [200:20:300]);
xlabel('latitude (deg N)'); ylabel('pressure (hPa)');
caxis([200 340]);
cb = colorbar('limits', [200 300], 'location', 'eastoutside');
ylabel(cb, 'temperature (K)');
% clabel(C, h, [220:10:280])
set(gcf, 'paperunits', 'inches', 'paperposition', ppos_wide);
set(gca, 'fontsize', fs, 'xtick', [-90:30:90], 'ylim', [200 1000], 'ydir', 'reverse', 'yscale', 'log');
print('jan_temp', '-dpng', '-r300');
close;
