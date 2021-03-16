clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

figure_params;
echam_info;

expr = 'rp000133';

if strcmp(expr, 'rp000172')
    prefix = sprintf('/project2/tas1/echam-aiv_rcc_6.1.00p1/experiments/%s', expr);
    yr_vec = 1:14;
else
    prefix = sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s', expr);
    yr_vec = 1:39;
end

foldername = sprintf('/project2/tas1/miyawaki/projects/002/figures/echam/%s/native/eq_check', expr);
if ~exist(foldername, 'dir');
    mkdir(foldername);
end


yr_list = strtrim(cellstr(num2str(yr_vec', '%04.f'))');
toa_net_ts = [];
ts_ts = [];

for y = 1:length(yr_list); yr = yr_list{y};

    lat = ncread(sprintf('%s/BOT_%s_%s.nc', prefix, expr, yr), 'lat');
    olr = ncread(sprintf('%s/BOT_%s_%s.nc', prefix, expr, yr), 'trad0');
    toa_sw = ncread(sprintf('%s/BOT_%s_%s.nc', prefix, expr, yr), 'srad0');
    ts = ncread(sprintf('%s/BOT_%s_%s.nc', prefix, expr, yr), 'tsurf');
    toa_net = olr + toa_sw;
    clear olr toa_sw

    toa_net = zonmean(toa_net, 1); % zonal mean
    ts = zonmean(ts, 1); % zonal mean

    clat = cosd(lat);
    clata = repmat(clat, [1 size(toa_net, 2)]);
    toa_net = mermean(toa_net, 1, clata, clat);
    toa_net_am = nanmean(toa_net);
    ts = mermean(ts, 1, clata, clat);
    ts_am = nanmean(ts);

    toa_net_ts = [toa_net_ts toa_net_am];
    ts_ts = [ts_ts ts_am];
    size(toa_net_ts)

end

figure(); clf; hold on; box on;
line([1 length(toa_net_ts)], [0 0], 'color', 'k', 'linewidth', 0.5);
plot([1:length(toa_net_ts)], toa_net_ts, 'k');
xlabel('yr'); ylabel('$R_T$ (W m$^{-2}$)');
title(par.echam.(expr));
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
set(gca, 'xlim', [1 length(toa_net_ts)], 'ylim', [-2 6], 'xminortick', 'on', 'yminortick', 'on')
print(sprintf('%s/net_toa.png', foldername), '-r300', '-dpng');
close;

figure(); clf; hold on; box on;
% line([1 length(toa_net_ts)], [0 0], 'color', 'k', 'linewidth', 0.5);
plot([1:length(ts_ts)], ts_ts, 'k');
xlabel('yr'); ylabel('$T_s$ (K)');
title(par.echam.(expr));
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
set(gca, 'xlim', [1 length(toa_net_ts)], 'xminortick', 'on', 'yminortick', 'on')
print(sprintf('%s/ts.png', foldername), '-r300', '-dpng');
close;

function varzm = zonmean(var0, idim)
    varzm = squeeze(nanmean(var0, idim));
end

function varmm = mermean(var0, idim, clata, clatv)
    varmm = nansum(clata.*var0, idim)/nansum(clatv);
end