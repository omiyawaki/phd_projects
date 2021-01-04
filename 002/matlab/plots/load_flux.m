function [flux_zt, vh, vh_mon, lat, plotdir] = load_flux(type, par)
    % load processed data/proc
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/flux_zt.mat', prefix_proc));
    load(sprintf('%s/vh.mat', prefix_proc));
    load(sprintf('%s/vh_mon.mat', prefix_proc));
end
