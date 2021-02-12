function plot_dmse_comp(dmse, dmse_std, dmse_era5c, dmse_merra, dmse_jra55, ylab, varname, leg, type, fw, par)

    figure(); clf; hold all; box on;
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    cmip5=plot([1:12],  circshift(dmse, par.shiftby, 2), 'color', 'k');
    mon2 = [1:12, fliplr(1:12)];
    dmse2 = [circshift(dmse+dmse_std, par.shiftby, 2), fliplr(circshift(dmse-dmse_std, par.shiftby, 2))];
    fill(mon2, dmse2, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
    era5c=plot([1:12],  circshift(dmse_era5c, par.shiftby, 2), 'color', par.blue);
    merra=plot([1:12],  circshift(dmse_merra, par.shiftby, 2), 'color', par.orange);
    jra55=plot([1:12],  circshift(dmse_jra55, par.shiftby, 2), 'color', par.green);
    
    if isfield(par, 'lat_center')
        make_title_lat(par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
    else
        make_title_lat(par.lat_bound, par.lat_pole, par);
    end
    
    % xlabel('Month');
    ylabel(ylab);
    
    if leg == ""
        legend([cmip5 era5c merra jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
    elseif leg == "_noleg"
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    end
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/0_mon_%s%s', par.folder, varname, leg), '-dpng', '-r300');
    close;

end
