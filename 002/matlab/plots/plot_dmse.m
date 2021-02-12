function plot_dmse(ra_lat, res_lat, lh_lat, sh_lat, leg, ymin, ymax, type, fw, par, tend_lat)

    if nargin == 11
        bool_tend = 1;
    else
        bool_tend = 0;
    end

    figure(); clf; hold all; box on;
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    ra=plot([1:12],  circshift(ra_lat, par.shiftby, 2), 'color', 0.5*[1 1 1]);
    if bool_tend
        tend=plot([1:12],  circshift(tend_lat, par.shiftby, 2), 'color', 'k');
    end
    res=plot([1:12], circshift(res_lat,par.shiftby,2), 'color', par.maroon);
    lhf=plot([1:12], circshift(lh_lat,par.shiftby,2), 'color', par.blue);
    shf=plot([1:12], circshift(sh_lat,par.shiftby,2), 'color', par.orange);
    
    if isfield(par, 'lat_center')
        make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
    else
        make_title_type_lat(type, par.lat_bound, par.lat_pole, par);
    end
    
    % xlabel('Month');
    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
    
    if leg == ""
        if bool_tend
            legend([tend, ra res lhf, shf], '$\partial_t h$', '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'orientation', 'vertical');
        else
            legend([ra res lhf, shf], '$R_a$', '$\partial_t h + \nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'orientation', 'vertical');
        end
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
    elseif leg == "_noleg"
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    end
    set(gca, 'ylim', [ymin ymax], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/0_mon_mse%s', par.folder, leg), '-dpng', '-r300');
    close;

end
