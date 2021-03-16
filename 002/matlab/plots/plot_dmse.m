% function plot_dmse(ra_lat, res_lat, lh_lat, sh_lat, leg, ymin, ymax, type, fw, par, tend_lat)
function plot_dmse(ra_lat, res_lat, lh_lat, sh_lat, leg, ymin, ymax, type, fw, par, ra_std, res_std, lh_std, sh_std, ra_std2, res_std2, lh_std2, sh_std2)

    % if nargin == 11
    %     bool_tend = 1;
    % else
    %     bool_tend = 0;
    % end
    bool_tend = 0;

    if nargin > 14
        drawspr = 'range';
        plotname = sprintf('%s/0_mon_mse%s_range', par.folder, leg);
    elseif nargin > 10
        drawspr = 'std';
        plotname = sprintf('%s/0_mon_mse%s_std', par.folder, leg);
    else
        drawspr = 'none';
        plotname = sprintf('%s/0_mon_mse%s', par.folder, leg);
    end

    figure(); clf; hold all; box on;

    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    ra=plot([1:12],  circshift(ra_lat, par.shiftby, 2), 'color', 0.5*[1 1 1]);
    if bool_tend
        tend=plot([1:12],  circshift(tend_lat, par.shiftby, 2), 'color', 'k');
    end
    res=plot([1:12], circshift(res_lat,par.shiftby,2), 'color', par.maroon);
    if strcmp(fw, 'ceresrad')
        stf=plot([1:12], circshift(lh_lat,par.shiftby,2), '--', 'color', par.blue);
        stf=plot([1:12], circshift(lh_lat,par.shiftby,2), ':', 'color', par.orange);
    else
        lhf=plot([1:12], circshift(lh_lat,par.shiftby,2), 'color', par.blue);
        shf=plot([1:12], circshift(sh_lat,par.shiftby,2), 'color', par.orange);
    end

    if ~strcmp(drawspr, 'none')
        mon2 = [1:12, fliplr(1:12)];

        if strcmp(drawspr, 'std')
            ra_std2 = [circshift(ra_lat+ra_std, par.shiftby, 2), fliplr(circshift(ra_lat-ra_std, par.shiftby, 2))];
            res_std2 = [circshift(res_lat+res_std, par.shiftby, 2), fliplr(circshift(res_lat-res_std, par.shiftby, 2))];
            lh_std2 = [circshift(lh_lat+lh_std, par.shiftby, 2), fliplr(circshift(lh_lat-lh_std, par.shiftby, 2))];
            sh_std2 = [circshift(sh_lat+sh_std, par.shiftby, 2), fliplr(circshift(sh_lat-sh_std, par.shiftby, 2))];
        elseif strcmp(drawspr, 'range')
            ra_std2 = [circshift(ra_std2, par.shiftby, 2), fliplr(circshift(ra_std, par.shiftby, 2))];
            res_std2 = [circshift(res_std2, par.shiftby, 2), fliplr(circshift(res_std, par.shiftby, 2))];
            lh_std2 = [circshift(lh_std2, par.shiftby, 2), fliplr(circshift(lh_std, par.shiftby, 2))];
            sh_std2 = [circshift(sh_std2, par.shiftby, 2), fliplr(circshift(sh_std, par.shiftby, 2))];
        end

        fill(mon2, ra_std2, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
        fill(mon2, res_std2, par.maroon, 'facealpha', 0.2, 'edgealpha', 0.2);
        fill(mon2, lh_std2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
        fill(mon2, sh_std2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
    end
    
    if ~strcmp(leg, "_legonly")
        if isfield(par, 'lat_center')
            make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
        else
            make_title_type_lat(type, par.lat_bound, par.lat_pole, par);
        end
    end
    
    % xlabel('Month');
    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
    
    if leg == ""
        if strcmp(fw, 'ceresrad')
            legend([ra res stf], '$R_a$', '$\partial_t h + \nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside', 'orientation', 'vertical');
        else
            if bool_tend
                legend([tend, ra res lhf, shf], '$\partial_t h$', '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'orientation', 'vertical');
            else
                legend([ra res lhf, shf], '$R_a$', '$\partial_t h + \nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'orientation', 'vertical');
            end
        end
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
    elseif leg == "_noleg"
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    end
    set(gca, 'ylim', [ymin ymax], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'yminortick', 'on', 'tickdir', 'out');

    if strcmp(leg, "_legonly")
        axis off;
        axis([10,11,10,11])
        legend([ra res lhf, shf], '$R_a$', '$\partial_t h + \nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'northwest', 'numcolumns', 2, 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 2.9 0.65])
    end

    print(plotname, '-dpng', '-r300');
    close;

end
