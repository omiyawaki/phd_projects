% plot dr1 and its decomposition
function plot_dr1(r1_var, r1_ann_var, dr1_var, comp1, comp2, leg, ymin, ymax, type, fw, par, dr1_std_var, comp1_std, comp2_std, dr1_std_var2, comp1_std2, comp2_std2)

    if nargin > 14
        drawspr = 'range';
        plotname = sprintf('%s/0_mon_dr1z_decomp%s_range', par.folder, leg);
    elseif nargin > 11
        drawspr = 'std';
        plotname = sprintf('%s/0_mon_dr1z_decomp%s_std', par.folder, leg);
    else
        drawspr = 'none';
        plotname = sprintf('%s/0_mon_dr1z_decomp%s', par.folder, leg);
    end

    var_text = '$\Delta R_1$';
    figure(); clf; hold all; box on;
    
    colororder({'k', 'k'});
    
    yyaxis left
    ylim_lo = min([r1_var,ymin]);
    ylim_up = max([r1_var,ymax]);
    tot=plot([1:12], circshift(r1_var,par.shiftby,2), 'k');
    disp(sprintf('Annual mean R_1 is %g', nanmean(r1_var)))
    ylabel(sprintf('$R_1$ (unitless)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    
    yyaxis right
    ylim_lo = -r1_ann_var(1)+ylim_lo;
    ylim_up = -r1_ann_var(1)+ylim_up;
    rcemax = par.ep-r1_ann_var(1);
    if rcemax > ylim_lo
        vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    end
    raemin = par.ga-r1_ann_var(1);
    if raemin < ylim_up
        vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    end
    
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    tot=plot([1:12], circshift(dr1_var,par.shiftby,2), 'k');
    if ~strcmp(drawspr, 'none')
        mon2 = [1:12, fliplr(1:12)];

        if strcmp(drawspr, 'std')
            dr1_std2 = [circshift(dr1_var+dr1_std_var, par.shiftby, 2), fliplr(circshift(dr1_var-dr1_std_var, par.shiftby, 2))];
            comp1_std2 = [circshift(comp1+comp1_std, par.shiftby, 2), fliplr(circshift(comp1-comp1_std, par.shiftby, 2))];
            comp2_std2 = [circshift(comp2+comp2_std, par.shiftby, 2), fliplr(circshift(comp2-comp2_std, par.shiftby, 2))];
        elseif strcmp(drawspr, 'range')
            dr1_std2 = [circshift(dr1_std_var, par.shiftby, 2), fliplr(circshift(dr1_std_var2, par.shiftby, 2))];
            comp1_std2 = [circshift(comp1_std, par.shiftby, 2), fliplr(circshift(comp1_std2, par.shiftby, 2))];
            comp2_std2 = [circshift(comp2_std, par.shiftby, 2), fliplr(circshift(comp2_std2, par.shiftby, 2))];
        end

        fill(mon2, dr1_std2, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
        fill(mon2, comp1_std2, par.maroon, 'facealpha', 0.2, 'edgealpha', 0.2);
        fill(mon2, comp2_std2, 0.5*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
    end
    res=plot([1:12], circshift(dr1_var,par.shiftby,2) - circshift(comp1+comp2,par.shiftby,2), '-.k');
    c1=plot([1:12],  circshift(comp1,par.shiftby,2), '-', 'color', par.maroon);
    c2=plot([1:12],  circshift(comp2,par.shiftby,2), '-', 'color', 0.5*[1 1 1]);
    
    if ~strcmp(leg, "_legonly")
        if isfield(par,'lat_center')
            make_title_type_lat(type, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound, par);
        else
            make_title_type_lat(type, par.lat_bound, par.lat_pole, par);
        end
    end

    ylabel(sprintf('$\\Delta R_1$ (unitless)'));
    
    if leg == ""
        if ~strcmp(fw, 'mse_old')
            ylabel(sprintf('$\\Delta R_1$ (unitless)'));
            legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\partial_t h + \nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\partial_t h + \nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside', 'orientation', 'vertical');
        else
            ylabel(sprintf('$\\Delta R_1^*$ (unitless)'));
            legend([tot res c1 c2], '$\Delta R_1^*$', 'Residual', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside', 'orientation', 'vertical');
        end
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
    elseif leg == "_noleg"
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    end
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');

    if strcmp(leg, "_legonly")
        axis off;
        axis([100 101 100 101])
        legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\partial_t h + \nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\partial_t h + \nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'northwest', 'numcolumns', 2, 'orientation', 'horizontal');
        set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 3.6 0.72])
    end
    
    print(plotname, '-dpng', '-r300');
    close;

end
