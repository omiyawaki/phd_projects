function plot_dr1_polar_line_topocomp(type, par)
    make_dirs(type, par)

    if strcmp(type, 'echam') & strcmp(par.echam.clim, 'echr0001')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        prefix_comp=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, 'echr0023');
        prefix_proc_comp=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, 'echr0023');
    else
        error('This function is only relevant for ECHAM AGCM simulation testing topography.');
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    tmp = load(sprintf('%s/%s/flux_z.mat', prefix_proc_comp, par.lat_interp)); flux_z_comp = tmp.flux_z; clear tmp; % load lat x mon RCAE data
    % load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [80];

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'L+O';
        elseif strcmp(land, 'l'); land_text = 'L';
        elseif strcmp(land, 'o'); land_text = 'O';
        end

        if strcmp(type, 'echam')
            if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                latn_pole = 90; latn = lat_bound:dlat:latn_pole; monlabel=par.monlabel;
                clatn = cosd(latn); % cosine of latnitude for cosine weighting
                clatn_mon = repmat(clatn', [1 12]);
                lats_pole = -90; lats = -lat_bound:-dlat:lats_pole; % monlabel=par.monlabelsh;
                clats = cosd(lats); % cosine of latsitude for cosine weighting
                clats_mon = repmat(clats', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_poleward_of_lat_-%g_topocomp', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                % fm and ra components
                fm_latn = interp1(grid.dim3.lat, flux_z.(land).res.(fw), latn);
                fm_latn = nansum(fm_latn.*clatn_mon)/nansum(clatn);
                fm_lats = interp1(grid.dim3.lat, flux_z.(land).res.(fw), lats);
                fm_lats = nansum(fm_lats.*clats_mon)/nansum(clats);
                comp_fm_latn = interp1(grid.dim3.lat, flux_z_comp.(land).res.(fw), latn);
                comp_fm_latn = nansum(comp_fm_latn.*clatn_mon)/nansum(clatn);
                comp_fm_lats = interp1(grid.dim3.lat, flux_z_comp.(land).res.(fw), lats);
                comp_fm_lats = nansum(comp_fm_lats.*clats_mon)/nansum(clats);

                ra_latn = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), latn);
                ra_latn = nansum(ra_latn.*clatn_mon)/nansum(clatn);
                ra_lats = interp1(grid.dim3.lat, flux_z.(land).ra.(fw), lats);
                ra_lats = nansum(ra_lats.*clats_mon)/nansum(clats);
                comp_ra_latn = interp1(grid.dim3.lat, flux_z_comp.(land).ra.(fw), latn);
                comp_ra_latn = nansum(comp_ra_latn.*clatn_mon)/nansum(clatn);
                comp_ra_lats = interp1(grid.dim3.lat, flux_z_comp.(land).ra.(fw), lats);
                comp_ra_lats = nansum(comp_ra_lats.*clats_mon)/nansum(clats);

                stf_latn = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), latn);
                stf_latn = nansum(stf_latn.*clatn_mon)/nansum(clatn);
                stf_lats = interp1(grid.dim3.lat, flux_z.(land).stf.(fw), lats);
                stf_lats = nansum(stf_lats.*clats_mon)/nansum(clats);
                comp_stf_latn = interp1(grid.dim3.lat, flux_z_comp.(land).stf.(fw), latn);
                comp_stf_latn = nansum(comp_stf_latn.*clatn_mon)/nansum(clatn);
                comp_stf_lats = interp1(grid.dim3.lat, flux_z_comp.(land).stf.(fw), lats);
                comp_stf_lats = nansum(comp_stf_lats.*clats_mon)/nansum(clats);

                [lh, sh] = rename_stf(type, flux_z, land);
                [comp_lh, comp_sh] = rename_stf(type, flux_z_comp, land);

                lh_latn = interp1(grid.dim3.lat, lh, latn);
                lh_latn = nansum(lh_latn.*clatn_mon)/nansum(clatn);
                lh_lats = interp1(grid.dim3.lat, lh, lats);
                lh_lats = nansum(lh_lats.*clats_mon)/nansum(clats);
                comp_lh_latn = interp1(grid.dim3.lat, comp_lh, latn);
                comp_lh_latn = nansum(comp_lh_latn.*clatn_mon)/nansum(clatn);
                comp_lh_lats = interp1(grid.dim3.lat, comp_lh, lats);
                comp_lh_lats = nansum(comp_lh_lats.*clats_mon)/nansum(clats);

                sh_latn = interp1(grid.dim3.lat, sh, latn);
                sh_latn = nansum(sh_latn.*clatn_mon)/nansum(clatn);
                sh_lats = interp1(grid.dim3.lat, sh, lats);
                sh_lats = nansum(sh_lats.*clats_mon)/nansum(clats);
                comp_sh_latn = interp1(grid.dim3.lat, comp_sh, latn);
                comp_sh_latn = nansum(comp_sh_latn.*clatn_mon)/nansum(clatn);
                comp_sh_lats = interp1(grid.dim3.lat, comp_sh, lats);
                comp_sh_lats = nansum(comp_sh_lats.*clats_mon)/nansum(clats);

                % REPLACE flat topo Ra with NH
                r1z_latn = fm_latn./ra_latn;
                r1z_lats = fm_lats./ra_lats;
                comp_r1z_latn = comp_fm_latn./comp_ra_latn;
                comp_r1z_lats = comp_fm_lats./comp_ra_lats;

                % r1z_lats_ra_n = (circshift(ra_latn,6,2) + lh_lats + sh_lats)./circshift(ra_latn,6,2);
                % r1z_lats_stf_n = (ra_lats + circshift(stf_latn,6,2))./ra_lats;
                % r1z_lats_lh_n = (ra_lats + circshift(lh_latn,6,2) + sh_lats)./ra_lats;
                % r1z_lats_sh_n = (ra_lats + circshift(sh_latn,6,2) + lh_lats)./ra_lats;

                r1z_lats_ra_f = (comp_ra_lats + lh_lats + sh_lats)./comp_ra_lats;
                r1z_lats_lh_f = (ra_lats + comp_lh_lats + sh_lats)./ra_lats;
                r1z_lats_sh_f = (ra_lats + lh_lats + comp_sh_lats)./ra_lats;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                ylim_lo = 0.2;
                ylim_up = 2;
                raemin = par.ga;
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                all_s=plot([1:12], circshift(r1z_lats,6,2), 'k');
                comp_all_s=plot([1:12], circshift(comp_r1z_lats,6,2), ':k');
                % all_n=plot([1:12], r1z_latn, '--k');
                % north=plot([1:12], r1z_latn, '--k');
                ra_f=plot([1:12], circshift(r1z_lats_ra_f,6,2), 'color', 0.5*[1 1 1]);
                % stf_n=plot([1:12], circshift(r1z_lats_stf_n,6,2), '-.k');
                lh_f=plot([1:12], circshift(r1z_lats_lh_f,6,2), 'color', par.blue);
                sh_f=plot([1:12], circshift(r1z_lats_sh_f,6,2), 'color', par.orange);
                xlabel('Month (relative to winter solstice)');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                title(sprintf('%s, AGCM, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound, lats_pole));
                % xlabel('Month');
                legend([all_s comp_all_s ra_f lh_f sh_f], '$R_1^C$', '$R_1^F$', '$\frac{R_a^F + \mathrm{NH}^C+\mathrm{SH}^C}{R_a^F}$', '$\frac{R_a^C+\mathrm{LH}^F+\mathrm{SH}^C}{R_a^C}$', '$\frac{R_a^C+\mathrm{LH}^C+\mathrm{SH}^F}{R_a^C}$', 'location', 'eastoutside', 'numcolumns', 1);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                print(sprintf('%s/0_mon_dr1z_repl', folder), '-dpng', '-r300');
                close;

                % decomp_ra = -(lh_lats+sh_lats)./ra_lats.^2.*(comp_ra_lats-ra_lats);
                % decomp_lh = (comp_lh_lats-lh_lats)./ra_lats;
                % decomp_sh = (comp_sh_lats-sh_lats)./ra_lats;
                decomp_ra = -(comp_lh_lats+comp_sh_lats)./comp_ra_lats.^2.*(ra_lats-comp_ra_lats);
                decomp_lh = (lh_lats-comp_lh_lats)./comp_ra_lats;
                decomp_sh = (sh_lats-comp_sh_lats)./comp_ra_lats;
                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                line([1,12],[0,0], 'linewidth', 0.5, 'color', 'k');
                % ylim_lo = 0.2;
                % ylim_up = 2;
                % raemin = par.ga;
                % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                dr1=plot([1:12], r1z_lats - comp_r1z_lats, 'k');
                res=plot([1:12], (r1z_lats-comp_r1z_lats)-(decomp_ra+decomp_lh+decomp_sh), '-.k');
                % plot([1:12], circshift(r1z_lats_alt,6,2), ':k');
                % % north=plot([1:12], r1z_latn, '--k');
                dra_n=plot([1:12], decomp_ra, 'color', 0.5*[1 1 1]);
                % % stf_n=plot([1:12], circshift(r1z_lats_stf_n,6,2), '-.k');
                dlh_n=plot([1:12], decomp_lh, 'color', par.blue);
                dsh_n=plot([1:12], decomp_sh, 'color', par.orange);
                xlabel('Month (relative to winter solstice)');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-0.2 0.6], 'yminortick', 'on', 'tickdir', 'out');
                title(sprintf('%s, AGCM, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound, lats_pole));
                % xlabel('Month');
                legend([dr1 res dra_n dlh_n dsh_n], '$R_1^C-R_1^F$', 'Residual', '$-\frac{\mathrm{LH}^F+\mathrm{SH}^F}{(R_a^F)^2}(R_a^C-R_a^F)$', '$\frac{1}{R_a^F}(\mathrm{LH}^C-\mathrm{LH}^F)$', '$\frac{1}{R_a^F}(\mathrm{SH}^C-\mathrm{SH}^F)$', 'location', 'eastoutside', 'numcolumns', 1);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                print(sprintf('%s/0_mon_dr1z_decomp', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
