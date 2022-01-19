function make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    for l = par.land_list; land = l{1};
        for plev_eval = [300:100:500]
            if ~exist(sprintf('%s/ma_diff/plev_%g/%s', plotdir, plev_eval, land), 'dir')
                mkdir(sprintf('%s/ma_diff/plev_%g/%s', plotdir, plev_eval, land));
            end
        end
        for si_eval = [0.8 0.85 0.9]
            if ~exist(sprintf('%s/inv_str/si_%g/%s', plotdir, si_eval, land), 'dir')
                mkdir(sprintf('%s/inv_str/si_%g/%s', plotdir, si_eval, land));
            end
        end
        for m = [1 4 6 7 8 9 10]; month = m(1);
            if ~exist(sprintf('%s/temp_zon_sel/%s/%s/%g', plotdir, num2str(par.ma_init), land, month), 'dir')
                mkdir(sprintf('%s/temp_zon_sel/%s/%s/%g', plotdir, num2str(par.ma_init), land, month));
            end
            if ~exist(sprintf('%s/temp_zon_sel_comp/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/temp_zon_sel_comp/%s/%g', plotdir, land, month));
            end
            if ~exist(sprintf('%s/thetaeq_zon_sel/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/thetaeq_zon_sel/%s/%g', plotdir, land, month));
            end
            if ~exist(sprintf('%s/mse_zon_sel/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/mse_zon_sel/%s/%g', plotdir, land, month));
            end
            if ~exist(sprintf('%s/dtdz_zon_sel/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/dtdz_zon_sel/%s/%g', plotdir, land, month));
            end
        end
        fw_vec = assign_fw(type, par);
        for f = fw_vec; fw = f{1};
            if ~exist(sprintf('%s/temp_binned_r1/%s/%s', plotdir, fw, land), 'dir')
                mkdir(sprintf('%s/temp_binned_r1/%s/%s', plotdir, fw, land));
            end
            if ~exist(sprintf('%s/ga_frac_binned_r1/%s/%s', plotdir, fw, land), 'dir')
                mkdir(sprintf('%s/ga_frac_binned_r1/%s/%s', plotdir, fw, land));
            end
            if ~exist(sprintf('%s/dtempsi_binned_r1/%s/%s', plotdir, fw, land), 'dir')
                mkdir(sprintf('%s/dtempsi_binned_r1/%s/%s', plotdir, fw, land));
            end
            if ~exist(sprintf('%s/eptlr_binned_r1/%s/%s', plotdir, fw, land), 'dir')
                mkdir(sprintf('%s/eptlr_binned_r1/%s/%s', plotdir, fw, land));
            end
            if ~exist(sprintf('%s/mselr_binned_r1/%s/%s', plotdir, fw, land), 'dir')
                mkdir(sprintf('%s/mselr_binned_r1/%s/%s', plotdir, fw, land));
            end
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/energy-flux/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux/%s/%s', plotdir, land, time));
            end
            if ~exist(sprintf('%s/energy-flux-comp/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux-comp/%s/%s', plotdir, land, time));
            end
            if ~exist(sprintf('%s/transport/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/transport/%s/%s', plotdir, land, time));
            end
            if ~exist(sprintf('%s/ga_diff/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/ga_diff/%s/%s', plotdir, land, time));
            end
            if par.do_surf; v_vec = {'p', 'z', 'si', 'pi'};
            else v_vec = {'p', 'z', 'si'}; end
            for v = v_vec; vert = v{1};
                if ~exist(sprintf('%s/temp_zon/%s/%s/%s', plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zon/%s/%s/%s', plotdir, land, time, vert));
                end
                if ~exist(sprintf('%s/temp_zonmean/%s/%s/%s', plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zonmean/%s/%s/%s', plotdir, land, time, vert));
                end
            end
            fw_vec = assign_fw(type, par);
            for f = fw_vec; fw = f{1};
                if ~exist(sprintf('%s/flux/%s/%s/%s', plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/flux/%s/%s/%s', plotdir, fw, land, time));
                elseif ~exist(sprintf('%s/dflux/%s/%s/%s', plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/dflux/%s/%s/%s', plotdir, fw, land, time));
                end
                if ~exist(sprintf('%s/flux-comp/%s/%s/%s', plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/flux-comp/%s/%s/%s', plotdir, fw, land, time));
                end
                if ~exist(sprintf('%s/dr2/%s/%s/%s', plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/dr2/%s/%s/%s', plotdir, fw, land, time));
                end
            end
            if ~exist(sprintf('%s/flag/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/flag/%s/%s', plotdir, land, time));
            end
        end
    end

    for vn = {'sice', 'olr_ts', 'divfm_lapt', 'va', 'trop', 'mid_mld', 'ga_lo', 'r1_lo', 'alb', 'sol', 'ps','tas', 'ts', 'sn', 'sftlf', 'legends', 'tend'}; varname = vn{1};
        if ~exist(sprintf('%s/%s', plotdir, varname), 'dir')
            mkdir(sprintf('%s/%s', plotdir, varname));
        end
    end
    
    if strcmp(type, 'echam')
        for vn = {'siced', 'friac', 'ahfres', 'ahfliac', 'ahfllac', 'ahflwac', 'ameltdepth', 'ameltfrac'}; varname = vn{1};
            if ~exist(sprintf('%s/%s', plotdir, varname), 'dir')
                mkdir(sprintf('%s/%s', plotdir, varname));
            end
        end
    end
end
