function make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    for l = {'lo', 'l', 'o'}; land = l{1};
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
        for m = [1 6 7]; month = m(1);
            if ~exist(sprintf('%s/temp_zon_sel/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/temp_zon_sel/%s/%g', plotdir, land, month));
            end
            if ~exist(sprintf('%s/thetaeq_zon_sel/%s/%g', plotdir, land, month), 'dir')
                mkdir(sprintf('%s/thetaeq_zon_sel/%s/%g', plotdir, land, month));
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
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/energy-flux/%s/%s', plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux/%s/%s', plotdir, land, time));
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

    for vn = {'va', 'trop', 'alb', 'sol', 'tas', 'ts', 'sn', 'siced', 'friac', 'sftlf', 'legends'}; varname = vn{1};
        if ~exist(sprintf('%s/%s', plotdir, varname), 'dir')
            mkdir(sprintf('%s/%s', plotdir, varname));
        end
    end
end
