function plot_temp_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    if par.ma
        if strcmp(par.ma_init, 'surf')
            load(sprintf('%s/ma_mon_lat_%s.mat', prefix_proc, par.ma_init));
            subplotdir = sprintf('temp_zon_sel/%s', par.ma_init);
        else
            load(sprintf('%s/ma_mon_lat_%g.mat', prefix_proc, par.ma_init));
            subplotdir = sprintf('temp_zon_sel/%g', par.ma_init);
        end
    end

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 4 6 7 10]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==4; mon_str = 'April';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July';
            elseif month==10; mon_str = 'October'; end;

            tasi_mon.(land) = squeeze(tasi.(land)(:,month,:));
            if par.ma
                masi_mon.(land) = squeeze(masi.(land)(:,month,:));

                % remove moist adiabat data below initialization level
                if ~strcmp(par.ma_init, 'surf')
                    masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
                end
            end

            tasi_sp(:,m) = interp1(lat, tasi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tasi_np(:,m) = interp1(lat, tasi_mon.(land), lat_pole); % sounding at lat_pole N
            tasi_smid(:,m) = interp1(lat, tasi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tasi_nmid(:,m) = interp1(lat, tasi_mon.(land), lat_mid); % sounding at lat_mid N
            tasi_eq(:,m) = interp1(lat, tasi_mon.(land), 0); % sounding at equator

            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                tasi_std_mon.(land) = squeeze(tasi_std.(land)(:,month,:));
                tasi_std_sp(:,m) = interp1(lat, tasi_std_mon.(land), -lat_pole); % sounding at -lat_pole S
                tasi_std_np(:,m) = interp1(lat, tasi_std_mon.(land), lat_pole); % sounding at lat_pole N
                tasi_std_smid(:,m) = interp1(lat, tasi_std_mon.(land), -lat_mid); % sounding at -lat_mid S
                tasi_std_nmid(:,m) = interp1(lat, tasi_std_mon.(land), lat_mid); % sounding at lat_mid N
                tasi_std_eq(:,m) = interp1(lat, tasi_std_mon.(land), 0); % sounding at equator
            end

            if par.ma
                masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
                masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
                masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
                masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
                masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator
            end

            % ALL
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                end
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
                end
            end
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, '--', 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, '--', 'color', 0.25*[1 1 1]);
            if par.ma
                h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            %if par.ma
            %    legend([h_np h_sp h_nmid h_smid], sprintf('%g N', lat_pole), sprintf('%g S', lat_pole), sprintf('%g N', lat_mid), sprintf('%g S', lat_mid), 'location', 'northeast');
            %end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/winter_summer', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;
            
            % winter summer 
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                end
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
                end
            end
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, '--', 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, '--', 'color', 0.25*[1 1 1]);
            if par.ma
                h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            %if par.ma
            %    legend([h_np h_sp h_nmid h_smid], sprintf('%g N', lat_pole), sprintf('%g S', lat_pole), sprintf('%g N', lat_mid), sprintf('%g S', lat_mid), 'location', 'northeast');
            %end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/winter_summer', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/np', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                end
            elseif m==6 | m == 7
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
                end
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, lat_mid, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/nhmid', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                tasi_np_u = tasi_np(:,m) + tasi_std_np(:,m);
                tasi_np_l = tasi_np(:,m) - tasi_std_np(:,m);
                tasi_nmid_u = tasi_nmid(:,m) + tasi_std_nmid(:,m);
                tasi_nmid_l = tasi_nmid(:,m) - tasi_std_nmid(:,m);
                
                tasi_np_2 = [tasi_np_l; flipud(tasi_np_u)];
                tasi_nmid_2 = [tasi_nmid_l; flipud(tasi_nmid_u)];
                si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                if m == 1
                    fill(tasi_np_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(tasi_nmid_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                elseif m==6 | m == 7
                    fill(tasi_np_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                    fill(tasi_nmid_2, si2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
                end
            end
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                text(tasi_np(30,m)-15, grid.dim3.si(30), sprintf('$\\mathbf{%g^\\circ}$\\textbf{N}', lat_pole), 'color', par.blue);
                text(tasi_nmid(50,m)+5, grid.dim3.si(50), sprintf('$\\mathbf{%g^\\circ}$\\textbf{N}', lat_mid), 'color', 0.25*[1 1 1]);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                end
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                text(tasi_np(30,m)-15, grid.dim3.si(30), sprintf('$\\mathbf{%g^\\circ}$\\textbf{N}', lat_pole), 'color', 0.25*[1 1 1]);
                text(tasi_nmid(50,m)+5, grid.dim3.si(50), sprintf('$\\mathbf{%g^\\circ}$\\textbf{N}', lat_mid), 'color', par.orange);
                if par.ma
                    h_np_ma = plot(masi_np(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
                end
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            %if par.ma
            %    legend([h_np h_nmid], sprintf('%g N', lat_pole), sprintf('%g N', lat_mid), 'location', 'northeast');
            %end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/nh_only', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % SH HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/sp', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID ONLY
            figure(); clf; hold all; box on;
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            if par.ma
                h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, -lat_mid, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/shmid', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % EQUATOR ONLY
            figure(); clf; hold all; box on;
            h_eq = plot(tasi_eq(:,m), grid.dim3.si, 'color', par.orange);
            if par.ma
                h_eq_ma = plot(masi_eq(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon_lat(type, mon_str, 0, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_eq(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/eq', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                tasi_sp_u = tasi_sp(:,m) + tasi_std_sp(:,m);
                tasi_sp_l = tasi_sp(:,m) - tasi_std_sp(:,m);
                tasi_smid_u = tasi_smid(:,m) + tasi_std_smid(:,m);
                tasi_smid_l = tasi_smid(:,m) - tasi_std_smid(:,m);
                
                tasi_sp_2 = [tasi_sp_l; flipud(tasi_sp_u)];
                tasi_smid_2 = [tasi_smid_l; flipud(tasi_smid_u)];
                si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                fill(tasi_sp_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                fill(tasi_smid_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
            end
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            text(tasi_sp(30,m)-15, grid.dim3.si(30), sprintf('$\\mathbf{%g^\\circ}$\\textbf{S}', lat_pole), 'color', par.blue);
            text(tasi_smid(50,m)+5, grid.dim3.si(50), sprintf('$\\mathbf{%g^\\circ}$\\textbf{S}', lat_mid), 'color', 0.25*[1 1 1]);
            if par.ma
                h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            %if par.ma
            %    legend([h_sp h_smid], sprintf('%g S', lat_pole), sprintf('%g S', lat_mid), 'location', 'northeast');
            %end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/sh_only', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                end
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                if par.ma
                    h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
                end
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/nh_only_vert', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            if par.ma
                h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/%s/%s/%g/sh_only_vert', plotdir, subplotdir, land, month), '-dpng', '-r300');
            close;

        end

        % NH HL
        figure(); clf; hold all; box on;
        h_np_wi = plot(tasi_np(:,1), grid.dim3.si, 'color', par.gray);
        h_np_sp = plot(tasi_np(:,4), grid.dim3.si, 'color', par.yellow);
        h_np_su = plot(tasi_np(:,7), grid.dim3.si, 'color', par.green);
        h_np_fa = plot(tasi_np(:,10), grid.dim3.si, 'color', par.orange);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_pole, par);
        axis('tight');
        legend([h_np_wi, h_np_sp, h_np_su, h_np_fa], 'Jan', 'Apr', 'Jul', 'Oct', 'location', 'eastoutside', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_hl', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % NH HL Jan Jun
        if strcmp(type, 'echam') & strcmp(par.echam.clim, 'rp000135')
            color_wi = par.gray;
        else
            color_wi = par.blue;
        end
        color_su = par.gray;
        figure(); clf; hold all; box on;
        h_np_wi = plot(tasi_np(:,1), grid.dim3.si, 'color', color_wi);
        h_np_su = plot(tasi_np(:,6), grid.dim3.si, 'color', color_su);
        text(tasi_np(30,1)-20, grid.dim3.si(30), '\textbf{January}', 'color', color_wi);
        text(tasi_np(50,6)+5, grid.dim3.si(50), '\textbf{June}', 'color', color_su);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_pole, par);
        axis('tight');
        % legend([h_np_wi, h_np_su], 'Jan', 'Jun', 'location', 'southwest', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_hl_1_6', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % NH HL Jan Jul
        if strcmp(type, 'echam') & strcmp(par.echam.clim, 'rp000135')
            color_wi = par.gray;
        else
            color_wi = par.blue;
        end
        color_su = par.gray;
        figure(); clf; hold all; box on;
        h_np_wi = plot(tasi_np(:,1), grid.dim3.si, 'color', color_wi);
        h_np_su = plot(tasi_np(:,7), grid.dim3.si, 'color', color_su);
        text(tasi_np(30,1)-20, grid.dim3.si(30), '\textbf{January}', 'color', color_wi);
        text(tasi_np(50,7)+5, grid.dim3.si(50), '\textbf{July}', 'color', color_su);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_pole, par);
        axis('tight');
        % legend([h_np_wi, h_np_su], 'Jan', 'Jun', 'location', 'southwest', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_hl_1_7', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % NH ML
        figure(); clf; hold all; box on;
        h_nmid_wi = plot(tasi_nmid(:,1), grid.dim3.si, 'color', par.gray);
        h_nmid_sp = plot(tasi_nmid(:,4), grid.dim3.si, 'color', par.yellow);
        h_nmid_su = plot(tasi_nmid(:,7), grid.dim3.si, 'color', par.green);
        h_nmid_fa = plot(tasi_nmid(:,10), grid.dim3.si, 'color', par.orange);
        h_nmid_wi_ma = plot(masi_nmid(:,1), grid.dim3.si, ':', 'color', par.gray);
        h_nmid_sp_ma = plot(masi_nmid(:,4), grid.dim3.si, ':', 'color', par.yellow);
        h_nmid_su_ma = plot(masi_nmid(:,7), grid.dim3.si, ':', 'color', par.green);
        h_nmid_fa_ma = plot(masi_nmid(:,10), grid.dim3.si, ':', 'color', par.orange);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_mid, par);
        axis('tight');
        legend([h_nmid_wi, h_nmid_sp, h_nmid_su, h_nmid_fa], 'Jan', 'Apr', 'Jul', 'Oct', 'location', 'eastoutside', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_ml', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % NH ML Jan Jun
        if strcmp(type, 'echam') & strcmp(par.echam.clim, 'rp000135')
            color_su = par.gray;
        else
            color_su = par.orange;
        end
        color_wi = par.gray;
        figure(); clf; hold all; box on;
        h_nmid_wi = plot(tasi_nmid(:,1), grid.dim3.si, 'color', color_wi);
        h_nmid_su = plot(tasi_nmid(:,6), grid.dim3.si, 'color', color_su);
        h_nmid_wi_ma = plot(masi_nmid(:,1), grid.dim3.si, ':', 'color', color_wi);
        h_nmid_su_ma = plot(masi_nmid(:,6), grid.dim3.si, ':', 'color', color_su);
        text(tasi_nmid(30,1)-30, grid.dim3.si(30), '\textbf{January}', 'color', color_wi);
        text(tasi_nmid(50,6)+5, grid.dim3.si(50), '\textbf{June}', 'color', color_su);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_mid, par);
        axis('tight');
        % legend([h_nmid_wi, h_nmid_su], 'Jan', 'Jun', 'location', 'southwest', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_ml_1_6', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % NH ML Jan Jul
        if strcmp(type, 'echam') & strcmp(par.echam.clim, 'rp000135')
            color_su = par.gray;
        else
            color_su = par.orange;
        end
        color_wi = par.gray;
        figure(); clf; hold all; box on;
        h_nmid_wi = plot(tasi_nmid(:,1), grid.dim3.si, 'color', color_wi);
        h_nmid_su = plot(tasi_nmid(:,7), grid.dim3.si, 'color', color_su);
        h_nmid_wi_ma = plot(masi_nmid(:,1), grid.dim3.si, ':', 'color', color_wi);
        h_nmid_su_ma = plot(masi_nmid(:,7), grid.dim3.si, ':', 'color', color_su);
        text(tasi_nmid(30,1)-30, grid.dim3.si(30), '\textbf{January}', 'color', color_wi);
        text(tasi_nmid(50,7)+5, grid.dim3.si(50), '\textbf{July}', 'color', color_su);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, lat_mid, par);
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/nh_ml_1_7', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

        % SH HL
        figure(); clf; hold all; box on;
        h_sp_wi = plot(tasi_sp(:,7), grid.dim3.si, 'color', par.gray);
        h_sp_sp = plot(tasi_sp(:,10), grid.dim3.si, 'color', par.yellow);
        h_sp_su = plot(tasi_sp(:,1), grid.dim3.si, 'color', par.green);
        h_sp_fa = plot(tasi_sp(:,4), grid.dim3.si, 'color', par.orange);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, -lat_pole, par);
        axis('tight');
        legend([h_sp_wi, h_sp_sp, h_sp_su, h_sp_fa], 'Jul', 'Oct', 'Jan', 'Apr', 'location', 'eastoutside', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/sh_hl', plotdir, subplotdir, land), '-dpng', '-r300');
        close;
        
        % SH ML
        figure(); clf; hold all; box on;
        h_smid_wi = plot(tasi_smid(:,7), grid.dim3.si, 'color', par.gray);
        h_smid_sp = plot(tasi_smid(:,10), grid.dim3.si, 'color', par.yellow);
        h_smid_su = plot(tasi_smid(:,1), grid.dim3.si, 'color', par.green);
        h_smid_fa = plot(tasi_smid(:,4), grid.dim3.si, 'color', par.orange);
        h_smid_wi_ma = plot(masi_smid(:,7), grid.dim3.si, 'color', par.gray);
        h_smid_sp_ma = plot(masi_smid(:,10), grid.dim3.si, 'color', par.yellow);
        h_smid_su_ma = plot(masi_smid(:,1), grid.dim3.si, 'color', par.green);
        h_smid_fa_ma = plot(masi_smid(:,4), grid.dim3.si, 'color', par.orange);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        make_title_type_lat_pt(type, -lat_mid, par);
        axis('tight');
        legend([h_smid_wi, h_smid_sp, h_smid_su, h_smid_fa], 'Jul', 'Oct', 'Jan', 'Apr', 'location', 'eastoutside', 'orientation', 'vertical');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/%s/%s/sh_ml', plotdir, subplotdir, land), '-dpng', '-r300');
        close;

    end
end
