% plot showing the seasonal amplitude of r1 as a function of mixed layer depth

clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

echam_info
figure_params

par.echam_clims =  par.echam.noice_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)

% par.Q = 200; % seasonality of insolation W m^-2
par.Q = 150; % seasonality of insolation W m^-2
par.Ra = -100; % annual mean radiative cooling W m^-2
par.rho = 1000; % density of water kg m^-3
par.cw = 4000; % specific heat capacity of water J kg^-1 K^-1
par.omega = 2*pi/(86400*365);
par.A = 201.4; % OLR at T = 273 K from North (1975)
par.B = 1.45; % OLR vs T W m^-2 K^-1 from North (1975)
% par.B = 20; % OLR vs T W m^-2 K^-1 from North (1975)
% par.de = 0.31;
% par.de = 0.5;
par.de = 1.2;
par.a0 = 0.68;
par.Qg = 340;
par.be = 23;
par.Tf = -10; % reference temperature for nondimensionalized annual mean temperature
par.s11 = -2*sind(par.be);
par.s20 = -5/16*(2-3*sind(par.be)^2);
par.s22 = 15/16*sind(par.be)^2;
par.ttype = 'tsurf'; % use tsurf or temp2?

% mixed layer depths
d.rp000046 = 50;
d.rp000149 = 45;
d.rp000135 = 40;
d.rp000147 = 35;
d.rp000131 = 30;
d.rp000145 = 25;
d.rp000133 = 20;
d.rp000141 = 15;
d.rp000034 = 10;

type = 'echam';
par.lat_interp = 'native';

for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    if k==1; load(sprintf('%s/grid.mat', prefix)); end;
    tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
    flux_z.(par.echam.clim) = tmp.flux_z;
    tmp=load(sprintf('%s/srfc.mat', prefix)); % load lat x mon RCAE data
    srfc.(par.echam.clim) = tmp.srfc;
    tmp=load(sprintf('%s/rad.mat', prefix)); % load lat x mon RCAE data
    rad.(par.echam.clim) = tmp.rad;

    % plot zonal temperature structure

    ann_ts = squeeze(nanmean(nanmean(srfc.(par.echam.clim).(par.ttype),1),3));
    if k==1;
        p2 = 1/2 * (3*sind(grid.dim2.lat).^2-1);
        pred_ts_lat = par.a0*par.Qg/(par.A+par.B*par.Tf)*(1 + (par.s20)/(1+6*par.de).*p2); % non dimensional
        pred_ts_lat = ((par.A+par.B*par.Tf)*pred_ts_lat-par.A)/par.B+273; % dimensional
    end
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, ann_ts, '--k');
    plot(grid.dim2.lat, pred_ts_lat, '-k');
    xlabel('latitude (deg)'); ylabel('$T_s$ (K)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
    set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', [-90 90], 'xtick', [-90:30:90]);
    print(sprintf('./figures_post/test/amp_r1_echam/lat_ts_echam_%g', d.(par.echam.clim)), '-dpng', '-r300');
    close;


    lat_bound_list = [15];

    for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
        dlat = 0.25; % step size for standard lat grid
        if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center;
        else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        % interpolate to midlat
        tsz_lat.(par.echam.clim) = interp1(grid.dim2.lat, squeeze(nanmean(srfc.(par.echam.clim).(par.ttype), 1)), lat);
        tsz_lat.(par.echam.clim) = nansum(tsz_lat.(par.echam.clim).*clat_mon)/nansum(clat);

        % annual mean
        tsz_ann.(par.echam.clim) = repmat(squeeze(nanmean(nanmean(srfc.(par.echam.clim).(par.ttype),1),3))', [1 12]);
        tsz_ann_lat.(par.echam.clim) = interp1(grid.dim2.lat, tsz_ann.(par.echam.clim), lat);
        tsz_ann_lat.(par.echam.clim) = nansum(tsz_ann_lat.(par.echam.clim).*clat_mon)/nansum(clat);

        % seasonality
        dtsz.(par.echam.clim) = squeeze(nanmean(srfc.(par.echam.clim).(par.ttype),1)) - tsz_ann.(par.echam.clim);
        dtsz_lat.(par.echam.clim) = interp1(grid.dim2.lat, dtsz.(par.echam.clim), lat);
        dtsz_lat.(par.echam.clim) = nansum(dtsz_lat.(par.echam.clim).*clat_mon)/nansum(clat);

        % range of seasonality
        range_dtsz.(par.echam.clim) = nanmax(dtsz_lat.(par.echam.clim)) - nanmin(dtsz_lat.(par.echam.clim));

        for l = {'lo'}; land = l{1};
            for f = {'mse'}; fw = f{1};
                % interpolate to midlat
                r1z_lat.(par.echam.clim) = interp1(grid.dim3.lat, flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw), lat);
                r1z_lat.(par.echam.clim) = nansum(r1z_lat.(par.echam.clim).*clat_mon)/nansum(clat);

                % annual mean
                r1z_ann.(par.echam.clim) = repmat(nanmean(flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat.(par.echam.clim) = interp1(grid.dim3.lat, r1z_ann.(par.echam.clim), lat);
                r1z_ann_lat.(par.echam.clim) = nansum(r1z_ann_lat.(par.echam.clim).*clat_mon)/nansum(clat);

                % seasonality
                dr1z.(par.echam.clim) = flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw) - r1z_ann.(par.echam.clim);
                dr1z_lat.(par.echam.clim) = interp1(grid.dim3.lat, dr1z.(par.echam.clim), lat);
                dr1z_lat.(par.echam.clim) = nansum(dr1z_lat.(par.echam.clim).*clat_mon)/nansum(clat);

                % range of seasonality
                range_dr1.(par.echam.clim) = nanmax(dr1z_lat.(par.echam.clim)) - nanmin(dr1z_lat.(par.echam.clim));

                % minimum of dR1
                min_dr1.(par.echam.clim) = nanmin(dr1z_lat.(par.echam.clim));

            end % fw

        end % land

    end %lat list

end

d_vec = linspace(10,50,100);
% pred_dts = par.Q./(par.B*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
% pred_dts = nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.a0*par.s11*par.Qg./(par.B*(1+2*par.de)*sqrt(1+((par.rho*par.cw*d_vec*par.omega/par.B)./(1+2*par.de)).^2));
pred_dts = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.a0*par.s11*par.Qg./(par.B*sqrt((1+2*par.de).^2+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
figure(); clf; hold all; box on;
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), range_dtsz.(par.echam.clim)/2, '*k');
end
plot(d_vec, pred_dts, '-k');
xlabel('d (m)'); ylabel('$\Delta T_s$ (K)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on');
print('./figures_post/test/amp_r1_echam/amp_ts_echam', '-dpng', '-r300');
close;

% pred_min_dr1 = par.Q./(par.Ra*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
pred_min_dr1 = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2*par.de./(par.Ra*sqrt((1+2*par.de).^2+(par.rho*par.cw*d_vec*par.omega/par.B).^2));

figure(); clf; hold all; box on;
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), min_dr1.(par.echam.clim), '*k');
end
plot(d_vec, pred_min_dr1, '-k');
xlabel('d (m)'); ylabel('$\min(\Delta R_1)$ (unitless)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on');
print('./figures_post/test/amp_r1_echam/amp_r1_echam', '-dpng', '-r300');
close;
