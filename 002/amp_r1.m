% plot showing the seasonal amplitude of r1 as a function of mixed layer depth

clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));
addpath(genpath('./matlab'));

echam_info
figure_params

par.echam_clims =  par.echam.noice_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)

par.ep = 0.1;

% par.Q = 200; % seasonality of insolation W m^-2
par.Q = 150; % seasonality of insolation W m^-2
par.Ra = -100; % annual mean radiative cooling W m^-2
% par.Ra = -110; % annual mean radiative cooling W m^-2
par.rho = 1000; % density of water kg m^-3
% par.rho = 997; % density of water kg m^-3
par.cw = 4000; % specific heat capacity of water J kg^-1 K^-1
% par.cw = 4179.6; % specific heat capacity of water J kg^-1 K^-1
par.omega = 2*pi/(86400*365);

% par.de = 0.31;
% par.de = 0.5;

% par.choice = "north";
% par.B = 1.45; % OLR vs T W m^-2 K^-1 from North (1975)
% par.A = 201.4; % OLR at T = 273 K from North (1975)
% par.de = 0.31; % Rose (2017) /North (1975)
% par.D = par.B*par.de;
% par.a0 = 0.68; % North (1975)

% par.choice = "koll";
% par.B = 2.218; % OLR vs T W m^-2 K^-1 from Koll & Cronin (2018)
% par.A = -339.647+par.B*273.15; % KC18
% par.de = 0.31; % Rose (2017) /North (1975)
% par.D = par.B*par.de;
% par.a0 = 0.68; % North (1975)

par.choice = "echam";
par.B = 2.33; % best fit to ECHAM 25 m
par.A = -409.75+par.B*273.15; % best fit to ECHAM 25 m
par.D = 0.90; % best fit to ECHAM 25 m
par.de = par.D/par.B; % best fit to ECHAM 25 m
par.a0 = 0.72; % best fit to ECHAM 25 m

% par.choice = "tuned";
% par.B = 2.33; % best fit to ECHAM 25 m
% par.A = -409.75+par.B*273.15; % best fit to ECHAM 25 m
% par.D = 2; % best fit to ECHAM 25 m
% par.de = par.D/par.B; % best fit to ECHAM 25 m
% par.a0 = 0.72; % best fit to ECHAM 25 m

% par.a0 = 0.68;
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
d.rp000086 = 5;
d.rp000172 = 3;

% best fit ECHAM coefficients for each simulation
B.rp000046 = 2.33;
B.rp000149 = 2.34;
B.rp000135 = 2.33;
B.rp000147 = 2.34;
B.rp000131 = 2.34;
B.rp000145 = 2.33;
B.rp000133 = 2.31;
B.rp000141 = 2.27;
B.rp000034 = 2.15;
B.rp000086 = 1.97;
B.rp000172 = 1.87;

D.rp000046 = 0.87;
D.rp000149 = 0.88;
D.rp000135 = 0.88;
D.rp000147 = 0.90;
D.rp000131 = 0.90;
D.rp000145 = 0.90;
D.rp000133 = 0.90;
D.rp000141 = 0.89;
D.rp000034 = 0.66;
D.rp000086 = 0.49;
D.rp000172 = 0.36;

xlim = [0 60];

%% B vs d
d_vec = linspace(0,60,100);
% exponential fit
darr = struct2array(d);
Barr = struct2array(B);
f_B = @(bcoef,darr) bcoef(1) .* exp(bcoef(2).*darr) + bcoef(3);
c_B = fminsearch(@(bcoef) norm(Barr - f_B(bcoef,darr)), [-1; -1; 1]);
figure(); clf; hold all; box on;
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), B.(par.echam.clim), '*k');
end
plot(d_vec, f_B(c_B, d_vec), '-k');
xlabel('$d$ (m)'); ylabel('$B$ (W m$^{-2}$ K$^{-1}$)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim)%, 'ylim', [0 50]);
print(sprintf('./figures_post/test/amp_r1_echam/B_echam_%s', par.choice), '-dpng', '-r300');
close;

%% D vs d
% exponential fit
darr = struct2array(d);
Darr = struct2array(D);
f_D = @(dcoef,darr) dcoef(1) .* exp(dcoef(2).*darr) + dcoef(3);
c_D = fminsearch(@(dcoef) norm(Darr - f_B(dcoef,darr)), [-1; -1; 1]);
figure(); clf; hold all; box on;
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), D.(par.echam.clim), '*k');
end
plot(d_vec, f_D(c_D, d_vec), '-k');
xlabel('$d$ (m)'); ylabel('$D$ (W m$^{-2}$ K$^{-1}$)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim)%, 'ylim', [0 50]);
print(sprintf('./figures_post/test/amp_r1_echam/D_echam_%s', par.choice), '-dpng', '-r300');
close;

% %% Ra vs d
% % exponential fit
% darr = struct2array(d);
% Barr = struct2array(B);
% f_B = @(bcoef,darr) bcoef(1) .* exp(bcoef(2).*darr) + bcoef(3);
% c_B = fminsearch(@(bcoef) norm(Barr - f_B(bcoef,darr)), [-1; -1; 1]);
% figure(); clf; hold all; box on;
% for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%     plot(d.(par.echam.clim), B.(par.echam.clim), '*k');
% end
% plot(d_vec, f_B(c_B, d_vec), '-k');
% xlabel('$d$ (m)'); ylabel('$B$ (W m$^{-2}$ K$^{-1}$)');
% set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
% set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', [3 50])%, 'ylim', [0 50]);
% print(sprintf('./figures_post/test/amp_r1_echam/B_echam_%s', par.choice), '-dpng', '-r300');
% close;

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

    % compute EBM prediction of div Fm (D * laplacian of ts)
    clat = cosd(grid.dim2.lat);
    dphi = deg2rad(grid.dim2.lat(2)-grid.dim2.lat(1)); % uniform lat grid
    dtdlat = nan(size(pred_ts_lat));
    dtdlat(2:end-1) = (pred_ts_lat(3:end) - pred_ts_lat(1:end-2))/(2*dphi); % central difference
    dtdlat(1) = (pred_ts_lat(2) - pred_ts_lat(1))/dphi; % forward diff
    dtdlat(end) = (pred_ts_lat(end) - pred_ts_lat(end-1))/dphi; % backward diff
    dtdlat = clat .* dtdlat;
    lapt = nan(size(pred_ts_lat));
    lapt(2:end-1) = (dtdlat(3:end) - dtdlat(1:end-2))/(2*dphi); % central difference
    lapt(1) = (dtdlat(2) - dtdlat(1))/dphi; % forward diff
    lapt(end) = (dtdlat(end) - dtdlat(end-1))/dphi; % backward diff
    lapt = lapt./clat;
    pred_r1 = -par.D*lapt/par.Ra;

    % globally averaged values
    rsdt = squeeze(nanmean(rad.(par.echam.clim).srad0d, 1));
    swtoa = flux_z.(par.echam.clim).lo.sw - flux_z.(par.echam.clim).lo.swsfc;
    swtoa_ann = nanmean(swtoa, 2);
    rsdt_ann = nanmean(rsdt, 2);
    swtoa_glb = nansum(swtoa_ann .* clat)/nansum(clat);
    rsdt_glb = nansum(rsdt_ann .* clat)/nansum(clat);
    % disp(sprintf('a for d=%g m is %g', d.(par.echam.clim) ,swtoa_glb/rsdt_glb))

    % average over midlatitude region
    lat_bound_list = [10];
    center = 50;

    for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
        dlat = 0.25; % step size for standard lat grid
        if lat_bound>0; lat_center=center; lat = [-lat_bound:dlat:lat_bound]+lat_center;
        else; lat_center=-center; lat = [-lat_bound:-dlat:lat_bound]+lat_center; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        % predicted r1
        pred_r1_avg = interp1(grid.dim2.lat, pred_r1, lat);
        pred_r1_avg = nansum(clat.*pred_r1_avg)/nansum(clat);
    
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
    
        % predicted temperature seasonality
        t_vec = [0.5:11.5]*365/12*86400; % monthly steps expressed in seconds
        par.ga = par.rho*par.cw*d.(par.echam.clim)*par.omega/par.B;
        pred_ts_mon = par.a0*par.s11*par.Qg*nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))/par.B*((1+2*par.de).^2+par.ga.^2).^(-1/2)*cos(par.omega*t_vec-atan(par.ga./(1+2*par.de)));
        
        % plot seasonality of temperature
        figure(); clf; hold all; box on;
        plot([1:12], dtsz_lat.(par.echam.clim), '--k');
        plot([1:12], pred_ts_mon, '-k');
        xlabel('latitude (deg)'); ylabel('$T_s$ (K)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
        set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', par.monlabelnh);
        print(sprintf('./figures_post/test/amp_r1_echam/mon_ts_echam_%g', d.(par.echam.clim)), '-dpng', '-r300');
        close;
    
        % range of seasonality
        range_dtsz.(par.echam.clim) = nanmax(dtsz_lat.(par.echam.clim)) - nanmin(dtsz_lat.(par.echam.clim));
    
        for l = {'lo'}; land = l{1};
            for f = {'mse_old'}; fw = f{1};
                % interpolate to midlat
                r1z_lat.(par.echam.clim) = interp1(grid.dim3.lat, flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw), lat);
                r1z_lat.(par.echam.clim) = nansum(r1z_lat.(par.echam.clim).*clat_mon)/nansum(clat);
                raz_lat.(par.echam.clim) = interp1(grid.dim3.lat, flux_z.(par.echam.clim).(land).ra.(fw), lat);
                raz_lat.(par.echam.clim) = nansum(raz_lat.(par.echam.clim).*clat_mon)/nansum(clat);
    
                % annual mean
                r1z_ann.(par.echam.clim) = repmat(nanmean(flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat.(par.echam.clim) = interp1(grid.dim3.lat, r1z_ann.(par.echam.clim), lat);
                r1z_ann_lat.(par.echam.clim) = nansum(r1z_ann_lat.(par.echam.clim).*clat_mon)/nansum(clat);
                raz_ann.(par.echam.clim) = repmat(nanmean(flux_z.(par.echam.clim).(land).ra.(fw), 2), [1 12]);
                raz_ann_lat.(par.echam.clim) = interp1(grid.dim3.lat, raz_ann.(par.echam.clim), lat);
                raz_ann_lat.(par.echam.clim) = nansum(raz_ann_lat.(par.echam.clim).*clat_mon)/nansum(clat);
                disp(sprintf('Ra for d=%g m is %g', d.(par.echam.clim), raz_ann_lat.(par.echam.clim)(1)))
    
                % seasonality
                dr1z.(par.echam.clim) = flux_z.(par.echam.clim).(land).res.(fw)./flux_z.(par.echam.clim).(land).ra.(fw) - r1z_ann.(par.echam.clim);
                dr1z_lat.(par.echam.clim) = interp1(grid.dim3.lat, dr1z.(par.echam.clim), lat);
                dr1z_lat.(par.echam.clim) = nansum(dr1z_lat.(par.echam.clim).*clat_mon)/nansum(clat);

                % plot r1 seasonality with EBM prediction
                t_vec = (15:30.5:365)*86400;
                pred_r1_lat.(par.echam.clim) = pred_r1_avg + nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2*par.D./(par.Ra*((par.B+2*par.D).^2+(par.rho*par.cw*d.(par.echam.clim)*par.omega).^2)).*( (par.B + 2*par.D )*cos(par.omega * t_vec) + par.rho*par.cw*d.(par.echam.clim)*par.omega*sin(par.omega*t_vec) );
                figure(); clf; hold all; box on;
                plot([1:12], r1z_lat.(par.echam.clim), 'k');
                plot([1:12], pred_r1_lat.(par.echam.clim), ':k');
                set(gca, 'xlim', [1 12]);
                print(sprintf('./figures_post/test/amp_r1_echam/r1_echam_%s_d%gm', par.choice, d.(par.echam.clim)), '-dpng', '-r300');

                % range of seasonality
                range_dr1.(par.echam.clim) = nanmax(dr1z_lat.(par.echam.clim)) - nanmin(dr1z_lat.(par.echam.clim));
    
                % minimum of R1
                min_r1.(par.echam.clim) = nanmin(r1z_lat.(par.echam.clim));

                % minimum of dR1
                min_dr1.(par.echam.clim) = nanmin(dr1z_lat.(par.echam.clim));
    
            end % fw
    
        end % land
    
    end %lat list

end

%% Annual mean R1
figure(); clf; hold all; box on;
line([xlim(1) xlim(2)], pred_r1_avg*[1 1], 'color', 'k');
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), r1z_ann_lat.(par.echam.clim)(1), '*k');
end
xlabel('$d$ (m)'); ylabel('$\overline{R_1}$ (unitless)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim);
print(sprintf('./figures_post/test/amp_r1_echam/r1_ann_echam_%s', par.choice), '-dpng', '-r300');
close;

%% T* vs d
% pred_dts = par.Q./(par.B*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
% pred_dts = nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.a0*par.s11*par.Qg./(par.B*(1+2*par.de)*sqrt(1+((par.rho*par.cw*d_vec*par.omega/par.B)./(1+2*par.de)).^2));
pred_dts = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.a0*par.s11*par.Qg./(par.B*sqrt((1+2*par.de).^2+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
figure(); clf; hold all; box on;
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    plot(d.(par.echam.clim), range_dtsz.(par.echam.clim)/2, '*k');
end
plot(d_vec, pred_dts, '-k');
xlabel('$d$ (m)'); ylabel('$T_s^*$ (K)');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim, 'ylim', [0 50]);
print(sprintf('./figures_post/test/amp_r1_echam/amp_ts_echam_%s', par.choice), '-dpng', '-r300');
close;

%% min(R1) vs d
% pred_min_dr1 = par.Q./(par.Ra*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
% pred_min_dr1 = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2*par.de./(par.Ra*sqrt((1+2*par.de).^2+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
pred_min_dr1 = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2*par.D./(par.Ra*sqrt((par.B+2*par.D).^2+(par.rho*par.cw*d_vec*par.omega).^2));

figure(); clf; hold all; box on;
rcemax = par.ep; 
ymin = -1.1;
ymax = 0.5;
vertices = [xlim(1) ymin; xlim(2) ymin; xlim(2) rcemax; xlim(1) rcemax];
patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_r1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_r1_avg + pred_min_dr1, '-k');
xlabel('$d$ (m)'); ylabel('$\min(R_1)$ (unitless)');
legend([pecham, ppred], 'AQUA', 'EBM', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim, 'ylim', [ymin ymax]);
print(sprintf('./figures_post/test/amp_r1_echam/amp_r1_echam_%s', par.choice), '-dpng', '-r300');
close;

figure(); clf; hold all; box on;
% line([3, 50], [-0.2 -0.2], 'color', 'r');
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_dr1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_min_dr1, '-k');
xlabel('$d$ (m)'); ylabel('$\min(\Delta R_1)$ (unitless)');
legend([pecham, ppred], 'AQUA', 'EBM', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim);
print(sprintf('./figures_post/test/amp_r1_echam/amp_dr1_echam_%s', par.choice), '-dpng', '-r300');
close;

%% Empirically-derived D(d)
%% min(R1) vs d
% pred_min_dr1 = par.Q./(par.Ra*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
pred_min_dr1_D = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2.*f_D(c_D, d_vec)./(par.Ra*sqrt((par.B+2.*f_D(c_D, d_vec)).^2+(par.rho*par.cw*d_vec*par.omega).^2));

figure(); clf; hold all; box on;
rcemax = par.ep; 
ymin = -1.1;
ymax = 0.5;
vertices = [3 ymin; 50 ymin; 50 rcemax; 3 rcemax];
patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_r1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_r1_avg + pred_min_dr1_D, '-k');
xlabel('$d$ (m)'); ylabel('$\min(R_1)$ (unitless)');
legend([pecham, ppred], 'AQUA', 'EBM', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim, 'ylim', [ymin ymax]);
print(sprintf('./figures_post/test/amp_r1_echam/amp_r1_D(d)_echam_%s', par.choice), '-dpng', '-r300');
close;

%% Empirically-derived B(d)
%% min(R1) vs d
% pred_min_dr1 = par.Q./(par.Ra*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
pred_min_dr1_B = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2*par.D./(par.Ra*sqrt((f_B(c_B, d_vec)+2*par.D).^2+(par.rho*par.cw*d_vec*par.omega).^2));

figure(); clf; hold all; box on;
rcemax = par.ep; 
ymin = -1.1;
ymax = 0.5;
vertices = [xlim(1) ymin; xlim(2) ymin; xlim(2) rcemax; xlim(1) rcemax];
patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_r1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_r1_avg + pred_min_dr1_B, '-k');
xlabel('$d$ (m)'); ylabel('$\min(R_1)$ (unitless)');
legend([pecham, ppred], 'AQUA', 'EBM', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim, 'ylim', [ymin ymax]);
print(sprintf('./figures_post/test/amp_r1_echam/amp_r1_B(d)_echam_%s', par.choice), '-dpng', '-r300');
close;

figure(); clf; hold all; box on;
% line([3, 50], [-0.2 -0.2], 'color', 'r');
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_dr1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_min_dr1, '-k');
ppred_B = plot(d_vec, pred_min_dr1_B, ':k');
xlabel('$d$ (m)'); ylabel('$\min(\Delta R_1)$ (unitless)');
legend([pecham, ppred, ppred_B], 'AQUA', 'EBM', 'EBM B(d)', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim);
print(sprintf('./figures_post/test/amp_r1_echam/amp_dr1_B(d)_echam_%s', par.choice), '-dpng', '-r300');
close;

%% Empirically-derived D(d) and B(d)
%% min(R1) vs d
% pred_min_dr1 = par.Q./(par.Ra*sqrt(1+(par.rho*par.cw*d_vec*par.omega/par.B).^2));
pred_min_dr1_DB = - nansum(cosd(lat).*sind(lat))/nansum(cosd(lat))*par.s11*par.Qg*par.a0*2.*f_D(c_D, d_vec)./(par.Ra*sqrt((f_B(c_B, d_vec)+2.*f_D(c_D, d_vec)).^2+(par.rho*par.cw*d_vec*par.omega).^2));

figure(); clf; hold all; box on;
rcemax = par.ep; 
ymin = -1.1;
ymax = 0.5;
vertices = [xlim(1) ymin; xlim(2) ymin; xlim(2) rcemax; xlim(1) rcemax];
patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    pecham = plot(d.(par.echam.clim), min_r1.(par.echam.clim), '*k');
end
ppred = plot(d_vec, pred_r1_avg + pred_min_dr1_DB, '-k');
xlabel('$d$ (m)'); ylabel('$\min(R_1)$ (unitless)');
legend([pecham, ppred], 'AQUA', 'EBM', 'location', 'southeast');
set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
set(gca, 'xminortick','on', 'yminortick', 'on', 'xlim', xlim, 'ylim', [ymin ymax]);
print(sprintf('./figures_post/test/amp_r1_echam/amp_r1_B(d)_D(d)_echam_%s', par.choice), '-dpng', '-r300');
close;
