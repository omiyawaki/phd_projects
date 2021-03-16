clear all
close all
% load physics
% rearth=ps.a;
% LV=ps.LV;
lat_int=[90:-.25:-90]';

lat_ceres=ncread('/project2/tas1/CERES/CERES_EBAF-TOA_Ed4.0_Subset_200003-201803.nc','lat');
Ra_toa_ceres=squeeze(mean(ncread('/project2/tas1/CERES/CERES_EBAF-TOA_Ed4.0_Subset_200003-201803.nc','toa_net_all_mon'),1));
Ra_sfc_ceres=squeeze(mean(ncread('/project2/tas1/CERES/CERES_EBAF-Surface_Ed4.0_Subset_200003-201803.nc','sfc_net_tot_all_mon'),1));

size(Ra_sfc_ceres)

return
Ra_ceres=Ra_toa_ceres-Ra_sfc_ceres;

for i=2:18
[10+1*(i-1)+11*(i-2),10+1*(i-1)+11+11*(i-2)];
Ra_ceres_month_all(:,:,i-1)=squeeze(Ra_ceres(:,10+1*(i-1)+11*(i-2):10+1*(i-1)+11+11*(i-2)));
end
Ra_ceres_month=squeeze(mean(Ra_ceres_month_all,3));

load radiation_dynamics_climatology.mat
ht=squeeze(mean(TETEN,3));
Fa=squeeze(mean(TEDIV,3));

for m=1:12
ei_int(:,m)=interp1(lat,squeeze(Fa(m,:)+ht(m,:)),lat_int,'spline');
Ra_int(:,m)=interp1(lat_ceres,squeeze(Ra_ceres_month(:,m)),lat_int,'spline');
end

rce=0.0*ei_int;
rae=0.0*ei_int;
for m=1:12
R1(:,m)=(ei_int(:,m)./(Ra_int(:,m)));
for i=1:numel(lat_int)

if(R1(i,m)<=0.3 & R1(i,m)>=-0.3)
rce(i,m) = 15.0;
end

if(R1(i,m)>=0.6)
rae(i,m) = 10.0;
end

end
end

rce(:,13)=rce(:,1);
rae(:,13)=rae(:,1);

rce_rae = rce+rae;
month=[1:1:13];

uscalefac=1.0;
cmax=10/uscalefac;
cintu=1/uscalefac;
cmp=colCog(2*cmax/cintu);
colormap(cmp)

set(gcf,'Renderer','painters')
figure(1)
contourf(month,lat_int,rce,[15 15])
hold on
contourf(month,lat_int,rae,[10 10])
shading flat
caxis([5 cmax+10]);
set(gca,'XLim',[1 13],'YLim',[-90 90], ...
    'XTick',[1:1:13], ...
    'Ytick',[-90:30:90]);
set(gca,'fontsize',40);
ylabel('latitude');
title('RCE -0.3<= R1 <= 0.3, RAE >= 0.6')
set(gca,'XTickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D';'J'],'fontsize',40)

set(gcf,'paperunits','inches');
set(gcf,'papersize',[17 8]);
set(gcf,'paperposition',[0 0 17 8]);
print('-depsc','-r300','rce_rae_erainterim_v2.eps');


