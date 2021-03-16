clear all
close all
load physics
rearth=ps.a;
LV=ps.LV;
lat_int=[90:-.05:-90]';
Z_int=[0:500:20e3]';

file_1x_atm=[char(strcat('/project2/tas1/echam6/piControl/ATM_piControl_r1i1p1-P_1900_1949.nc'))];
file_1x_atm_rh=[char(strcat('/project2/tas1/echam6/piControl/ATM_piControl_r1i1p1-P_1900_1949.rhum.nc'))];
file_1x=[char(strcat('/project2/tas1/echam6/piControl/BOT_piControl_r1i1p1-P_1900_1949.nc'))];
file_sb=[char(strcat('/project2/tas1/arejaygraham/6HourSims/rjg_20170908/BOT_rjg_20170908_0078_0083_seas.nc'))];

lat=ncread(file_1x,'lat');
lon=ncread(file_1x,'lon');
lev=ncread(file_1x_atm,'lev');
w=squeeze(ncread(file_1x_atm,'omega'));
T=squeeze(ncread(file_1x_atm,'t'));
Z=squeeze(ncread(file_1x_atm,'geopoth'));
RH=squeeze(ncread(file_1x_atm_rh,'rhum'));
aps=squeeze(ncread(file_1x,'aps'));

toa_1x=squeeze(nanmean(ncread(file_1x,'trad0')+ncread(file_1x,'srad0'),1));
ei_1x=squeeze(nanmean(ncread(file_1x,'trad0')+ncread(file_1x,'srad0')-(ncread(file_1x,'srads')+ncread(file_1x,'trads')+ncread(file_1x,'ahfs')+ncread(file_1x,'ahfl')),1));
ei_dse_1x=squeeze(nanmean(LV*ncread(file_1x,'precip')+ncread(file_1x,'trad0')+ncread(file_1x,'srad0')-(ncread(file_1x,'srads')+ncread(file_1x,'trads')+ncread(file_1x,'ahfs')),1));
swabs_1x=squeeze(nanmean(ncread(file_1x,'srad0')-(ncread(file_1x,'srads')),1));
ghe_1x=squeeze(nanmean(ncread(file_1x,'trad0')-(ncread(file_1x,'trads')),1));
shf_1x=-squeeze(nanmean(ncread(file_1x,'trads')+ncread(file_1x,'ahfs')+ncread(file_1x,'ahfl'),1));
hf_1x=-squeeze(nanmean(ncread(file_1x,'ahfs')+ncread(file_1x,'ahfl'),1));
ahfl_1x=-squeeze(nanmean(ncread(file_1x,'ahfl'),1));
ahfs_1x=-squeeze(nanmean(ncread(file_1x,'ahfs'),1));
albedo_1x=squeeze(nanmean(ncread(file_1x,'albedo'),1));
trads_1x=-squeeze(nanmean(ncread(file_1x,'trads'),1));
pme_1x=squeeze(nanmean(LV*ncread(file_1x,'precip')+(ncread(file_1x,'ahfl')),1));
Lp_1x=LV*squeeze(nanmean(ncread(file_1x,'precip'),1));
Lpc_1x=LV*squeeze(nanmean(ncread(file_1x,'aprc'),1));
Lpl_1x=LV*squeeze(nanmean(ncread(file_1x,'aprl'),1));

for m=1:12

Ra_1x(:,m)=ghe_1x(:,m)+swabs_1x(:,m);
R1_1x(:,m)=(ei_1x(:,m)./(Ra_1x(:,m)));

R1_4050N(m) = sum(cos(lat(22:27)*pi/180).*R1_1x(22:27,m),1)./sum(cos(lat(22:27)*pi/180));
Ra_4050N(m) = sum(cos(lat(22:27)*pi/180).*Ra_1x(22:27,m),1)./sum(cos(lat(22:27)*pi/180));
divF_4050N(m) = sum(cos(lat(22:27)*pi/180).*ei_1x(22:27,m),1)./sum(cos(lat(22:27)*pi/180));
R1_4050S(m) = sum(cos(lat(70:75)*pi/180).*R1_1x(70:75,m),1)./sum(cos(lat(70:75)*pi/180));
Ra_4050S(m) = sum(cos(lat(70:75)*pi/180).*Ra_1x(70:75,m),1)./sum(cos(lat(70:75)*pi/180));
divF_4050S(m) = sum(cos(lat(70:75)*pi/180).*ei_1x(70:75,m),1)./sum(cos(lat(70:75)*pi/180));
end

DeltaR1_4050N=R1_4050N-mean(R1_4050N);
DeltaRa_4050N=-(Ra_4050N-mean(Ra_4050N)).*(mean(divF_4050N)/(mean(Ra_4050N).^2));
DeltadivF_4050N=(divF_4050N-mean(divF_4050N))/(mean(Ra_4050N));
DeltaR1_4050S=R1_4050S-mean(R1_4050S);
DeltaRa_4050S=-(Ra_4050S-mean(Ra_4050S)).*(mean(divF_4050S)/(mean(Ra_4050S).^2));
DeltadivF_4050S=(divF_4050S-mean(divF_4050S))/(mean(Ra_4050S));

figure(1)
plot(DeltaR1_4050N,'k','LineWidth',3)
hold on
plot(DeltadivF_4050N,'b','LineWidth',3)
plot(DeltaRa_4050N,'r','LineWidth',3)
plot(DeltaR1_4050N-(DeltaRa_4050N+DeltadivF_4050N),'k--','LineWidth',3)
legend('\DeltaR_1','\DeltadivF','\DeltaRa','Location','SouthEast')
plot(0.0*DeltaR1_4050S,'k')
set(gca,'xlim',[1 12]);
set(gca,'XTick',[1:1:20]);
set(gca,'ylim',[-.5 .5]);
set(gca,'YTick',[-0.5:.25:0.5]);
set(gca,'fontsize',30);
title('40-50N')
set(gca,'XTickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'],'fontsize',30)


figure(2)
plot(circshift(DeltaR1_4050S,[0 6]),'k','LineWidth',3)
hold on
plot(circshift(DeltadivF_4050S,[0 6]),'b','LineWidth',3)
plot(circshift(DeltaRa_4050S,[0 6]),'r','LineWidth',3)
plot(circshift(DeltaR1_4050S-(DeltaRa_4050S+DeltadivF_4050S),[0 6]),'k--','LineWidth',3)
legend('\DeltaR_1','\DeltadivF','\DeltaRa')
plot(0.0*DeltaR1_4050S,'k')
set(gca,'xlim',[1 12]);
set(gca,'XTick',[1:1:20]);
set(gca,'ylim',[-.5 .5]);
set(gca,'YTick',[-0.5:.25:0.5]);
set(gca,'fontsize',30);
title('40-50S')
set(gca,'XTickLabel',['J';'A';'S';'O';'N';'D';'J';'F';'M';'A';'M';'J'],'fontsize',30)
