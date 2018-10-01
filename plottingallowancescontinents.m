%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: plottingallowancescontinents.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
addpath('/home/christos/csirolib/');
% Default options
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')
%IPCC = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_gtsr_waves_clean.mat');
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_150_clean.mat')
% name change when choosing GTSR+wave allowances and filtering of very high allowances
% muislonlat_sort = gtsrwaves_sort;	
%  allow_normal_muis100 = allow_normal_gtsr_waves; 
% ids = find(allow_normal_muis100>8);
%  allow_normal_muis100(ids) = [];
% muislonlat_sort(ids,:) = [];
 
for i = 1:length(muislonlat_sort(:,1))
    if muislonlat_sort(i,1) > 180
        muislonlat_sort(i,1) = muislonlat_sort(i,1)-360;
    end
end
%% plotting continents on map

figure

c = allow_normal_muis100(:,1);
%c = allow_normal_PW(:,1);
scatter(muislonlat_sort(:,1),muislonlat_sort(:,2),30,c,'filled');
title('GTSR Allowances (m)','Fontsize',18)
ylabel('Latitude','Fontsize',16)
xlabel('Longitude','Fontsize',16)
;hold on;

draw_land(0, [0.8 0.8 0.8], 'k');
grid on
;hold on;
a = rectangle('Position',[-15 34 35+15 72-34], 'EdgeColor','blue');
a.LineWidth = 2;
%b = rectangle('Position',[-20 -35 55+20 34+35],'EdgeColor','k');
%b.LineWidth = 2;
plot([-20 -20 35 35 55 55 -20], [-35 34 34 5 5 -35 -35], 'color', 'k', 'Linewidth',2)

c = rectangle('Position',[35 5 180-35 90-5],'EdgeColor','y');
c.LineWidth = 2;
d = rectangle('Position',[90 -48 90 48],'EdgeColor','green'); %oceania
d.LineWidth = 2;
e = rectangle('Position',[-168 13 -20+168 90-13],'EdgeColor','r');
e.LineWidth = 2;
f = rectangle('Position',[-92 -57 -32+92 13+57],'EdgeColor','r');
f.LineWidth = 2;
g = rectangle('Position',[-180 -90 360 30],'EdgeColor','c'); %antarctica
g.LineWidth = 2;
xlim([-180 180]);
ylim([-90 90]);
colorbar
caxis([0 2.5])
%%
% figure
% 
% c = europeallow(:,1);
% %c = allow_normal_PW(:,1);
% scatter(europelonlat(:,1),europelonlat(:,2),30,c,'filled');
% title('GTSR Allowances, Europe')
% ylabel('Latitude')
% ;hold on;
% 
% draw_land(0, [0.8 0.8 0.8], 'k');
% colorbar
%% Europe
figure
ideurope = muislonlat_sort(:,2)>34 & muislonlat_sort(:,2)<72 & muislonlat_sort(:,1)>-15 & muislonlat_sort(:,1)<35;
europelonlat = muislonlat_sort(ideurope,1:2);
europeallow = allow_normal_muis100(ideurope,1);
scatter(europelonlat(:,2),europeallow,'filled')

ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
P = polyfit(europelonlat(:,2),europeallow,1);
yfit = P(1)*europelonlat(:,2)+P(2);
hold on;
plot(europelonlat(:,2),yfit,'r-.','Linewidth',2);
str = sprintf('GTSR and wave height maxima change, Europe, slope=%.3f',round(P(1),3))
title(str,'Fontsize',18)


%% Africa
clear idafrica africalonlat africaallow
figure
idafrica = muislonlat_sort(:,2)<34 & muislonlat_sort(:,2)>-35 & muislonlat_sort(:,1)>-20 & muislonlat_sort(:,1)<55;
africalonlat = muislonlat_sort(idafrica,1:2);
africaallow = allow_normal_muis100(idafrica,1);
ind = find(africalonlat(:,1)>35 &  africalonlat(:,2)>5);
africalonlat(ind,:) = []; 
africaallow(ind) = []; 

% partial trends
africa = horzcat(africalonlat, africaallow);
africa = sortrows(africa,2);
index = find(africa(:,2) > 0,1);
africalonlat = africa(:,1:2);
africaallow = africa(:,3);


scatter(africalonlat(:,2),africaallow,'filled')
%title('GTSR, Africa')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
P = polyfit(africalonlat(1:index,2),africaallow(1:index),1);
yfit = P(1).*africalonlat(1:index,2)+P(2);
hold on;
plot(africalonlat(1:index,2),yfit,'r-.','Linewidth',2);

P2 = polyfit(africalonlat(index:length(africaallow),2),africaallow(index:length(africaallow)),1);
yfit2 = P2(1).*africalonlat(index:length(africaallow),2)+P2(2);
hold on;
plot(africalonlat(index:length(africaallow),2),yfit2,'r-.','Linewidth',2);
str = sprintf('GTSR, Africa, slope SH = %.3f, slope NH = %.3f',P(1), P2(1))
title(str,'Fontsize',18)
%% Asia
figure
idasia = muislonlat_sort(:,2)<90 & muislonlat_sort(:,2)>5 & muislonlat_sort(:,1)>35 & muislonlat_sort(:,1)<180;
asialonlat = muislonlat_sort(idasia,1:2);
asiaallow = allow_normal_muis100(idasia,1);
scatter(asialonlat(:,2),asiaallow,'filled')
%title('GTSR, Asia')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
xlim([0 90])
P = polyfit(asialonlat(:,2),asiaallow,1);
yfit = P(1)*asialonlat(:,2)+P(2);
hold on;
plot(asialonlat(:,2),yfit,'r-.','Linewidth',2);
str = sprintf('GTSR and wave height maxima change, Asia, slope=%.3f',P(1))
title(str,'Fontsize',18)
%% oceania
figure
idoceania = muislonlat_sort(:,2)<0 & muislonlat_sort(:,2)>-48 & muislonlat_sort(:,1)>90 & muislonlat_sort(:,1)<180;
oceanialonlat = muislonlat_sort(idoceania,1:2);
oceaniaallow = allow_normal_muis100(idoceania,1);
scatter(oceanialonlat(:,2),oceaniaallow,'filled')
%title('GTSR, Oceania')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
P = polyfit(oceanialonlat(:,2),oceaniaallow,1);
yfit = P(1)*oceanialonlat(:,2)+P(2);
hold on;
plot(oceanialonlat(:,2),yfit,'r-.','Linewidth',2);
str = sprintf('GTSR and wave height maxima change, Oceania, slope=%.3f',P(1))
title(str,'Fontsize',18)
%% northam
figure
idnortham = muislonlat_sort(:,2)<90 & muislonlat_sort(:,2)>13 & muislonlat_sort(:,1)>-168 & muislonlat_sort(:,1)<-20;
northamlonlat = muislonlat_sort(idnortham,1:2);
northamallow = allow_normal_muis100(idnortham,1);
scatter(northamlonlat(:,2),northamallow,'filled')
%title('GTSR, North America')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
P = polyfit(northamlonlat(:,2),northamallow,1);
yfit = P(1)*northamlonlat(:,2)+P(2);
hold on;
plot(northamlonlat(:,2),yfit,'r-.','Linewidth',2);
str = sprintf('GTSR and wave height maxima change, North America, slope=%.3f',P(1))
title(str,'Fontsize',18)
%% southam
figure
idsoutham = muislonlat_sort(:,2)<13 & muislonlat_sort(:,2)>-56 & muislonlat_sort(:,1)>-92 & muislonlat_sort(:,1)<-32;
southamlonlat = muislonlat_sort(idsoutham,1:2);
southamallow = allow_normal_muis100(idsoutham,1);
southam = horzcat(southamlonlat, southamallow);
southam = sortrows(southam,2);
index = find(southam(:,2) > 0,1);
southamlonlat = southam(:,1:2);
southamallow = southam(:,3);
scatter(southamlonlat(:,2),southamallow,'filled')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
%title('GTSR, South America')
P = polyfit(southamlonlat(1:index,2),southamallow(1:index),1);
yfit = P(1).*southamlonlat(1:index,2)+P(2);
P2 = polyfit(southamlonlat(index:length(southamallow),2),southamallow(index:length(southamallow)),1);
yfit2 = P2(1).*southamlonlat(index:length(southamallow),2)+P2(2);
hold on;
plot(southamlonlat(1:index,2),yfit,'r-.','Linewidth',2);
%hold on;
%plot(southamlonlat(index:length(southamallow),2),yfit2,'r-.','Linewidth',2);
str = sprintf('GTSR and wave height maxima change, South America, slope=%.3f',P(1))
title(str,'Fontsize',18)
%% antarct
figure
idantarct = muislonlat_sort(:,2)<-60 & muislonlat_sort(:,2)>-90 & muislonlat_sort(:,1)>-180 & muislonlat_sort(:,1)<180;
antarctlonlat = muislonlat_sort(idantarct,1:2);
antarctallow = allow_normal_muis100(idantarct,1);
scatter(antarctlonlat(:,2),antarctallow,'filled')
%title('GTSR, Antarctica')
ylabel('Allowance (m)','Fontsize',16)
xlabel('Latitude','Fontsize',16)
P = polyfit(antarctlonlat(:,2),antarctallow,1);
yfit = P(1)*antarctlonlat(:,2)+P(2);
hold on;
plot(antarctlonlat(:,2),yfit,'r-.','Linewidth',2);
str = sprintf('GTSR, Antarctica, slope=%.3f',P(1))
title(str,'Fontsize',18)
