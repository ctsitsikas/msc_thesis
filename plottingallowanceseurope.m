%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: plottingallowanceseurope.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting europe
clear all
close all
% Default options
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')
addpath('/home/christos/csirolib/');

%GESLA2 allowances
%load('/home/christos/Desktop/for_Roderik_from_Aimee/figure_plot_files/nanidsort.mat')
%IPCC_PW = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_9nov.mat');
%IPCC_PW.G2_PW_lonlat_sort(nanidxsort,:)=[];
%IPCC_PW.allow_normal_PW(nanidxsort,:)=[];

%GTSR+wave allowances
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_waves_150_clean.mat');
IPCC_PW.G2_PW_lonlat_sort = gtsrwaves_sort;
IPCC_PW.allow_normal_PW = allow_normal_gtsr_waves;

%GTSR allowances
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_150_clean.mat');
IPCC_PW.G2_PW_lonlat_sort = gtsrwaves_sort;
IPCC_PW.allow_normal_PW = allow_normal_gtsr_waves;

%%
for i = 1:length(IPCC_PW.G2_PW_lonlat_sort(:,1))
    if IPCC_PW.G2_PW_lonlat_sort(i,1) > 330
        IPCC_PW.G2_PW_lonlat_sort(i,1) = IPCC_PW.G2_PW_lonlat_sort(i,1)-360;
    end
end
for i = 1:length(muislonlat_sort(:,1))
    if muislonlat_sort(i,1) > 330
        muislonlat_sort(i,1) = muislonlat_sort(i,1)-360;
    end
end

%%
%Plotting allowances
figure
subplot(2,1,2)
c = allow_normal_muis150(:,1);
%c = allow_mix_PW(:,1);
scatter(muislonlat_sort(:,1),muislonlat_sort(:,2),30,c,'filled');
title('GTSR Allowances, Europe (m)','Fontsize',18)
ylabel('Latitude','Fontsize',16)
;hold on;
draw_land(0, [0.8 0.8 0.8], 1, 'k');
xlim([  -15 35]);
    ylim([34 72]);
%xlim([0 360])
%ylim([-90 90])
col = colorbar
col.FontSize=12
caxis([0 1.2])
subplot(2,1,1)
c = IPCC_PW.allow_normal_PW(:,1);
scatter(IPCC_PW.G2_PW_lonlat_sort(:,1),IPCC_PW.G2_PW_lonlat_sort(:,2),30,c,'filled');
%title('GTSR with waves height maxima change Allowances with 150km limit, Europe (m)','Fontsize',18)
title('GESLA-2 Allowances, Europe (m)','Fontsize',18)
;hold on;
draw_land(0, [0.8 0.8 0.8], 1, 'k','crd');
xlim([  -15 35]);
    ylim([34 72]);
%xlim([0 360])
%ylim([-90 90])
col = colorbar
col.FontSize=12
caxis([0 1.2])
xlabel('Longitude','Fontsize',16)
ylabel('Latitude','Fontsize',16)