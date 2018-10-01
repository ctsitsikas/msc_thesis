%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: plottingallowances_waves.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%close all
addpath('/home/christos/csirolib/');
% Default options
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')

% IPCC_PW = load('../compute_allowances/allowances_ipcc_GTSR.mat');
IPCC_PW = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_waves_150_clean.mat');
%IPCC_PW = load('../compute_allowances/allowances_ipcc_9nov.mat');
IPCC = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_150_clean.mat');

%IPCC.allowances = horzcat(IPCC.muislonlat_sort,IPCC.allow_mix_muis100(:,1));
%IPCC.allowances = IPCC.allowances(~any(isnan(IPCC.allowances),2),:);
%IPCC.muislonlat_sort = IPCC.allowances(:,1:2);
%IPCC.allow_mix_PW = IPCC.allowances(:,3);
%allowgtsr = horzcat(IPCC.muislonlat_sort(:,1),IPCC.muislonlat_sort(:,2),IPCC.allow_normal_muis100(:,1));
%allowwaves = horzcat(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),IPCC_PW.allow_normal_gtsr_waves);
over = find(IPCC_PW.allow_normal_gtsr_waves>10);
IPCC_PW.gtsrwaves_sort(over,:)=[];
IPCC_PW.allow_normal_gtsr_waves(over)=[];
IPCC_PW.term1(over)=[];
IPCC_PW.term2(over)=[];
%clear IPCC.muislonlat_sort IPCC_PW
%IPCC.muislonlat_sort(:,1) = allowgtsr(:,1);
%IPCC.muislonlat_sort(:,2) = allowgtsr(:,2);
%IPCC.allow_normal_muis100 = allowgtsr(:,3);
%IPCC_PW.gtsrwaves_sort(:,1) = allowwaves(:,1);
%IPCC_PW.gtsrwaves_sort(:,2) = allowwaves(:,2);
%IPCC_PW.allow_normal_gtsr_waves = allowwaves(:,3);
% %% clean from weird points
% load('/home/christos/Dropbox/thesis/weird_points.mat');
% check = weirdpoints(:,1:2);
% [common ida idb] = intersect(allowgtsr(:,1:2), check, 'rows', 'stable');
% [common2 id2a id2b] = intersect(allowwaves(:,1:2), check, 'rows', 'stable');
% weirdgtsr = allowgtsr(ida,:);
% weirdwaves = allowwaves(id2a,:);
% 
% allowgtsr(ida,:) = [];
% allowwaves(id2a,:) = [];
% IPCC.muislonlat_sort(:,1) = allowgtsr(:,1);
% IPCC.muislonlat_sort(:,2) = allowgtsr(:,2);
% IPCC.allow_normal_muis100 = allowgtsr(:,3);
% IPCC_PW.gtsrwaves_sort(:,1) = allowwaves(:,1);
% IPCC_PW.gtsrwaves_sort(:,2) = allowwaves(:,2);
% IPCC_PW.allow_normal_gtsr_waves = allowwaves(:,3);
%%
%Plotting allowances
figure
subplot(2,1,1)
c1 = IPCC.allow_normal_muis100(:,1);
%c = IPCC.allow_mix_PW(:,1);
scatter(IPCC.muislonlat_sort(:,1),IPCC.muislonlat_sort(:,2),30,c1,'filled');
title('GTSR Allowances (m)','Fontsize',18)
ylabel('Latitude','Fontsize',16)
;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
col = colorbar
col.FontSize = 12
caxis([-0.1 2.5])
subplot(2,1,2)
c2 = IPCC_PW.allow_normal_gtsr_waves;
%c = IPCC_PW.allow_normal_PW(:,1);
%scatter(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),30,c,'filled');
scatter(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),30,c2,'filled');
title('GTSR including waves height maxima change Allowances (m)','Fontsize',18)

;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
col = colorbar
col.FontSize = 12
caxis([-0.1 2.5])
xlabel('Longitude','Fontsize',16)
ylabel('Latitude','Fontsize',16)
%%
%Plotting allowances
figure
subplot(2,1,1)
c1 = weirdgtsr(:,3);
%c = IPCC.allow_mix_PW(:,1);
scatter(weirdgtsr(:,1),weirdgtsr(:,2),30,c1,'filled');
title('GTSR Allowances (m) - repeated points')
ylabel('Latitude')
;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
colorbar
caxis([-0.1 2.5])
subplot(2,1,2)
c2 = weirdwaves(:,3);
%c = IPCC_PW.allow_normal_PW(:,1);
%scatter(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),30,c,'filled');
scatter(weirdwaves(:,1),weirdwaves(:,2),30,c2,'filled');
title('GTSR including waves Allowances (m) - repeated points')

;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
colorbar
caxis([-0.1 2.5])
xlabel('Longitude')
ylabel('Latitude')



%% differences
for i=1:length(IPCC.muislonlat_sort(:,1))
    if IPCC.muislonlat_sort(i,1) > 180
                IPCC.muislonlat_sort(i,1) = IPCC.muislonlat_sort(i,1) - 360;
    end
end
for i=1:length(IPCC_PW.gtsrwaves_sort(:,1))
    if IPCC_PW.gtsrwaves_sort(i,1) > 180
                IPCC_PW.gtsrwaves_sort(i,1) = IPCC_PW.gtsrwaves_sort(i,1) - 360;
    end
end
clear C ia ib

[C,ia,ib] = intersect(IPCC_PW.gtsrwaves_sort(:,1:2),IPCC.muislonlat_sort(:,1:2),'rows');
figure
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
% Apply the colormap.
colormap(colorMap);
scatter(IPCC_PW.gtsrwaves_sort(ia,1),IPCC_PW.gtsrwaves_sort(ia,2),30,c2(ia) - c1(ib),'filled'); hold on;

col = colorbar
caxis([-0.5 0.5])
col.FontSize=12
;hold on;
draw_land(0, [0.8 0.8 0.8], 1, 'k');
xlim([  -15 35]);
    ylim([34 72]);
xlabel('Longitude','Fontsize',16)
ylabel('Latitude','Fontsize',16)
title('Difference between GTSR allowances with and without wave change' ,'Fontsize',18)
%%
%Plotting slc
gamma = 0.5772;

mean_waves = IPCC_PW.wavesatgtsr_sort(:,3) + IPCC_PW.wavesatgtsr_sort(:,4).*gamma;
std_waves = (pi./sqrt(6)).*IPCC_PW.wavesatgtsr_unique(:,4);


figure
subplot(2,1,1)
c1 = (IPCC.std_SLR.^2 + std_waves.^2)./(2.*IPCC.muisscale_unique);
%c = IPCC.allow_mix_PW(:,1);
scatter(IPCC.muislonlat_unique(:,1),IPCC.muislonlat_unique(:,2),30,c1,'filled');
title('Allowances variability term')
ylabel('Latitude')
;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
colorbar
caxis([-0.5 1])
subplot(2,1,2)
c2 = IPCC_PW.mean_SLR + mean_waves;
%c = IPCC_PW.allow_normal_PW(:,1);
%scatter(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),30,c,'filled');
scatter(IPCC_PW.gtsrwaves_unique(:,1),IPCC_PW.gtsrwaves_unique(:,2),30,c2,'filled');
title('Allowances mean term')

;hold on;
coast('k')
%xlim([  -15 35]);
%    ylim([34 72]);
xlim([0 360])
ylim([-90 90])
colorbar
caxis([0 1])
xlabel('Longitude')
ylabel('Latitude')