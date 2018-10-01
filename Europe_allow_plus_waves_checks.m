%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Europe_allow_plus_waves_checks.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
% Default options for fontsizes
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')
%addpath('/home/christos/csirolib/');
% load allowances
load('compute_allowances/allowances_ipcc_GTSR_150_clean.mat');
load('compute_allowances/allowances_ipcc_GTSR_waves_150_clean.mat');
% find common points between sets
[C,ia,ib] = intersect(gtsrwaves_sort(:,1:2),muislonlat_sort(:,1:2),'rows','stable');
allowgtsr = allow_normal_muis150;
allowwave = allow_normal_gtsr_waves;
% find mean and std of SLR at all points and then keep the common ones 
for i = 1:length(idx_out)
    mSLR(i) = mean_SLR(idx_out(i));
    sSLR(i) = std_SLR(idx_out(i));
end
mSLR = mSLR(ia)';
sSLR = sSLR(ia)';

Da = allowwave(ia)-allowgtsr(ib);
mean_waves = mean_waves(ia)'; % mean and std of wave distribution at commont points
std_waves = std_waves(ia)';

id0 = find(Da~=0); %find indices of points with change due to waves and keep only these
Da = Da(id0);
mSLR = mSLR(id0);
sSLR = sSLR(id0);
mean_waves = mean_waves(ia);
std_waves = std_waves(ia);
mean_waves = mean_waves(id0);
std_waves = std_waves(id0);
gtsrwaves_sort = gtsrwaves_sort(ia,1:2);
muisscale_sort = muisscale_sort(ia,1);
muislonlat_sort = muislonlat_sort(ib,1:2);
gtsrwaves_sort = gtsrwaves_sort(id0,1:2);
muislonlat_sort = muislonlat_sort(id0,1:2);
muisscale_sort = muisscale_sort(id0);
% part1 = abs(mean_waves(id0))./abs(Da);
% part2 = std_waves(id0)./abs(Da);

data = horzcat(gtsrwaves_sort,Da); 

data = sortrows(data,2);


for i = 1:length(gtsrwaves_sort(:,1))
    if data(i,1) > 330
        data(i,1) = data(i,1)-360;
    end
end


%%

figure %choose to plot Da or contribution by mean or std

%scatter(data(:,1),data(:,2),30,(abs(mean_waves)./(abs(mean_waves)+(std_waves./muisscale_sort)))*100,'filled');
scatter(data(:,1),data(:,2),30,data(:,3),'filled')
%scatter(data(:,1),data(:,2),30,((std_waves./muisscale_sort)./(abs(mean_waves)+(std_waves./muisscale_sort)))*100,'filled');
%title('Contribution of mean wave maxima change to allowances, Europe (%)','Fontsize',18)

;hold on;
%coast('k')
draw_land(0, [0.8 0.8 0.8], 1, 'k');
xlim([  -35 35]); %limit to Europe
ylim([34 72]);
% xlim([0 360])
% ylim([-90 90])
col = colorbar
col.FontSize=12
caxis([0 1])
xlabel('Longitude')
ylabel('Latitude')

%%
dataeur = data(data(:,1)>-35 & data(:,1)<35,:); %limit data to European longitudes
scatter(dataeur(:,2),dataeur(:,3),'filled')
xlim([34 72]);
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]

line(get(axh,'XLim'), [0 0], 'Color', color); %plot line at y=0
ylim([-1 3.5])
xlabel('Latitude')
ylabel('Allowance (m)')
title('Contribution of wave change to allowances against latitude, Europe (m)')
