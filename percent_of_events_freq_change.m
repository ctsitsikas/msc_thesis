%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: percent_of_events_freq_change.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%choose to upload GESLA2 (allowances_ipcc_9nov) where first needs to remove
%non gumbel points(nanidsort), GTSR, or GTSR+waves

%GESLA-2
%load('/home/christos/Desktop/for_Roderik_from_Aimee/figure_plot_files/nanidsort.mat')
%IPCC_PW = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_9nov.mat');
%IPCC_PW.G2_PW_lonlat_sort(nanidxsort,:)=[];
%IPCC_PW.G2_PW_scale_sort(nanidxsort,:)=[];
%IPCC_PW.allow_normal_PW(nanidxsort,:)=[];

%GTSR+waves
%IPCC_PW = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_waves_150_clean.mat');

%GTSR
IPCC = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_150_clean.mat');
IPCC_PW.gtsrwaves_sort = IPCC.muislonlat_sort;
IPCC_PW.allow_normal_gtsr_waves = IPCC.allow_normal_muis150;
IPCC_PW.muisscale_sort = IPCC.muisscale_sort;

%% remove GTSR+waves allowances higher than x meters
allowwaves = horzcat(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),IPCC_PW.allow_normal_gtsr_waves,IPCC_PW.muisscale_sort(:,1));
limit = 8;
over = find(allowwaves(:,3)>limit);
allowwaves(over,:)=[];
IPCC_PW.gtsrwaves_sort(:,1) = allowwaves(:,1);
IPCC_PW.gtsrwaves_sort(:,2) = allowwaves(:,2);
IPCC_PW.allow_normal_gtsr_waves = allowwaves(:,3);
IPCC_PW.muisscale_sort = allowwaves(:,4);

%% if need to cut out pacific 
%load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/gtsr_waves_allowances_no_weirds.mat')
% %A
% ind = IPCC_PW.gtsrwaves_sort(:,2)>-24 & IPCC_PW.gtsrwaves_sort(:,2)<25 & IPCC_PW.gtsrwaves_sort(:,1)<226 & IPCC_PW.gtsrwaves_sort(:,1)>152;
% IPCC_PW.gtsrwaves_sort(ind,:) = []; 
% IPCC_PW.allow_normal_gtsr_waves(ind,:) = [];
% IPCC_PW.muisscale_sort(ind,:) = [];
% 
% %B
% ind = IPCC_PW.gtsrwaves_sort(:,2)>-6.5 & IPCC_PW.gtsrwaves_sort(:,2)<25 & IPCC_PW.gtsrwaves_sort(:,1)<152 & IPCC_PW.gtsrwaves_sort(:,1)>150;
% IPCC_PW.gtsrwaves_sort(ind,:) = []; 
% IPCC_PW.allow_normal_gtsr_waves(ind,:) = [];
% IPCC_PW.muisscale_sort(ind,:) = [];
% %C
% ind = IPCC_PW.gtsrwaves_sort(:,2)>-6.5 & IPCC_PW.gtsrwaves_sort(:,2)<10 & IPCC_PW.gtsrwaves_sort(:,1)<150 & IPCC_PW.gtsrwaves_sort(:,1)>130;
% IPCC_PW.gtsrwaves_sort(ind,:) = []; 
% IPCC_PW.allow_normal_gtsr_waves(ind,:) = [];
% IPCC_PW.muisscale_sort(ind,:) = [];
% %D
% ind = IPCC_PW.gtsrwaves_sort(:,2)>-6.5 & IPCC_PW.gtsrwaves_sort(:,2)<5 & IPCC_PW.gtsrwaves_sort(:,1)<130 & IPCC_PW.gtsrwaves_sort(:,1)>116;
% IPCC_PW.gtsrwaves_sort(ind,:) = []; 
% IPCC_PW.allow_normal_gtsr_waves(ind,:) = [];
% IPCC_PW.muisscale_sort(ind,:) = [];

%% GTSR
scale = IPCC.muisscale_sort(:,1);
allowances = IPCC.allow_normal_muis150;

N = exp(allowances./scale);
scatter(IPCC.muislonlat_sort(:,1),IPCC.muislonlat_sort(:,2),15,N) %plot change of frequency for all points
colorbar

% choose fold number of frequency (here 10^3) to find percent of points
% affected
percent = length(N(N>10^3))/length(N);

%% GESLA2
scalegesla = IPCC_PW.G2_PW_scale_sort(:,1);
allowancesgesla = IPCC_PW.allow_normal_PW;

Ng = exp(allowancesgesla./scalegesla);

percentgesla = length(Ng(Ng>10^3))/length(Ng);

%% GTSR+waves
clear scale allowances percent N
scale = IPCC_PW.muisscale_sort(:,1);
allowances = IPCC_PW.allow_normal_gtsr_waves;

N = exp(allowances./scale);
scatter(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),20,N)
colorbar

percent = length(N(N>10^3))/length(N);