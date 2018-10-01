%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: allowances_with_waves.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%cd Desktop/for_Roderik_from_Aimee/compute_allowances/
%%This script computes allowances for a NORMAL/GAUSSIAN distribution
fprintf('-----> Start \n')
%%
load('/home/christos/Dropbox/thesis/gtsr_waves_near_clean.mat');
%  muislat = gtsr_clean(:,2);
%  muislon = gtsr_clean(:,1);
%  loc = gtsr_clean(:,3);
%  muisscale = gtsr_clean(:,4);
gtsr = gtsratnearestwaves_clean;
% gtsr = horzcat(muislon,muislat,loc, muisscale);
gtsr(any(isnan(gtsr), 2), :) = [];

muislon = gtsr(:,1);
muislat = gtsr(:,2);
loc = gtsr(:,3);
muisscale = gtsr(:,4);

% Step 0 (only once): determine grid points closest to TG location%
%COMMENT THIS STEP WHEN IT IS DONE
%max distance kept

limit = 150000;
%This step finds the closest ocean grid point for each tide gauge
%Load original Tide gauge data (locations, gumbel scaling, names, number of years)
%load('Step_1/muis.mat'); %this is the GESLA2 tide gauge data from Phil Woodworth.

for i=1:length(muislon)
    if muislon(i) < 0
        muislon(i) = muislon(i) + 360;
    end
end
%muis = horzcat(gtsrwaves,muisscale,loc);
%muis = muis(~any(isnan(muis),2),:);

%muis = muis(idx_sw,:);


gtsrwaves = horzcat(muislon,muislat);
% muisscale = muis(:,3);
% loc = muis(:,4);
%gtsrwaves(:,2) = gtsrwaves(:,2) -90; %correct latitudes to be similar to GESLA
%load land-ocean mask from sea-level fields
temp = ncread('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/sealevel-change-fields-at-2100-slangen.nc','SEtotalB');
temp1 = ncread('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/sealevel-change-fields-at-2100-slangen.nc','totalB');
for i=1:360
    for j =1:180
        if temp(i,j)==-9999 || temp1(i,j)==-9999
            mask_in(i,j) = NaN;
        else
            mask_in(i,j) = temp(i,j);
        end
    end
end


mask = ~isnan(mask_in);
lon = ncread('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/sealevel-change-fields-at-2100-slangen.nc','lon');
lat = ncread('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/sealevel-change-fields-at-2100-slangen.nc','lat');
masklon = mask .* repmat(lon,[1 180]);
masklon(masklon == 0) = NaN;
masklat = mask .* repmat(lat',[360 1]);
masklat(masklat == 0) = NaN;
clear lon lat
a = 6371221;
d2r = 3.1415927/180;
% now matching the tide gauge lat-lon to map locations
for station_nr = 1:length(gtsrwaves(:,1));
    %londiff(:,:) = abs(repmat(gtsrwaves(station_nr,1),[360 180])-masklon(:,:)).*cosd(gtsrwaves(station_nr,2)).*(a*(3.14159/180)); %compute longitude distances
    %latdiff(:,:) = abs(repmat(gtsrwaves(station_nr,2),[360 180])-masklat(:,:)).*(a.*(3.14159/180)); %compute latitude distances
    totdiff = a.*acos(cos(d2r.*repmat(gtsrwaves(station_nr,2),[360 180])).*cos(d2r.*masklat(:,:)).*cos(d2r.*(repmat(gtsrwaves(station_nr,1),[360 180])-masklon(:,:)))+sin(d2r.*repmat(gtsrwaves(station_nr,2),[360 180])).*sin(d2r.*masklat(:,:)));
    %totdiff = sqrt(londiff.^2 + latdiff.^2); %compute total distance
    [num idx] = min(totdiff(:)); % find minimum distance
    [lonlat_map(station_nr,1) lonlat_map(station_nr,2)] = ind2sub(size(totdiff),idx); %write array location of minimum
    gtsrwaves(station_nr,3) = masklon(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lon
    gtsrwaves(station_nr,4) = masklat(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lat
    gtsrwaves(station_nr,5) = num;
    muisscale(station_nr,2) = num;
    loc(station_nr,2) = num;
end

lonlat_map = horzcat(lonlat_map,gtsrwaves(:,5));
lonlat_map(lonlat_map(:,3) > limit, :) = [];
lonlat_map = lonlat_map(:,1:2);
gtsrwaves(gtsrwaves(:,5) > limit, :) = [];
muisscale(muisscale(:,2) > limit, :) = [];
loc(loc(:,2) > limit, :) = [];
%manual corrections (Hudson Bay)
%lonlat_map(535,2) = 127;
[a b c] =intersect(gtsrwaves(:,1:2),wavesatgtsr(:,1:2),'rows','stable');
wavesatgtsr=wavesatgtsr(c,:);
       %% check locations in map to see if they match
        %surface(mask)
        %hold on
        %scatter(lonlat_map(:,2),lonlat_map(:,1),'og')
        %scatter(gtsrwaves(:,2)+90,gtsrwaves(:,1),'rx')
        %shading interp

        %Because there are sometime more tide gauges in the same location, we sort the lat-lons to reduce computing time later
        [lonlat_map_sort, index] = sortrows(lonlat_map,[2 1]);
        gtsrwaves_sort = gtsrwaves(index,:);
        wavesatgtsr_sort = wavesatgtsr(index,:);
        %G2_PW_nyr_sort = G2_PW_nyr(index,:);
        muisscale_sort = muisscale(index,:);
        loc_sort = loc(index,:);
        %G2_PW_filen_sort = G2_PW_filen(index,:);
        clear i idx index j latdiff londiff num station_nr totdiff
        
        %FIND UNIQUE GRIDPOINTS FOR RUNNING seawise and reducing analysis
        [lonlat_map_unique, idx_in, idx_out] = unique(lonlat_map_sort,'rows','stable');
        gtsrwaves_unique = gtsrwaves_sort(idx_in,:);
        wavesatgtsr_unique = wavesatgtsr_sort(idx_in,:);
        %G2_PW_nyr_unique = G2_PW_nyr_sort(idx_in,:);
        muisscale_unique = muisscale_sort(idx_in,:);
        %G2_PW_filen_unique = G2_PW_filen_sort(idx_in,:);
        % Distance filter:
%         lonlat_map_sort = horzcat(lonlat_map_sort,gtsrwaves_sort(:,5));
%         lonlat_map_unique = horzcat(lonlat_map_unique,gtsrwaves_unique(:,5));
%         idx_in = horzcat(idx_in,gtsrwaves_unique(:,5));
%         idx_out = horzcat(idx_out,gtsrwaves_sort(:,5));
%     gtsrwaves_sort(gtsrwaves_sort(:,5) > limit, :) = [];
%     gtsrwaves_unique(gtsrwaves_unique(:,5) > limit, :) = [];
%     muisscale_sort(muisscale_sort(:,2) > limit, :) = [];
%     muisscale_unique(muisscale_unique(:,2) > limit, :) = [];
%     lonlat_map_unique(lonlat_map_unique(:,3) > limit, :) = [];
%     lonlat_map_sort(lonlat_map_sort(:,3) > limit, :) = [];
%     idx_in(idx_in(:,2) > limit, :) = [];
%     idx_out(idx_out(:,2) > limit, :) = [];
    % save this sorted data because we will use it in the rest of this script
    save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_waves_150_data_clean_sorted.mat
    % save the tide gauge locations in a text file so that they can be read in SEAWISE  
    lonlat_map_unique = lonlat_map_unique(:,1:2);
    idx_in = idx_in(:,1);
    idx_out = idx_out(:,1);
     save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_waves_150_unique_clean.txt lonlat_map_unique -ASCII
     save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_waves_150_sort_clean.txt gtsrwaves_sort -ASCII
%  

%% Step 1: Load tide gauges & Gumbel values 
clear all
fprintf('----- Step 1: load TG data \n')
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_waves_150_data_clean_sorted.mat'); %this is the file created in step 0



%% Step 3: loop over PDF's, compute CDFs, save
fprintf('----- Step 3: compute CDFs \n')
gamma = 0.5772;

for station_nr = 1:length(gtsrwaves_unique) %loop over all TG stations

%point to YOUR location of seawise; 
SL_prob(:,1) = double(ncread(['/home/christos/seawise/svn/project.seawise/branches/seawise/seawise-test-1-',num2str(lonlat_map_unique(station_nr,1)),'-',num2str(lonlat_map_unique(station_nr,2)),'-r196-gtsrclean/combined-probability-distributions.nc'],'x-axis'));
SL_prob(:,2) = double(ncread(['/home/christos/seawise/svn/project.seawise/branches/seawise/seawise-test-1-',num2str(lonlat_map_unique(station_nr,1)),'-',num2str(lonlat_map_unique(station_nr,2)),'-r196-gtsrclean/combined-probability-distributions.nc'],'P_combined_3'));
SL_prob(:,3) = SL_prob(:,2)/sum(SL_prob(:,2));
% Compute the CDF from the PDF
SL_CDF = zeros(length(SL_prob),2);
for i = 1:length(SL_prob)
    SL_CDF(i,1) = sum(SL_prob(1:i,3));
    SL_CDF(i,2) = SL_prob(i,1);
end

% Bring down to list of 56 values only for use in fortran script (in step 4)
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_3/distrib_CDF.mat')
test(:,2) = 0;
for i = 2:length(test)
    I = find(round(SL_CDF(:,1),3) < test(i,1));
    test(i,2) = SL_CDF(max(I)+1,2); %takes first value that equals test(i,1)
end
    I = find(round(SL_CDF(:,1),3) == test(1,1));
    test(1,2) = SL_CDF(max(I),2);%takes last value that ==0;
%save CDF per TG station
save(['/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_3/CDF_normal_christos_',num2str(lonlat_map_unique(station_nr,1)),'-',num2str(lonlat_map_unique(station_nr,2)),'.txt'],'test','-ASCII')
clear i I SL_CDF SL_pro test SL_prob


%% Step4 run fortran-script to compute mix of Normal distributions & read into Matlab
% This script runs from MATLAB except for the compiling (ONE TIME only)
% open terminal and find Step_4 directory
% type: gfortran fmn1682.f -o fmn1682

mixt_nr = 1; %1 = normal distribution
%fprintf('----- Step 4: compute mixture of normals with %i normal for station %i \n',mixt_nr,station_nr)
unix(['cat /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_3/CDF_normal_christos_',num2str(lonlat_map_unique(station_nr,1)),'-',num2str(lonlat_map_unique(station_nr,2)),'.txt | ./Step_4/fmn1682 ',num2str(mixt_nr)]); % saves in file: 'Step_4/mixed_normals.txt'

filename = '/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_4/mixed_normals.txt';% load mixed_normals.txt %mean, std, weight
delimiter = ' ';
startRow = 9;
endRow = startRow+mixt_nr-1;
formatSpec = '%f%f%f%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
mixednormals = [dataArray{1:end-1}]; %mean, std, weight
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;

%store output per station
one_nrm(:,:,station_nr) = mixednormals;
clear mixednormals

%% Step5: compute Allowances for Normal Distributions of SLR
%fprintf('----- Step 5: compute allowances \n'); %this can be more than 1 station per grid point

% 1 normal distribution
mean_SLR(station_nr)        = one_nrm(:,1,station_nr);
std_SLR(station_nr)         = one_nrm(:,2,station_nr);

idx_in((length(gtsrwaves_unique(:,1))+1),1) = (length(gtsrwaves(:,1))+1); %to make sure the last tide gauge is included
for TG_stat= idx_in(station_nr):idx_in(station_nr+1)-1
fprintf('----- Step 5a: actual station nr %i \n',TG_stat)
gumbel_PW = muisscale_sort(TG_stat);
mean_waves(TG_stat) = wavesatgtsr_sort(TG_stat,3) + wavesatgtsr_sort(TG_stat,4).*gamma;
std_waves(TG_stat) = (pi./sqrt(6)).*wavesatgtsr_sort(TG_stat,4);
term1(TG_stat) = (mean_SLR(station_nr)+mean_waves(TG_stat));
term2(TG_stat) = ((std_SLR(station_nr).^2 + std_waves(TG_stat).^2)./(2.*gumbel_PW));
allow_normal_gtsr_waves(TG_stat,1) = (mean_SLR(station_nr)+mean_waves(TG_stat)) + ((std_SLR(station_nr).^2 + std_waves(TG_stat).^2)./(2.*gumbel_PW));% a = mean + (std^2 / 2*gumbel);


clear gumbel_* 
end %end loop over stations for each map value
clear i mixednormals mixt_nr
end %end loop over map values
%%
mean_waves = wavesatgtsr_sort(:,3) + wavesatgtsr_sort(:,4).*gamma;
std_waves = (pi./sqrt(6)).*wavesatgtsr_sort(:,4);
for k = 1:length(7629)
    mean_SLR_all(k,1) = mean_SLR(idx_out(k));
    std_SLR_all(k,1) = std_SLR(idx_out(k));
end

allowances = mean_SLR_all + mean_waves + ((std_SLR_all.^2 + std_waves.^2)./muisscale_sort(:,1));
%%
clear mix_nrm one_nrm station_nr

save allowances_ipcc_GTSR_waves_150_clean.mat