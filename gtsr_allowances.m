%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: gtsr_allowances.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%cd Desktop/for_Roderik_from_Aimee/compute_allowances/
%%This script computes allowances for a NORMAL/GAUSSIAN distribution
fprintf('-----> Start \n')
%%
% gtsrlat = ncread('/home/christos/Desktop/THESIS/GTSRGUMBELparv3.nc','latitude');
% gtsrlon = ncread('/home/christos/Desktop/THESIS/GTSRGUMBELparv3.nc','longitude');
% loc = ncread('/home/christos/Desktop/THESIS/GTSRGUMBELparv3.nc','loc');
% gtsrscale = ncread('/home/christos/Desktop/THESIS/GTSRGUMBELparv3.nc','shape');
load('/home/christos/Dropbox/thesis/gtsrclean.mat');


gtsrlon = gtsr(:,1);
gtsrlat = gtsr(:,2);
loc = gtsr(:,3);
gtsrscale = gtsr(:,4);

% Step 0 (only once): determine grid points closest to TG location%   
%COMMENT THIS STEP WHEN IT IS DONE
%max distance kept

limit = 150000;
%This step finds the closest ocean grid point for each gtsr point

for i=1:length(gtsrlon)
    if gtsrlon(i) < 0
                gtsrlon(i) = gtsrlon(i) + 360;
    end
end



gtsrlonlat = horzcat(gtsrlon,gtsrlat);
% gtsrscale = gtsr(:,3);
% loc = gtsr(:,4);
%gtsrlonlat(:,2) = gtsrlonlat(:,2) -90; %correct latitudes to be similar to GESLA 
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
    for station_nr = 1:length(gtsrlonlat(:,1));
        %londiff(:,:) = abs(repmat(gtsrlonlat(station_nr,1),[360 180])-masklon(:,:)).*cosd(gtsrlonlat(station_nr,2)).*(a*(3.14159/180)); %compute longitude distances
        %latdiff(:,:) = abs(repmat(gtsrlonlat(station_nr,2),[360 180])-masklat(:,:)).*(a.*(3.14159/180)); %compute latitude distances
        totdiff = a.*acos(cos(d2r.*repmat(gtsrlonlat(station_nr,2),[360 180])).*cos(d2r.*masklat(:,:)).*cos(d2r.*(repmat(gtsrlonlat(station_nr,1),[360 180])-masklon(:,:)))+sin(d2r.*repmat(gtsrlonlat(station_nr,2),[360 180])).*sin(d2r.*masklat(:,:)));
        %totdiff = sqrt(londiff.^2 + latdiff.^2); %compute total distance
        [num idx] = min(totdiff(:)); % find minimum distance
        [lonlat_map(station_nr,1) lonlat_map(station_nr,2)] = ind2sub(size(totdiff),idx); %write array location of minimum
        gtsrlonlat(station_nr,3) = masklon(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lon
        gtsrlonlat(station_nr,4) = masklat(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lat
        gtsrlonlat(station_nr,5) = num;
        gtsrscale(station_nr,2) = num;
        loc(station_nr,2) = num;
        dist(station_nr) = num;
    end
   
    lonlat_map = horzcat(lonlat_map,gtsrlonlat(:,5));
    lonlat_map(lonlat_map(:,3) > limit, :) = [];
    lonlat_map = lonlat_map(:,1:2);
    gtsrlonlat(gtsrlonlat(:,5) > limit, :) = [];
    gtsrscale(gtsrscale(:,2) > limit, :) = [];
    loc(loc(:,2) > limit, :) = [];

        
        %% check locations in map to see if they match
        %surface(mask)
        %hold on
        %scatter(lonlat_map(:,2),lonlat_map(:,1),'og')
        %scatter(gtsrlonlat(:,2)+90,gtsrlonlat(:,1),'rx')
        %shading interp

        %Because there are sometime more gtsr points in the same location, we sort the lat-lons to reduce computing time later
        [lonlat_map_sort, index] = sortrows(lonlat_map,[2 1]);
        gtsrlonlat_sort = gtsrlonlat(index,:);

        gtsrscale_sort = gtsrscale(index,:);
        loc_sort = loc(index,:);

        clear i idx index j latdiff londiff num station_nr totdiff
        
        %FIND UNIQUE GRIDPOINTS FOR RUNNING seawise and reducing analysis
        [lonlat_map_unique, idx_in, idx_out] = unique(lonlat_map_sort,'rows','stable');
        gtsrlonlat_unique = gtsrlonlat_sort(idx_in,:);

        gtsrscale_unique = gtsrscale_sort(idx_in,:);
    % save this sorted data because we will use it in the rest of this script
    save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_data_150_clean_sorted.mat
    % save the gtsr locations in a text file so that they can be read in SEAWISE  
    lonlat_map_unique = lonlat_map_unique(:,1:2);
    idx_in = idx_in(:,1);
    idx_out = idx_out(:,1);
     save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_unique_150_clean.txt lonlat_map_unique -ASCII
     save /home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_sort_150_clean.txt gtsrlonlat_sort -ASCII
%  

%% Step 1: Load tide gauges & Gumbel values 
clear all
fprintf('----- Step 1: load TG data \n')
load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/Step_1/gtsr_data_150_clean_sorted.mat'); %this is the file created in step 0


%% Step 2: Run Seawise for EACH tide gauge location to get the PDF of sea-level change
% THIS IS NOT DONE IN MATLAB!
fprintf('----- {Step 2: make sure seawise has been run} \n')

%% Step 3: loop over PDF's, compute CDFs, save
fprintf('----- Step 3: compute CDFs \n')

for station_nr = 1:length(gtsrlonlat_unique) %loop over all TG stations
station_nr
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
fprintf('----- Step 4: compute mixture of normals with %i normal for station %i \n',mixt_nr,station_nr)
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
fprintf('----- Step 5: compute allowances \n'); %this can be more than 1 station per grid point

% 1 normal distribution
mean_SLR(station_nr)        = one_nrm(:,1,station_nr);
std_SLR(station_nr)         = one_nrm(:,2,station_nr);

idx_in((length(gtsrlonlat_unique(:,1))+1),1) = (length(gtsrlonlat(:,1))+1); %to make sure the last tide gauge is included
for TG_stat= idx_in(station_nr):idx_in(station_nr+1)-1
fprintf('----- Step 5a: actual station nr %i \n',TG_stat)
gumbel_PW = gtsrscale_sort(TG_stat);
allow_normal_gtsr100(TG_stat,1) = mean_SLR(station_nr) + ((std_SLR(station_nr).^2)./(2.*gumbel_PW));% a = mean + (std^2 / 2*gumbel);

clear gumbel_* 
end %end loop over stations for each map value
clear i mixednormals mixt_nr
end %end loop over map values
%%

clear mix_nrm one_nrm station_nr

save allowances_ipcc_GTSR_150_clean.mat

%% DONE
fprintf('<----- Finish \n')

