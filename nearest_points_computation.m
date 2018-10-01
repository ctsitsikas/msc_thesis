%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: nearest_points_computation.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% loads data and transposes arrays to columns if any
load('/home/christos/Dropbox/thesis/gtsr_at_waves.mat');
lon = muislon;
lat = muislat;
muisgtlonlat = horzcat(muislon, muislat);
muisgtscale = scalecomb';
loc = loccomb';
clear loccomb scalecomb
%max distance kept in meters
limit = 150000;

        %load land-ocean mask from sea-level fields
    temp = ncread('Step_1/sealevel-change-fields-at-2100-slangen.nc','SEtotalB');
    temp1 = ncread('Step_1/sealevel-change-fields-at-2100-slangen.nc','totalB');
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
    
    %lon and lat of SLR grid
    lon = ncread('Step_1/sealevel-change-fields-at-2100-slangen.nc','lon');
    lat = ncread('Step_1/sealevel-change-fields-at-2100-slangen.nc','lat');
    masklon = mask .* repmat(lon,[1 180]); 
    masklon(masklon == 0) = NaN;
    masklat = mask .* repmat(lat',[360 1]);
    masklat(masklat == 0) = NaN;
    clear lon lat
    a = 6371221;    % earth radius in meters
    d2r = 3.1415927/180;    %degrees to radians converter
 % now matching the tide gauge lat-lon to map locations
    for station_nr = 1:length(muisgtlonlat(:,1));
        totdiff = a.*acos(cos(d2r.*repmat(muisgtlonlat(station_nr,2),[360 180])).*cos(d2r.*masklat(:,:)).*cos(d2r.*(repmat(muisgtlonlat(station_nr,1),[360 180])-masklon(:,:)))+sin(d2r.*repmat(muisgtlonlat(station_nr,2),[360 180])).*sin(d2r.*masklat(:,:)));

        [num idx] = min(totdiff(:)); % find minimum distance
        [lonlat_map(station_nr,1) lonlat_map(station_nr,2)] = ind2sub(size(totdiff),idx); %write array location of minimum
        muisgtlonlat(station_nr,3) = masklon(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lon
        muisgtlonlat(station_nr,4) = masklat(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lat
        muisgtlonlat(station_nr,5) = num;   %keeps also the distance to apply the limit
        muisgtscale(station_nr,2) = num;    
        loc(station_nr,2) = num;
    end
   
    lonlat_map = horzcat(lonlat_map,muisgtlonlat(:,5));
    lonlat_map(lonlat_map(:,3) > limit, :) = [];
    lonlat_map = lonlat_map(:,1:2);
    muisgtlonlat(muisgtlonlat(:,5) > limit, :) = [];
    muisgtscale(muisgtscale(:,2) > limit, :) = [];
    loc(loc(:,2) > limit, :) = [];

        
        %% check locations in map to see if they match


        %Because there are sometime more tide gauges in the same location, we sort the lat-lons to reduce computing time later
        [lonlat_map_sort, index] = sortrows(lonlat_map,[2 1]);
        muisgtlonlat_sort = muisgtlonlat(index,:);
        
        muisgtscale_sort = muisgtscale(index,:);
        loc_sort = loc(index,:);
        
        clear i idx index j latdiff londiff num station_nr totdiff
        
        %FIND UNIQUE GRIDPOINTS FOR RUNNING seawise and reducing analysis
        [lonlat_map_unique, idx_in, idx_out] = unique(lonlat_map_sort,'rows','stable');
        muisgtlonlat_unique = muisgtlonlat_sort(idx_in,:);
        
        muisgtscale_unique = muisgtscale_sort(idx_in,:);
        

    % save this sorted data because we will use it in the rest of this script
    save Step_1/muisgt_data_sorted.mat
    % save the tide gauge locations in a text file so that they can be read in SEAWISE  
    lonlat_map_unique = lonlat_map_unique(:,1:2);
    idx_in = idx_in(:,1);
    idx_out = idx_out(:,1);
     save Step_1/muisgt_unique_100.txt lonlat_map_unique -ASCII
     save Step_1/muisgt_sort_100.txt muisgtlonlat_sort -ASCII