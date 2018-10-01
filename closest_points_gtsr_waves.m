%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: closest_points_gtsr_waves.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
load('/home/christos/Dropbox/thesis/P2.mat');
load('/home/christos/Dropbox/thesis/gtsrclean.mat');


%max distance kept
limit = 100000;

for i=1:length(gtsr(:,1))
    if gtsr(i,1) < 0
                gtsr(i,1) = gtsr(i,1) + 360;
    end
end

for i=1:length(lon)
    if lon(i) < 0
                lon(i) = lon(i) + 360;
    end
end





gtsrlonlat = gtsr(:,1:2);
gtsrscale = gtsr(:,4);
gtsrloc = gtsr(:,3);
%load waves grid
    
temp = loc2;
temp1 = scale2;
    for i=1:360
        for j =1:160
            if temp(j,i)==NaN || temp1(j,i)==NaN
                mask_scale(j,i) = NaN;
            else
                mask_scale(j,i) = temp1(j,i);
                mask_loc(j,i) = temp(j,i);
            end
        end
    end
     mask_scale = mask_scale';
     mask_loc = mask_loc';

    mask = ~isnan(mask_scale);

    masklon = mask .* repmat(lon,[1 160]); 
    masklon(masklon == 0) = NaN;
    masklat = mask .* repmat(lat',[360 1]);
    masklat(masklat == 0) = NaN;
    
    
    a = 6371221;
    d2r = 3.1415927/180;

    
    % now matching the tide gauge lat-lon to map locations
    for station_nr = 1:length(gtsrlonlat(:,1));
        totdiff = a.*acos(cos(d2r.*repmat(gtsrlonlat(station_nr,2),[360 160])).*cos(d2r.*masklat(:,:)).*cos(d2r.*(repmat(gtsrlonlat(station_nr,1),[360 160])-masklon(:,:)))+sin(d2r.*repmat(gtsrlonlat(station_nr,2),[360 160])).*sin(d2r.*masklat(:,:)));
        [num idx] = min(totdiff(:)); % find minimum distance
        [lonlat_map(station_nr,1) lonlat_map(station_nr,2)] = ind2sub(size(totdiff),idx); %write array location of minimum
        gtsrlonlat(station_nr,3) = masklon(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lon
        gtsrlonlat(station_nr,4) = masklat(lonlat_map(station_nr,1),lonlat_map(station_nr,2));%write closest lat
        gtsrlonlat(station_nr,5) = num;
        gtsrscale(station_nr,2) = mask_scale(lonlat_map(station_nr,1),lonlat_map(station_nr,2));
        gtsrloc(station_nr,2) = mask_loc(lonlat_map(station_nr,1),lonlat_map(station_nr,2));
        gtsrscale(station_nr,3) = num;
        gtsrloc(station_nr,3) = num;
    end
   
    lonlat_map = horzcat(lonlat_map,gtsrlonlat(:,5));
    lonlat_map(lonlat_map(:,3) > limit, :) = [];
    lonlat_map = lonlat_map(:,1:2);
    gtsrlonlat(gtsrlonlat(:,5) > limit, :) = [];
    gtsrscale(gtsrscale(:,3) > limit, :) = [];
    gtsrloc(gtsrloc(:,3) > limit, :) = [];
    
    clear num
    %% save nearest points of Delta-wave and GTSR
    %clear gtsratnearestwaves
    gtsratnearestwaves_clean = horzcat(gtsrlonlat(:,1:2),gtsrloc(:,1),gtsrscale(:,1));
    wavesatgtsr = horzcat(gtsrlonlat(:,1:2),gtsrloc(:,2),gtsrscale(:,2));
    save gtsr_waves_near_clean.mat gtsratnearestwaves_clean wavesatgtsr