%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: downloadhist.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Downloads historical maxima from CSIRO
clear all
models = strings;
models = ["ACCESS1.0","BCC-CSM1.1","CNRM-CM5","GFDL-CM3","HadGEM2-ES","INMCM4","MIROC5","MRI-CGCM3"];
mods = {'ACC', 'BCC', 'CNRM', 'GFDL','HADG', 'INMCM', 'MIROC', 'MRI'};

dates = strings;

%creates dates as part of file names
for y = 1979:2005
    for m = 1:12
        ms = num2str(m, '%02i');
        
        dates = horzcat(dates,(strcat(string(y),ms)));
        
    end
end

dates=dates(2:end);

clear ms m y
%%
for j = 1  %length(models)
    
    models(j)
    for i = 1:24 %(length(dates))
        i
        
        
        url = sprintf('http://tds-mel.csiro.au/thredds/dodsC/Global_wave_projections/HISTORICAL/CMIP5/%s/ww3_outf_%s.nc',models(j),dates(i));
       
        a = double(ncread(url, 'Hs'));
        
        
        temp = max(a,[],3);
        temp2(:,:,i) = squeeze(temp(:,:,1));
        temp2(temp2==-999) = NaN;
        clear a
        
    end

    
maxes = zeros(360, 160, (length(temp2(1,1,:))/12));
for k = 1 : (length(temp2(1,1,:))/12)
    z1 = (k - 1) * 12 + 1;
    z2 = z1 + 11;
    maxes(:, :, k) = max(temp2(:,:, z1:z2), [], 3);
end
maxhist(j).anmax = maxes;

end
%     mid(j).gumbel = nan(2,160,360);
%
%     for lon = 1:360
%         for lat = 1:160
%             A=-squeeze(anmax(lon,lat,:));
%             A(isnan(A))=[];
%             [gumbel(:,lat,lon)] = evfit(A);
%         end
%     end
%     mid(j).gumbel(1,:,:) = -mid(j).gumbel(1,:,:);
clear temp temp2


save hemer_data_hist.mat
