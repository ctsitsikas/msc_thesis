%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: wave_data_gumbel_models.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%load maxima of models and coordinates from random model data file 
load('C:/Users/christos/Dropbox/thesis/hemer_data_end.mat');
load('C:/Users/christos/Dropbox/thesis/hemer_data_hist.mat');
load('C:/Users/christos/Desktop/THESIS/data/MID21c/RCP85/MIROC5/ww3_monthly_analysis.mat', 'lon');
load('C:/Users/christos/Desktop/THESIS/data/MID21c/RCP85/MIROC5/ww3_monthly_analysis.mat', 'lat');

%% check models seperately

for model = 1:8
    model
    endc(model).gumbel = nan(2,160,360); %gumbel fit parameters
    enderr(model).confgumbel = nan(4,160,360); %confidence of gumbel fit
    for lo = 1:length(lon)
        for la = 1:length(lat)
            
            
            A=-maxend(model).anmax(lo,la,:);
            
            A(isnan(A))=0;
            
            A = squeeze(A(1,1,:));
            
            %Hs_m(isnan(Hs_m,1),limnanlat(n),limnanlon(n))=[];
            [endc(model).gumbel(:,la,lo) parmci] = evfit(A);
            enderr(model).confgumbel(:,la,lo) = reshape(parmci,[1, 4]);
            
        end
    end
    endc(model).gumbel(1,:,:) = -endc(model).gumbel(1,:,:); %loc parameter has to be - for maxima distribution
end

%% check models seperately - hist

for model = 1:8
    model
hist(model).gumbel = nan(2,160,360);
    histerr(model).confgumbel = nan(4,160,360);
    for lo = 1:length(lon)
        for la = 1:length(lat)
            
            
            A=-maxhist(model).anmax(lo,la,:);
            
            A(isnan(A))=0;
            
            A = squeeze(A(1,1,:));
            
            %Hs_m(isnan(Hs_m,1),limnanlat(n),limnanlon(n))=[];
            [hist(model).gumbel(:,la,lo) parmci] = evfit(A);
            histerr(model).confgumbel(:,la,lo) = reshape(parmci,[1, 4]);
            
        end
    end
hist(model).gumbel(1,:,:) = -hist(model).gumbel(1,:,:);
end

%% plot models
models = ["ACCESS1","BCCCSM1","CNRMCM5","GFDLCM3","HADGEM2ES","INMCM4","MIROC5","MRICGCM3"];
f = figure
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Scale parameter of end century maxima Gumbel fit'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 18;
p.FontWeight = 'bold';

for i = 1:8
    temp = endc(i).gumbel(2,:,:);
    temp(temp==0)=NaN;
    subplot(3,3,i, 'Parent',p)
    colormap(jet)
    pcolor(lon,lat,squeeze(squeeze(temp(1,:,:))))          
    
    shading('interp')
    
    colorbar
    
    caxis([0 2.5])
    str = sprintf(models(i));
    title(str,'FontSize',18);
    if i == 1
        xlabel('Longitude', 'Fontsize',16)
        ylabel('Latitude', 'Fontsize',16)
    end
    
end


%% plot confidence of fit as difference of range over central value
models = ["ACCESS1","BCCCSM1","CNRMCM5","GFDLCM3","HADGEM2ES","INMCM4","MIROC5","MRICGCM3"];
f = figure
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Scale parameter 95% confidence percentage of end century maxima gumbel fit'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 16;
p.FontWeight = 'bold';

for i = 1:8
    temp = endc(i).gumbel(2,:,:);
    temper = abs(enderr(i).confgumbel(1,:,:)-enderr(i).confgumbel(2,:,:));
    
    temp(temp==0)=NaN;
    temper(temper==0)=NaN;
    error = temper./temp.*100;
    subplot(3,3,i, 'Parent',p)
    colormap(jet)
    pcolor(lon,lat,squeeze(squeeze(error(1,:,:))))          
    
    shading('interp')
    
    colorbar
    
    caxis([0 30])
    str = sprintf(models(i));
    title(str, 'Fontsize', 14);
    if i == 1
        xlabel('Longitude', 'Fontsize',14)
        ylabel('Latitude', 'Fontsize',14)
    end
end

%% change temp for end/hist ensemble

sum19 = zeros(360,160,19);
sum20 = zeros(360,160);
for i = 1:8
    
    tempend = maxend(i).anmax;
    
    sum19 = sum19(:,:,1:19) + tempend(:,:,1:19);
    
    if i ~= 2 && i ~= 5
           sum20 = sum20 + tempend(:,:,20);
    end
    
    clear temp19
end
avg19 = sum19./8;
avg20 = sum20./6;

ensemblemean_end = cat(3,avg19,avg20);

sum27 = zeros(360,160,27);
for i = 1:8
    temphist = maxhist(i).anmax;
    sum27 = sum27 + temphist;
    
    clear temp
end

ensemblemean_hist = sum27./8;
%% GUMBEL WITH NAN=0 for hist

gumbelhist = nan(2,160,360);
ensemblemean_hist(isnan(ensemblemean_hist)) = 0 ;
for i = 1:length(lon)
    for j = 1:length(lat)
        A=-ensemblemean_hist(i,j,:);
        A(isnan(A))=[];
        
        [gumbelhist(:,j,i) parmci] = evfit(squeeze(A(1,1,:)));
        confgumbelhist(:,j,i) = reshape(parmci,[1, 4]);
    end
end
gumbelhist(1,:,:) = -gumbelhist(1,:,:);
gumbelhist(gumbelhist == 0) = NaN;
confgumbelhist(confgumbelhist == 0) = NaN;
%% GUMBEL WITH NAN=0 for end
clear A
gumbelend = nan(2,160,360);
ensemblemean_end(isnan(ensemblemean_end)) = 0 ;
for i = 1:length(lon)
    for j = 1:length(lat)
        A=-ensemblemean_end(i,j,:);
        A(isnan(A))=[];
        
        [gumbelend(:,j,i) parmci] = evfit(squeeze(A(1,1,:)));
        confgumbelend(:,j,i) = reshape(parmci,[1, 4]);
    end
end
gumbelend(1,:,:) = -gumbelend(1,:,:);
gumbelend(gumbelend == 0) = NaN;
confgumbelend(confgumbelend == 0) = NaN;
%%
figure
pcolor(lon,lat,squeeze(gumbelhist(1,:,:)))
shading('interp')
title('Location parameter of historical maxima ensemble mean (m)','Fontsize',18)
colorbar
caxis([0 14])
xlabel('Longitude', 'Fontsize',16)
ylabel('Latitude', 'Fontsize',16)
figure
pcolor(lon,lat ,squeeze(gumbelhist(2,:,:)))
shading('interp')
title('Scale parameter of historical maxima ensemble mean','Fontsize',18)
colorbar
caxis([0 1])
xlabel('Longitude', 'Fontsize',16)
ylabel('Latitude', 'Fontsize',16)
%%
figure
temploc = gumbelhist(1,:,:);
temperloc = abs(confgumbelhist(1,:,:)-confgumbelhist(2,:,:));
temploc(temploc==0)=NaN;
temperloc(temperloc==0)=NaN;
errorloc = temperloc./temploc.*100;
pcolor(lon,lat ,squeeze(errorloc(1,:,:)))
shading('interp');
title('Location parameter 95% confidence percentage of historical maxima ensemble mean');
colorbar;
caxis([0 15])

figure
tempscale = gumbelhist(2,:,:);
temperscale = abs(confgumbelhist(3,:,:)-confgumbelhist(4,:,:));
tempscale(tempscale==0)=NaN;
temperscale(temperscale==0)=NaN;
errorscale = temperscale./tempscale.*100;
pcolor(lon,lat ,squeeze(errorscale(1,:,:)))
shading('interp');
title('Scale parameter 95% confidence percentage of historical maxima ensemble mean')
colorbar;
caxis([40 70])

%%
figure

% Create colormap that is green for negative, red for positive,
% and a chunk inthe middle that is black.
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
% Apply the colormap.
colormap(colorMap);

locdiff = gumbelend(1,:,:)-gumbelhist(1,:,:);
pcolor(lon,lat ,squeeze(locdiff(1,:,:)))
shading('interp');
title('Location parameter difference between end century and historical period for ensemble mean (m)','Fontsize', 18);
c = colorbar;
caxis([-0.8 0.8])
c.FontSize = 12;
xlabel('Longitude', 'Fontsize',16)
ylabel('Latitude', 'Fontsize',16)
%%
clear all
%load the ensemble values, saved from upper part
load('gumbelmaxendhist.mat');

locnames = ["San Fransisco", "New York", "Ijmuiden", "Bay of Bengal", "Kanmen, China", "Le Conquet", "Mar del Plata", "Fremantle","Pago Pago"];
lons = [-122.5 -72.5 4.5 87.5 122.5 -4.5 -57.5 -115.5 -170.5];
lats = [37.5 40.5 52.5 13.5 28.5 48.5 -38.5 -32.5 -14.5];
f = figure
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Gumbel PDF with parameters and maxima histogram for end century in selected locations'; 
%p.Title= 'Gumbel PDF of historical and end century period at selected locations';
p.TitlePosition = 'centertop'; 
p.FontSize = 16;
p.FontWeight = 'bold';
for i = 1:length(lons)
    
    subplot(3,3,i, 'Parent',p)
    x = 1:0.01:10;
    idlat(i) = find(round(lat,1)==lats(i),1);
    idlon(i) = find(round(lon,1) == lons(i),1);
    mu(i) = gumbelend(1,idlat(i),idlon(i));
    lambda(i) = gumbelend(2,idlat(i),idlon(i));

    
    % historical

    muh(i) = gumbelhist(1,idlat(i),idlon(i));
    lambdah(i) = gumbelhist(2,idlat(i),idlon(i));
    muupper(i) = -confgumbelhist(2,idlat(i),idlon(i));
	lambdaupper(i) = confgumbelhist(4,idlat(i),idlon(i));
	mulower(i) = -confgumbelhist(1,idlat(i),idlon(i));
	lambdalower(i) = confgumbelhist(3,idlat(i),idlon(i));
    
    
    plot(x,PDF(x,muh(i),lambdah(i)), 'Linewidth', 2, 'Color', 'r')  ;hold on;
    plot(x,PDF(x,mu(i),lambda(i)), 'Linewidth', 2, 'Color', 'b') ;hold on;
    
%     plot(x,PDF(x,muupper(i),lambda(i))) ;hold on;
%     plot(x,PDF(x,mulower(i),lambda(i))) ;hold on;
%     plot(x,PDF(x,mu(i),lambdaupper(i))) ;hold on;
%     plot(x,PDF(x,mu(i),lambdalower(i))) ;hold on;
    if i == 1

        ylabel('P(x)')
        
    end
    %yyaxis right
    %histogram(ensemblemean_end(idlon(i),idlat(i),:),5, 'FaceAlpha', 0.1)

    xlim([1 10])
    %ylim([0 4.5])
    
    
    str = sprintf("%s", locnames(i))
    
    str1 = sprintf("location = %4.2f m", mu(i));
    str2 = sprintf("scale = %4.2f m", lambda(i));
    annotation('textbox','String',{str1;str2});
    title(str)
    grid on;
%     if i == 2
%         legend('\mu, \lambda','\mu_{upper}, \lambda','\mu_{lower}, \lambda','\mu, \lambda_{upper}','\mu, \lambda_{lower}')
%         %legend('Historical','End Century');
%     end
    %if i == 1
     %   xlabel('x (m)')
      %  ylabel('counts')
        
    end
    %lgnd = legend(str1);
     %   set(lgnd,'color','none');
    %end
end
