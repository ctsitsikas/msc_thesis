%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: subtract_distributions.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First part subtracts end-hist distributions, second part adds gtsr+waves
clear all

%load the wave Gumbel distributions
load('/home/christos/Dropbox/thesis/gumbelmaxendhist.mat', 'gumbelhist', 'gumbelend', 'lon','lat');

gamma = 0.5772; %Euler constant
xtest = -3:0.1:14; % values for x
countreal = 0;
countimag = 0;
for i = 1:length(lon)
    
    for j = 1:length(lat)
                

        loc1 = gumbelhist(1,j,i);    
        scale1 = gumbelhist(2,j,i);  
        loc = gumbelend(1,j,i);      
        scale = gumbelend(2,j,i);    
        meanP1 = loc1 + scale1*gamma; %mean and std of Gumbel 
        meanP = loc + scale*gamma;
        sigma1 = (pi./sqrt(6)).*scale1;
        sigma = (pi./sqrt(6)).*scale;
        
        a(j,i) = meanP1-6*sigma1;
        b(j,i) = meanP1+6*sigma1;
        
        % c = (mean-10*sigma)-a;
        % d = (mean+10*sigma)-b;
        
        e(j,i) = (meanP-6*sigma);
        f(j,i) = (meanP+6*sigma);
        

%         x1 = a:0.1:b; 
%         x = e:0.1:f; 
        P1 = PDF(xtest,loc1, scale1);
        P = PDF(xtest,loc, scale);
        %covariance if needed (in case two PDFs are correlated)
        %covPP1 = (1./length(P)).*sum((P-mean(P)).*(P1-mean(P1)));
        
        mean2 = meanP-meanP1;
        sigma2 = sqrt(sigma^2-sigma1^2);
        
        c(j,i) = (mean2-6*sigma2);
        d(j,i) = (mean2+6*sigma2);
        
        %computes Gumbel parameters of subtraction, if std_future<std_hist gives
        %0 values
        if isreal((sqrt(6)./pi).*sigma2) == 1
            scale2(j,i) = (sqrt(6)./pi).*sigma2;
            loc2(j,i) = mean2 - scale2(j,i).*gamma;
        else
            scale2(j,i) = 0;
            loc2(j,i) = 0;
        end
        
        clear P P1
    end
end

%% checking with seawise output
clear dist mu lambda mu1 lambda1 mu2 lambda2
locnames = ["San Fransisco", "New York", "Ijmuiden", "Bay of Bengal", "Kanmen, China", "Brest", "Mar del Plata", "Fremantle","Pago Pago"];
lons = [-122.5 -72.5 4.5 87.5 122.5 -4.5 -57.5 115.5 -170.5];
lats = [37.5 40.5 52.5 13.5 28.5 48.5 -38.5 -32.5 -14.5];
% swout  = xlsread('location1.xls',i);
% swout(swout == NaN) = 0;

fig = figure

p = uipanel('Parent',fig,'BorderType','none'); 
p.Title = 'End Century maxima distributions approximated by MATLAB in selected locations'; 
%p.Title= 'Gumbel PDF of historical and end century period at selected locations';
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
for i = 1:length(lons)
    idlat(i) = find(round(lat,1)==lats(i));
    idlon(i) = find(round(lon,1) == lons(i));
    file = sprintf('/home/christos/seawise/svn/project.seawise/branches/seawise/seawise-synthetic-no-fingerprint-r196-loc%d/fort.45',i);
    
    if i == 3
        idlon(i) = 185;
        idlat(i) = 137; 
    end
    
    sw = importdata(file,' ');
     
    swout = sw(:,1:4);
    swout(any(isnan(swout), 2), :) = [];
    
    idsw1 = find(round(swout(:,1),3)==round(e(idlat(i), idlon(i)),3));
    idsw2 = find(round(swout(:,1),3)==round(f(idlat(i), idlon(i)),3));
    
    swout = swout(idsw1:idsw2,:);
    
    
    
    
    
    
    subplot(3,3,i, 'Parent',p)
    %x = 1:0.01:10;

    mu2(i) = loc2(idlat(i),idlon(i));
    lambda2(i) = scale2(idlat(i),idlon(i));
    mu1(i) = gumbelhist(1,idlat(i),idlon(i));     %gumbelhist(1,50,50);
    lambda1(i) = gumbelhist(2,idlat(i),idlon(i));   %gumbelhist(2,50,50);
    mu(i) = gumbelend(1,idlat(i),idlon(i));      %muisdata(5103,3);
    lambda(i) = gumbelend(2,idlat(i),idlon(i)); 

    str = sprintf('point%s',string(i));
    dist1.(str)(:,1) = a(idlat(i), idlon(i)):0.001:b(idlat(i), idlon(i));
    dist2.(str)(:,1) = c(idlat(i), idlon(i)):0.001:d(idlat(i), idlon(i));
    dist.(str)(:,1) = e(idlat(i), idlon(i)):0.001:f(idlat(i), idlon(i));
    
    dist2.(str)(:,2) = PDF(dist2.(str)(:,1),mu2(i), lambda2(i));
    dist1.(str)(:,2) = PDF(dist1.(str)(:,1),mu1(i), lambda1(i));
    
    dist.(str)(:,2) = PDF(dist.(str)(:,1),mu(i), lambda(i));
    
    % normalizing
    dist.(str)(:,2) = dist.(str)(:,2)./max(dist.(str)(:,2)); 
    dist1.(str)(:,2) = dist1.(str)(:,2)./max(dist1.(str)(:,2));
    dist2.(str)(:,2) = dist2.(str)(:,2)./max(dist2.(str)(:,2));

    swout(:,2:4) = swout(:,2:4)./max(swout(:,2:4));
    
    swcdf = cumsum(swout(:,4),1);
    swcdf = swcdf./max(swcdf);
    matcdf = cumsum(dist.(str)(:,2),1);
    matcdf = matcdf./max(matcdf);
    
    prcsw = swout((find(round(swcdf,3)==0.950)),1);
    prcsw = prcsw(1);
    prcmat = swout((find(round(matcdf,3)==0.950)),1);
    prcmat = prcmat(1);
    
    
    plot(dist.(str)(:,1), dist.(str)(:,2),'Color','blue')
    

    ;hold on; plot(swout(:,1),swout(:,4),'Color','red')
       ;hold on; line([prcmat prcmat], [0 1],'Color','blue'); 
       ;hold on; line([prcsw prcsw], [0 1],'Color','red'); 

    ylim([0 1.1])
    xlabel('x (m)')
    ylabel('P(x)')


  
    str1 = sprintf('\b%s\nmatlab: 95prc = %0.3f \ndata: 95prc = %0.3f', locnames(i),prcmat,prcsw);
    title(str1)
end
%% plotting location of P2 (change distribution)
scale2(scale2 == 6999) = NaN;
pcolor(lon, lat, loc2)
shading('interp')
colorbar
caxis([-0.5 0.8])
title('End century and historical maxima distributions difference distribution P_{change} location parameter')
