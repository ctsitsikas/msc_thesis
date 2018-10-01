%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: gtsr_gesla2_comparison.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%load data to find nearest points
load('/home/christos/Dropbox/thesis/gtsrclean.mat');
gesla = xlsread('/home/christos/Desktop/GESLA2_PW','A1:H658');
addpath('/home/christos/csirolib/');
%%

gtsrlonlat = gtsr(:,1:2);
gtsrscale = gtsr(:,4);
gtsrloc = gtsr(:,3);


%gesla2sort = sortrows(gesla2,2);
G2_PW_lonlat = gesla(:,1:2);
G2_PW_scale = gesla(:,3);
G2_PW_loc = gesla(:,6);
limit = 100000;
for i=1:length(gtsrlonlat(:,1))
    if gtsrlonlat(i,1) < 0
                gtsrlonlat(i,1) = gtsrlonlat(i,1) + 360;
    end
end
a = 6371221;
    d2r = 3.1415927/180;
for i = 1:length(G2_PW_lonlat(:,1))
    
        replon = transpose(repmat(G2_PW_lonlat(i,1),1,length(gtsrlonlat(:,1))));
        replat = transpose(repmat(G2_PW_lonlat(i,2),1,length(gtsrlonlat(:,1))));
        %londiff = abs(replon-gtsrlonlat(:,1)).*cosd(G2_PW_lonlat(i,2)).*(a.*(3.14159/180)); %compute longitude distances
        %latdiff = abs(replat-gtsrlonlat(:,2)).*(a.*(3.14159/180));
        %totdiff = sqrt(londiff.^2 + latdiff.^2);
        totdiff = a.*acos(cos(d2r.*replat).*cos(d2r.*gtsrlonlat(:,2)).*cos(d2r.*(replon-gtsrlonlat(:,1)))+sin(d2r.*replat).*sin(d2r.*gtsrlonlat(:,2)));
        [num idx] = min(totdiff(:)); % find minimum distance
        %[lonlat_map(i,1) lonlat_map(i,2)] = ind2sub(size(totdiff),idx);
        G2_PW_lonlat(i,3) = gtsrlonlat(idx,1);%write closest lon
        G2_PW_lonlat(i,4) = gtsrlonlat(idx,2);%write closest lat
        %G2_PW_lonlat(i,5) = num;
        G2_PW_scale(i,2) = gtsrscale(idx);
        G2_PW_loc(i,2) = gtsrloc(idx);
        distance(i) = num;
end
tidesgtsr = horzcat(G2_PW_lonlat,G2_PW_scale,G2_PW_loc,distance');
%tidesgtsr(tidesgtsr(:,9) > limit, :) = [];
% G2_PW_lonlat(G2_PW_lonlat(:,5) > limit, :) = [];
% G2_PW_scale(G2_PW_scale(:,3) > limit, :) = [];
% G2_PW_scale(G2_PW_scale(:,3) > limit, :) = [];
%tidesgtsr = horzcat(G2_PW_lonlat(:,1:4),G2_PW_scale,G2_PW_loc);
tidesgtsr = tidesgtsr(~any(isnan(tidesgtsr),2),:);

tidesgtsrsort = sortrows(tidesgtsr,2);
G2_PW_lonlat_sort = tidesgtsrsort(:,1:2);
G2_PW_scale_sort = tidesgtsrsort(:,5);
G2_PW_loc_sort = tidesgtsrsort(:,7);
%%
figure
coast('k') ;hold on;
scatter(gtsrlonlat(:,1),gtsrlonlat(:,2),20, gtsrloc,'filled') 
title('GTSR location parameter')
colorbar
caxis([0 3])
xlabel('Longitude')
ylabel('Latitude')
figure
coast('k') ;hold on;
scatter(G2_PW_lonlat_sort(:,1),G2_PW_lonlat_sort(:,2),20, G2_PW_loc_sort,'filled')
colorbar
title('GESLA2 location parameter')
caxis([0 3])
xlabel('Longitude')
ylabel('Latitude')
%%
coast('k') ;hold on;
scatter(tidesgtsr(:,3),tidesgtsr(:,4),'red') ;hold on;
scatter(tidesgtsr(:,1),tidesgtsr(:,2),'x','blue')

%%
scalediff = tidesgtsr(:,5)-tidesgtsr(:,6);

locdiff = tidesgtsr(:,7)-tidesgtsr(:,8);
ratioscale = tidesgtsr(:,5)./tidesgtsr(:,6);
ratioloc = tidesgtsr(:,7)./tidesgtsr(:,8);
rmslocnew = sqrt((sum((tidesgtsr(:,7)-tidesgtsr(:,8)).^2))/length(tidesgtsr(:,1)));
rmsesc = sqrt(immse(tidesgtsr(:,6),tidesgtsr(:,5)));
rmseloc = sqrt(immse(tidesgtsr(:,8),tidesgtsr(:,7)));
meanloc = mean(tidesgtsr(:,7));
meansc = mean(tidesgtsr(:,5));
scale = meansc;
loc = meanloc;
deltaloc = rmseloc;
deltasc = rmsesc;

figure(1)
x = [0:0.01:4];
plot(x,(1./meansc).*exp(-(((x-meanloc)./meansc)+exp(-(x-meanloc)./meansc))),'Linewidth',2)
;hold on;

%figure(2)

plot(x,(1./meansc).*exp(-(((x-(meanloc+deltaloc))./meansc)+exp(-(x-(meanloc+deltaloc))./meansc))),'Linewidth',2)

;hold on;


plot(x,(1./meansc).*exp(-(((x-(meanloc-deltaloc))./meansc)+exp(-(x-(meanloc-deltaloc))./meansc))),'Linewidth',2)


;hold on;

%plot(x,(1./meansc).*exp(-(((x-meanloc)./meansc)+exp(-(x-meanloc)./meansc))),'Linewidth',2)


;hold on;

plot(x,(1./(meansc+deltasc)).*exp(-(((x-meanloc)./(meansc+deltasc))+exp(-(x-meanloc)./(meansc+deltasc)))),'Linewidth',2) ;hold on;
plot(x,(1./(meansc-deltasc)).*exp(-(((x-meanloc)./(meansc-deltasc))+exp(-(x-meanloc)./(meansc-deltasc)))),'Linewidth',2)
%title('$\bar{\mu_s}+\delta_{loc}$,$\bar{\lambda_s}$','interpreter','latex', 'FontSize', 14)
h=legend('$\bar{\mu_s}, \bar{\lambda_s}$','$\bar{\mu_s}+\delta_{loc}$, $\bar{\lambda_s}$','$\bar{\mu_s}-\delta_{loc}$, $\bar{\lambda_s}$','$\bar{\mu_s}$, $\bar{\lambda_s}+\delta_{sc}$','$\bar{\mu_s}$, $\bar{\lambda_s}-\delta_{sc}$')
set(h,'Interpreter','latex','FontSize',12);
title('Gumbel PDF of mean GESLA parameters with RMSE of GESLA-gtsr')
%%
figure(1)
coast('k') ;hold on;
%draw_land(1, 'k','crd');
scatter(tidesgtsr(:,1),tidesgtsr(:,2),20, ratioloc,'filled')
colorbar
caxis([0.5 1.5])
%;hold on;
%scatter(G2_PW_lonlat(:,3),G2_PW_lonlat(:,4),'blue')
xlim([0 360])
ylim([-90 90])
xlabel('Longitude')
ylabel('Latitude')
title('location_{GESLA2}/location_{GTSR}')
%legend('Tide gauges','gtsr');
%%
figure(2)
coast('k') ;hold on;
%draw_land(1, 'k','crd');
scatter(tidesgtsr(:,1),tidesgtsr(:,2),20, ratioscale,'filled')
colorbar
caxis([0 4])
%;hold on;
%scatter(G2_PW_lonlat(:,3),G2_PW_lonlat(:,4),'blue')
xlim([0 360])
ylim([-90 90])
xlabel('Longitude')
ylabel('Latitude')
title('scale_{GESLA2}/scale_{GTSR}')
%legend('Tide gauges','gtsr');
%%
coast('k') ;hold on;
%draw_land(1, 'k','crd');
scatter(tidesgtsr(:,1),tidesgtsr(:,2),20, locdiff,'filled')
colorbar
caxis([-1 1])
%;hold on;
%scatter(G2_PW_lonlat(:,3),G2_PW_lonlat(:,4),'blue')
xlim([0 360])
ylim([-90 90])
xlabel('Longitude')
ylabel('Latitude')
title('location_{GESLA2}-location_{GTSR}')
%legend('Tide gauges','gtsr');
%%
%tidesgtsr(tidesgtsr(:,5) > 2, :) = [];
scatter(tidesgtsr(:,7), tidesgtsr(:,8),'filled','blue') ;hold on;
xlabel('GESLA2 location parameter')
ylabel('GTSR location parameter')
%xlim([0 8])
%ylim([0 8])
P = polyfit(tidesgtsr(:,7),tidesgtsr(:,8),1);
yfit = P(1)*tidesgtsr(:,7)+P(2);
hold on;
plot(tidesgtsr(:,7),yfit,'r-.');
str = sprintf('linear fit y=ax+b with a=%f and b=%f',P(1),P(2))
title(str)
%%
%tidesgtsr(tidesgtsr(:,5) > 2, :) = [];
scatter(tidesgtsr(:,5), tidesgtsr(:,6),'filled','blue') ;hold on;
xlabel('GESLA2 scale parameter')
ylabel('GTSR scale parameter')
%xlim([0 0.6])
%ylim([0 0.6])
P = polyfit(tidesgtsr(:,5),tidesgtsr(:,6),1);
yfit = P(1)*tidesgtsr(:,5)+P(2);
hold on;
plot(tidesgtsr(:,5),yfit,'r-.');
str = sprintf('linear fit y=ax+b with a=%f and b=%f',P(1),P(2))
title(str)
%%

scatter(tidesgtsr(:,2),(tidesgtsr(:,5)./tidesgtsr(:,6)), 'blue','filled')
xlabel('Latitude')
ylabel('$scale_{GESLA2}/scale_{GTSR}$','interpreter','latex', 'FontSize', 14)
%ylim([0 3.5])

