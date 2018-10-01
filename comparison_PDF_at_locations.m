%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: comparison_PDF_at_locations.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
addpath('/home/christos/csirolib/');

% Default options
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')

%load data
load('/home/christos/Dropbox/thesis/compute_allowances/gtsrdata.mat');
load('/home/christos/Dropbox/thesis/gesla.mat');

%change format of longitude
for i=1:length(gtsr(:,1))
    if gtsr(i,1) < 0
                gtsr(i,1) = gtsr(i,1) + 360;
    end
end

sanfr = [237.6000 37.8000];
ny = [286.1000 40.5000];
ijmuid = [4.6000 52.5000];
bengal = [87 20.7000];
kanmen = [121.7000 28.3000];
brest = [355.3000 48.1000];
plata = [302.8000 -37.6000];
freman = [115.8000 -31.9000];
pago = [189.3000 -14.3000];
ipcclocs = [sanfr; ny; ijmuid; bengal; kanmen; brest; plata; freman; pago];
idgesla = [390 420 547 236 271 499 43 74 134]; %idx of locations at gesla2 dataset(found manually)
locnames = ["San Fransisco", "New York", "Ijmuiden", "Bay of Bengal", "Kanmen, China", "Le Conquet", "Mar del Plata", "Fremantle","Pago Pago"];
%round data coordinates to compare with locations
muis1 = round(gtsr(:,1),1);
muis2 = round(gtsr(:,2),1);

for i = 1:length(ipcclocs(:,1))
    
    ids = find(muis1==ipcclocs(i,1) & muis2==ipcclocs(i,2),1);
    
    ipccids(i) = ids;
end


f = figure
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Gumbel Probability Density Function of GTSR and GESLA-2'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 16;
p.FontWeight = 'bold';
x = 0:0.01:4.5;
ipccids(7) = 404;
for i = 1:9
    
    subplot(3,3,i, 'Parent',p)
    
    plot(x,PDF(x,gtsr(ipccids(i),3),gtsr(ipccids(i),4))./max(PDF(x,gtsr(ipccids(i),3),gtsr(ipccids(i),4))), 'Linewidth', 2, 'Color', 'g')   ;hold on;        % line plot
    plot(x,PDF(x,gesla(idgesla(i),3),gesla(idgesla(i),4))./max(PDF(x,gesla(idgesla(i),3),gesla(idgesla(i),4))), 'Linewidth', 2, 'Color', 'r')

    str = sprintf(locnames(i))
    title(str)
    if i==2
    legend('GTSR PDF','GESLA2 PDF')
    end
    scale(i,1)= gesla(idgesla(i),4);
    loc(i,1)=gesla(idgesla(i),3);
    scale(i,2)= gtsr(ipccids(i),4);
    loc(i,2)=gtsr(ipccids(i),3);
     
    

if i==1
    xlabel('water height (m)')
    ylabel('Probability')
end
xlim([0 4.5])
ylim([0 1.2])
end

%% plot the nine selected locations on map

coast('k') ;hold on;
colors = [1 0 0 ; 0 1 0 ; 0 0 1 ; 0.4 0.4 0.4 ; 1 1 0 ; 1 0 1 ; 0 0 0.5 ; 0 1 1 ; 0.9 0.5 0];
for i = 1:9
    sc(i) = scatter(gesla(idgesla(i),1),gesla(idgesla(i),2),50, colors(i,:),'filled') ;hold on;
    colormap(jet);
end
str = cellstr(locnames);
legend([sc],str);

ylim([-90 90])
xlim([0 360])