%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: slr_gtsr_waves_at_ipcc_locations.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot all allowances and SLC PDF for nine locations
clear all
% Default options for fontsize
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultTextFontSize',16)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')

%load GTSR+waves allowances
IPCC_PW = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_waves_150_clean.mat');

%load GTSR allowances
IPCC = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_GTSR_150_clean.mat');

%load GESLA2 allowances
load('/home/christos/Desktop/for_Roderik_from_Aimee/figure_plot_files/nanidsort.mat');
IPCC_gesla = load('/home/christos/Desktop/for_Roderik_from_Aimee/compute_allowances/allowances_ipcc_9nov.mat');
IPCC_gesla.G2_PW_lonlat_sort(nanidxsort,:)=[];
IPCC_gesla.allow_normal_PW(nanidxsort,:)=[];


%% remove GTSR+waves allowances higher than 10 meters
allowwaves = horzcat(IPCC_PW.gtsrwaves_sort(:,1),IPCC_PW.gtsrwaves_sort(:,2),IPCC_PW.allow_normal_gtsr_waves,IPCC_PW.muisscale_sort(:,1));
clear IPCC_PW
over = find(allowwaves(:,3)>10);
allowwaves(over,:)=[];
IPCC_PW.gtsrwaves_sort(:,1) = allowwaves(:,1);
IPCC_PW.gtsrwaves_sort(:,2) = allowwaves(:,2);
IPCC_PW.allow_normal_gtsr_waves = allowwaves(:,3);
IPCC_PW.muisscale_sort = allowwaves(:,4);

%%
% coordinates of nine locations
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
ipcclocsg = [387 422 543 231 263 511 34 70 133]; %idx of locations for GESLA2 allowances (found manually)

locnames = ["San Fransisco", "New York", "Ijmuiden", "Bay of Bengal", "Kanmen, China", "Le Conquet", "Mar del Plata", "Fremantle","Pago Pago"];

%%  round data coordinates to compare with locations
    muis1 = round(IPCC.muislonlat_sort(:,1),1);
    muis2 = round(IPCC.muislonlat_sort(:,2),1);
    waves1 = round(IPCC_PW.gtsrwaves_sort(:,1),1);
    waves2 = round(IPCC_PW.gtsrwaves_sort(:,2),1);

for i = 1:length(ipcclocs(:,1)) %find idx of allowances at locations

    ids = find(muis1==ipcclocs(i,1) & muis2==ipcclocs(i,2),1);
    ids2 = find(waves1==ipcclocs(i,1) & waves2==ipcclocs(i,2),1);

    ipccids(i) = ids;
    ipccidsg(i) = ids2;

    %read SLC PDF from seawise output
    x(:,i) = ncread(['/home/christos/seawise/svn/project.seawise/branches/seawise/seawise-test-1-',num2str(IPCC.lonlat_map_sort(ids,1)),'-',num2str(IPCC.lonlat_map_sort(ids,2)),'-r196-gtsrclean/combined-probability-distributions.nc'],'x-axis');

    pcomb(:,i) = ncread(['/home/christos/seawise/svn/project.seawise/branches/seawise/seawise-test-1-',num2str(IPCC.lonlat_map_sort(ids,1)),'-',num2str(IPCC.lonlat_map_sort(ids,2)),'-r196-gtsrclean/combined-probability-distributions.nc'],'P_combined_3');
end

%%
allowgesla = IPCC_gesla.allow_normal_PW;
allowwave = IPCC_PW.allow_normal_gtsr_waves;
allowgtsr = IPCC.allow_normal_muis150(:,1);
f = figure
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'SLR PDF with allowances from GESLA-2, GTSR and GTSR combined with wave maxima change'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 16;
p.FontWeight = 'bold';
for i = 1:9
    
    subplot(3,3,i, 'Parent',p)
    plot(x(:,i),pcomb(:,i), 'Linewidth', 2)           % SLC PDF plot
    xlim([-0.5 2.5])
    ylim([0 2*10^(-3)])
    str = sprintf(locnames(i))
    title(str,'Fontsize',18)
    %allowances lines plots
    line([allowgesla(ipcclocsg(i)) allowgesla(ipcclocsg(i))],[0 2*10^(-3)],'Color','red', 'Linewidth', 2)
    line([allowgtsr(ipccids(i)) allowgtsr(ipccids(i))],[0 2*10^(-3)],'Color','green', 'Linewidth', 2)
    line([allowwave(ipccidsg(i)) allowwave(ipccidsg(i))],[0 2*10^(-3)],'Color','black', 'Linewidth', 2)
    xt = get(gca,'XTick');
    set(gca, 'XTick', 0:0.5:2.5)
    if i == 1
        xlabel('water height (m)','Fontsize',16)
        ylabel('Probability','Fontsize',16)
    end
        if i == 3
        
        lgnd = legend('SLR', 'GESLA-2', 'GTSR', ['GTSR+' char(10) 'waves'],'Location','northwest','Box','off')
        set(lgnd,'color','none');
        lgnd.FontSize = 12;
        end

end

