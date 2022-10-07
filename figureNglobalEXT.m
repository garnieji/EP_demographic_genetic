%% global and regional population for all demographic scenario combined:
% 4 extreme scenarios (remove the case without extreme) and
% all dispersive scneario
load('metapop_finalnsim20_Jan25.mat')

%%
time=2009:2100;
% combined extreme and dispersal demo scenario
nscenario = 5;
nscenEXT  = 1;%5;
nscenDISP = 9;
ic = 1;
for scenario=1:nscenario
    SCEN(scenario).Ncombined=zeros(nsimulation*nscenDISP*nscenEXT,nt,nc);
    for c=1:nc
        i1=1;i2=nsimulation;
        x=1;
        for scenEXT=1:nscenEXT
            for scenDISP=1:nscenDISP
                x=x+1;
                SCEN(scenario).Ncombined(i1:i2,:,c)=...
                    SC(scenario).EXT(scenEXT).DISP(scenDISP).col(c).Ncol;
                i1=i2+1;
                i2=x*nsimulation;
            end
            
        end
    end
end

% regional and global ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
nsimTOT=nsimulation*nscenDISP*(nscenEXT);
for scenario=1:nscenario
    SCEN(scenario).Nglobal=zeros(nsimTOT,nt);
    SCEN(scenario).NWED=zeros(nsimTOT,nt);
    SCEN(scenario).NIO=zeros(nsimTOT,nt);
    SCEN(scenario).NWPO=zeros(nsimTOT,nt);
    SCEN(scenario).NRS=zeros(nsimTOT,nt);
    SCEN(scenario).NBAS=zeros(nsimTOT,nt);
    
    
    for c=1:nc
        SCEN(scenario).Nglobal=SCEN(scenario).Nglobal+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    
    for c=1:13
        x=0;
        SCEN(scenario).NWED=SCEN(scenario).NWED+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    
    for c=14:27
        x=0;
        SCEN(scenario).NIO=SCEN(scenario).NIO+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    for c=28:38
        x=0;
        SCEN(scenario).NWPO=SCEN(scenario).NWPO+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    
    for c=39:45
        x=0;
        SCEN(scenario).NRS=SCEN(scenario).NRS+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    
    for c=46:54
        x=0;
        SCEN(scenario).NBAS=SCEN(scenario).NBAS+ ...
            SCEN(scenario).Ncombined(:,:,c);
    end
    
    SCEN(scenario).HIglobal=prctile(SCEN(scenario).Nglobal,90);
    SCEN(scenario).LOWglobal=prctile(SCEN(scenario).Nglobal,5);
    SCEN(scenario).MEDglobal=prctile(SCEN(scenario).Nglobal,50);
    
    SCEN(scenario).CEreg(:,:,1)=...
        [prctile(SCEN(scenario).NWED,95)' ...
        prctile(SCEN(scenario).NWED,50)'...
        prctile(SCEN(scenario).NWED,5)'];
    
    SCEN(scenario).CEreg(:,:,2)=...
        [prctile(SCEN(scenario).NIO,95)' ...
        prctile(SCEN(scenario).NIO,50)'...
        prctile(SCEN(scenario).NIO,5)'];
    
    SCEN(scenario).CEreg(:,:,3)=...
        [prctile(SCEN(scenario).NWPO,95)' ...
        prctile(SCEN(scenario).NWPO,50)'...
        prctile(SCEN(scenario).NWPO,5)'];
    
    SCEN(scenario).CEreg(:,:,4)=...
        [prctile(SCEN(scenario).NRS,95)' ...
        prctile(SCEN(scenario).NRS,50)'...
        prctile(SCEN(scenario).NRS,5)'];
    
    SCEN(scenario).CEreg(:,:,5)=...
        [prctile(SCEN(scenario).NBAS,95)' ...
        prctile(SCEN(scenario).NBAS,50)'...
        prctile(SCEN(scenario).NBAS,5)'];
end

%% each extrem scenario

for scenario=1:nscenario
    for scenEXT=1:nscenEXT %scenEXT=1:nscenEXT
        SCEN(scenario).EXT(scenEXT).Ncombineddisp=...
            zeros(nsimulation*nscenDISP,nt,nc);% with no extreme event
        for c=1:nc
            i1=1;i2=nsimulation;
            x=1;
            for scenDISP=1:nscenDISP
                x=x+1;
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(i1:i2,:,c)=...
                    SC(scenario).EXT(scenEXT).DISP(scenDISP).col(c).Ncol;
                i1=i2+1;
                i2=x*nsimulation;
            end
            
        end
    end
end

% regional and global ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
nsimTOT=nsimulation*nscenDISP;

for scenario=1:nscenario
    for scenEXT=1:nscenEXT
        SCEN(scenario).EXT(scenEXT).Nglobal=zeros(nsimTOT,nt);
        SCEN(scenario).EXT(scenEXT).NWED=zeros(nsimTOT,nt);
        SCEN(scenario).EXT(scenEXT).NIO=zeros(nsimTOT,nt);
        SCEN(scenario).EXT(scenEXT).NWPO=zeros(nsimTOT,nt);
        SCEN(scenario).EXT(scenEXT).NRS=zeros(nsimTOT,nt);
        SCEN(scenario).EXT(scenEXT).NBAS=zeros(nsimTOT,nt);
        
        
        for c=1:nc
            SCEN(scenario).EXT(scenEXT).Nglobal=...
                SCEN(scenario).EXT(scenEXT).Nglobal+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
        
        
        
        for c=1:13
            x=0;
            SCEN(scenario).EXT(scenEXT).NWED=...
                SCEN(scenario).EXT(scenEXT).NWED+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
        
        for c=14:27
            x=0;
            SCEN(scenario).EXT(scenEXT).NIO=...
                SCEN(scenario).EXT(scenEXT).NIO+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
        for c=28:38
            x=0;
            SCEN(scenario).EXT(scenEXT).NWPO=...
                SCEN(scenario).EXT(scenEXT).NWPO+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
        
        for c=39:45
            x=0;
            SCEN(scenario).EXT(scenEXT).NRS=...
                SCEN(scenario).EXT(scenEXT).NRS+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
        
        for c=46:54
            x=0;
            SCEN(scenario).EXT(scenEXT).NBAS=...
                SCEN(scenario).EXT(scenEXT).NBAS+ ...
                SCEN(scenario).EXT(scenEXT).Ncombineddisp(:,:,c);
        end
    end
end
%%
%CT=cbrewer('seq', 'Greys', 5)
% CT = cbrewer('seq', 'Oranges', 5)
CT = [128,0,0;...
         210,0,0;...
         210,85,0;...
         255,153,85;...
         255,214,180]./255;
fontsi=36;Fo=fontsi;

%arraysc=[1 2 4 5]; titre={'RCP8.5','Special','Paris 2','Paris1.5'};
%arraysc=1:4; titre={'RCP8.5','Special','RCP4.5','Paris 2'};
arraysc=[1 2 3 5]; titre={'RCP8.5','Special','RCP4.5','Paris 1.5'};

figure
for sc=1%:4
%     subplot(2,2,sc) % GLOBAL

    scenario=arraysc(sc);
    HIglobal=prctile(SCEN(scenario).Nglobal,97.5);
    LOWglobal=prctile(SCEN(scenario).Nglobal,2.5);
    MEDglobal=prctile(SCEN(scenario).Nglobal,50);
    
    hold on
    
    area(time,HIglobal,'LineStyle','none','FaceColor',[0.6 0.6 0.6],'FaceAlpha',.3)
    area(time,LOWglobal,'LineStyle','none','FaceColor','w')
    plot(time,HIglobal,'color',[0.5 0.5 0.5],'LineWidth',2)
    plot(time,LOWglobal,'color',[0.5 0.5 0.5],'LineWidth',2)
    
    % Plot PopTot
    plot(time,MEDglobal,'-k', 'linewidth',5)
    plot(time,prctile(SCEN(scenario).EXT(1).Nglobal,50),'color',CT(1,:),'linewidth',3)

    for scenEXT=2:nscenEXT
        plot(time,prctile(SCEN(scenario).EXT(scenEXT).Nglobal,50),'color',CT(scenEXT,:),'linewidth',5)
    end
    ylabel('Total population size','fontsize',fontsi)%,'interpreter','latex');
    xlabel('Years','fontsize',fontsi)%,'interpreter','latex');
    %legend('Without Dispersion','RD / pm=0.1 ','RD / pm=0.9','INF / pm=0.1','INF / pm=0.9')
    title('Global Antarctic population','fontsize',fontsi)
    set(gca,'FontSize',fontsi, ...
        'FontName'   , 'Helvetica',...
        'FontWeight' , 'bold',...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'YTick',...
        [0 50000 100000 150000 200000 250000 300000],...
        'YTickLabel',...
        {'0','50000','100000','150000','200000','250000','300000'},...
        'Box'         , 'off'     , ...% Define box around whole figure
        'LineWidth'   , 3);
    xlim([2009 2100])
     title(titre(sc),'FontSize',48);
end
%%
     legend('','',' ', ' ', 'No extreme events',...
    'Historical extreme on reproduction',...
    'Historical extreme on reproduction and survival',...
    'Sea ice related extreme on reproduction',...
    'Sea ice related extreme on reproduction and survival',...
'FontSize',48,'FontWeight' , 'bold')
%%

%CT=cbrewer('seq', 'Greys', 5)
% CT=cbrewer('seq', 'Oranges', 5)
CT = [128,0,0;...
         210,0,0;...
         210,85,0;...
         255,153,85;...
         255,214,180]./255;
fontsi=36;Fo=fontsi;

%arraysc=[1 2 4 5]; titre={'RCP8.5','Special','Paris 2','Paris1.5'};
%arraysc=1:4; titre={'RCP8.5','Special','RCP4.5','Paris 2'};
arraysc=1:5; titre={'RCP8.5','Special','RCP4.5','Paris 2','Paris1.5'};

figure
for sc=1:5
    subplot(3,3,sc) % GLOBAL

    scenario=arraysc(sc);
    HIglobal=prctile(SCEN(scenario).Nglobal,97.5);
    LOWglobal=prctile(SCEN(scenario).Nglobal,2.5);
    MEDglobal=prctile(SCEN(scenario).Nglobal,50);
    
    hold on
    
    area(time,HIglobal,'LineStyle','none','FaceColor',[0.6 0.6 0.6],'FaceAlpha',.3)
    area(time,LOWglobal,'LineStyle','none','FaceColor','w')
    plot(time,HIglobal,'color',[0.5 0.5 0.5],'LineWidth',2)
    plot(time,LOWglobal,'color',[0.5 0.5 0.5],'LineWidth',2)
    
    % Plot PopTot
    plot(time,MEDglobal,'-k', 'linewidth',5)
    plot(time,prctile(SCEN(scenario).EXT(1).Nglobal,50),'color',CT(1,:),'linewidth',2)

    for scenEXT=2:nscenEXT
        plot(time,prctile(SCEN(scenario).EXT(scenEXT).Nglobal,50),'color',CT(scenEXT,:),'linewidth',5)
    end
    ylabel('Total population size','fontsize',fontsi)%,'interpreter','latex');
    xlabel('Years','fontsize',fontsi)%,'interpreter','latex');
    %legend('Without Dispersion','RD / pm=0.1 ','RD / pm=0.9','INF / pm=0.1','INF / pm=0.9')
    title('Global Antarctic population','fontsize',fontsi)
    set(gca,'FontSize',fontsi, ...
        'FontName'   , 'Helvetica',...
        'FontWeight' , 'bold',...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'YTick',...
        [0 50000 100000 150000 200000 250000 300000],...
        'YTickLabel',...
        {'0','50000','100000','150000','200000','250000','300000'},...
        'Box'         , 'off'     , ...% Define box around whole figure
        'LineWidth'   , 3);
    xlim([2009 2100])
     title(titre(sc),'FontSize',48);
end
%%
     legend('','',' ', ' ', 'No extreme events',...
    'Historical extreme on reproduction',...
    'Historical extreme on reproduction and survival',...
    'Sea ice related extreme on reproduction',...
    'Sea ice related extreme on reproduction and survival',...
'FontSize',48,'FontWeight' , 'bold')
