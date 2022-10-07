clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];

Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS',{'A-B', 'seas'}};

%% Growth rate per colony / Carrying capacity and initial condition
% load('Growth_final_nsim20_Jan25.mat')
SCENARIO = 1:5; % Scenario climatic: 1= RCP8.5 | 2= 2.6°C | 3 = 2.4°C RCP4.5 | 4= Paris 2 | 5= Paris 1.5
ns = length(SCENARIO);
scenEXT = 1;  % Scenario extreme event: 1=none

% Choice = 1 SI
% load('post_proc_EP_7pm_informed_jd0','D_post','Pm_post','D_mean','Pm_est','D_est')
ncol = 54;
%% Sub pop
i_subgroup = [1,5,8,21,25,39,46,55];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-20 Dawson to Kloa Point
% 21-24 Fold Island to Cape Darnley
% 25-38 Amanda Bay Point Geologie Davis Bay
% 39-45 Ross Sea
% 46-54 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;
% Pm   = zeros(size(Pm_post,1),ncol);
% Pm_e = zeros(1,ncol);
% for j = 1:n_subgroup
%     Pm(:,i_subgroup(j):i_subgroup(j+1)-1) = Pm_post(:,j)*ones(1,i_subgroup(j+1)-i_subgroup(j));
%     Pm_e(1,i_subgroup(j):i_subgroup(j+1)-1) = Pm_est(j)*ones(1,i_subgroup(j+1)-i_subgroup(j));
% end
% Pm_mean = median(Pm,1);


Tf = 2100;
T0 = 2009;
tt = T0:Tf;
nt = length(tt);

% nsimu = length(D_post);
% nsimu = nsimulation;

% Ntot_disp   = zeros(nsimu,nt,ns);
% Ntot_nodisp = zeros(nsimu,nt,ns);
% is = 1;
% for scenario = SCENARIO
%     rmedian = IC(scenario).EXT(scenEXT).Rmedian(:,:);
%     R       = permute(IC(scenario).EXT(scenEXT).R, [3 2 1]); % IC(scenario).R: sim*nt*col
%     NN_disp   = zeros(nsimu,nt);
%     NN_nodisp = zeros(nsimu,nt);
%     parfor i = 1:nsimu
%         ir = mod(i,nsimulation)+nsimulation.*(mod(i,nsimulation)==0);
%         r = R(:,:,ir);
%         pm = Pm(i,:)';
% %         [N_nodisp,N_disp] = EP_project_informed(Tf,T0,D_post(i),pm,0,r,rmedian);
% %         [N_nodisp,N_disp] = EP_project_informed(Tf,T0,D_mean,Pm_mean',0,r,rmedian);
%         [N_nodisp,N_disp] = EP_project_informed(Tf,T0,D_est,Pm_e',0,r,rmedian);
% %         [N_nodisp,N_disp] = EP_project_informed(Tf,T0,D_post(i),pm,0,r,r);
%         NN_disp(i,:) = sum(N_disp);
%         NN_nodisp(i,:) = sum(N_nodisp);
%     end
%     
%     Ntot_nodisp(:,:,is) = NN_nodisp;
%     Ntot_disp(:,:,is)   = NN_disp;
%     is = is+1;
% end
% Ntot_nodisp_med = quantile(Ntot_nodisp,[0.05,0.5,0.95],1);
% Ntot_disp_med = quantile(Ntot_disp,[0.05,0.5,0.95],1);
% save('post_proc_projection_SI_mean.mat','-v7.3','Ntot_disp','Ntot_nodisp','Ntot_disp_med','Ntot_nodisp_med')

%% Figure
% load('post_proc_projection_SI.mat','Ntot_disp','Ntot_disp_med')
% load('post_proc_projection_SI_r.mat')
load('post_proc_projection_SI_mean.mat') %,'Ntot_disp','Ntot_disp_med')

% Ntot_global = readmatrix('Jimmy_Nglobal_noext_alldispersion.xlsx'); % count for each colony in 2009 from Fretwell et al. 
% Ntot_global = Ntot_global(:,2:end);
% Ntot_global = Ntot_global(2:end,:);

load('metapop_finalsim20_Jan25_Nglobal.mat')
Ntot_global = permute(Ntot_global_med,[3,2,1]);
Diff = zeros(2000,nt,5);
for ic = SCENARIO
    Diff(:,:,ic) = (Ntot_disp(:,:,ic)-Ntot_global(ic,:))./Ntot_global(ic,:)*100;
end
Diff_med = quantile(Diff,[0.05,0.5,0.95],1);

Diff_disp = (Ntot_disp-Ntot_nodisp)./Ntot_nodisp*100;
Diff_disp_med = quantile(Diff_disp,[0.05,0.5,0.95],1);

Diff_climate = (Ntot_global(2:end,:)-Ntot_global(1,:))./Ntot_global(1,:)*100;
% Diff_climate_med = quantile(Diff_climate,[0.05,0.5,0.95],1);
% Color = colormap(hot(5));
Color = [128,0,0;...
         210,0,0;...
         210,85,0;...
         255,153,85;...
         255,214,180]./255;
figure(1)
clf
hold on
ttt = [tt, fliplr(tt)];
for ic = SCENARIO
    % ic = 1;
    inBetween = [Ntot_disp_med(1,:,ic), fliplr(Ntot_disp_med(3,:,ic))];
    f = fill(ttt, inBetween,[0.8,0.8,0.8]);
    set(f,'EdgeColor','none','FaceAlpha', 0.4)
    plot(tt,Ntot_global(ic,:),'--','color',[0.5,0.5,0.5],'linewidth',2)
    plot(tt,Ntot_disp_med(2,:,ic),'-','color',Color(ic,:),'linewidth',1.5)
    
    % ic = 2;
    % inBetween = [Ntot_nodisp_med(1,:,ic), fliplr(Ntot_nodisp_med(3,:,ic))];
    % f = fill(ttt, inBetween,Color(3,:));
    % set(f,'EdgeColor','none','FaceAlpha', 0.2)
%     plot(tt,Ntot_nodisp_med(2,:,ic),'-.','color',[0.5,0.5,0.5],'linewidth',2)
    % plot(tt,Ntot_disp_med(2,:,ic),'-','color',Color(3,:),'linewidth',1.5)
    
end

ylim([0,3e5])
xlim([T0,Tf])
xlabel('Time (years)','Interpreter','latex','FontSize',15)
ylabel({'Total population size','over Antarctic continent'},'fontsize',15,'interpreter','latex');

%%% Diff

figure(2)
clf
hold on
line([T0,Tf],[0,0],'color','k','linewidth',2)

for ic = SCENARIO
    
    inBetween = [Diff_med(1,:,ic), fliplr(Diff_med(3,:,ic))];
    f = fill(ttt, inBetween,[0.8,0.8,0.8]);
    set(f,'EdgeColor','none','FaceAlpha', 0.2)
    
    plot(tt,Diff_med(2,:,ic),'-','color',Color(ic,:),'linewidth',3)
end
axis([T0,Tf,-15,30])
xlabel('Time (years)','Interpreter','latex','FontSize',15)
ylabel({'Percentage difference in', 'total population size' },'fontsize',15,'interpreter','latex');

figure(21)
clf
hold on
line([T0,Tf],[0,0],'color','k','linewidth',2)

for ic = SCENARIO
    
    inBetween = [Diff_disp_med(1,:,ic), fliplr(Diff_disp_med(3,:,ic))];
    f = fill(ttt, inBetween,[0.8,0.8,0.8]);
    set(f,'EdgeColor','none','FaceAlpha', 0.2)
    
    plot(tt,Diff_disp_med(2,:,ic),'-','color',Color(ic,:),'linewidth',3)
end
axis([T0,Tf,-10,10])
xlabel('Time (years)','Interpreter','latex','FontSize',15)
ylabel({'Percentage difference in', 'total population size' },'fontsize',15,'interpreter','latex');


figure(3)
clf
hold on
line([T0,Tf],[0,0],'color','k','linewidth',2)

for ic = 1:4
%     
%     inBetween = [Diff_climate_med(1,:,ic), fliplr(Diff_climate_med(3,:,ic))];
%     f = fill(ttt, inBetween,[0.8,0.8,0.8]);
%     set(f,'EdgeColor','none','FaceAlpha', 0.2)
    
    plot(tt,Diff_climate(ic,:),'-','color',Color(ic+1,:),'linewidth',3)
end
xlim([T0,Tf])
xlabel('Time (years)','Interpreter','latex','FontSize',15)
ylabel({'Percentage difference in', 'total population size' },'fontsize',15,'interpreter','latex');











