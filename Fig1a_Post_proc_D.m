clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];
Ir = [7,3,2];     % Number of iteration process
Ir_rand = 2; % Number of iteration process parallel
ir = 0;
D_mean = [];
D_POST = [];
Reg_D  = [];
D_pos = {'Random','Semi-Informed','Informed'};
y_d = [0.002,0.01,0.02];
for choice = 1 %0:2 % 0 = full random / 1 = informed disp - random search / 2 = full informed
    
    Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS','AtoBe'};
    if (choice == 0)
        load('post_proc_EP_7pm_random.mat')
    elseif(choice == 1)
        load('post_proc_EP_7pm_informed_jd0.mat')
    else
        load('post_proc_EP_7pm_informed_jd1.mat')        
    end
  
    d_med  = median(D_post);
    D_mean = [D_mean,d_med];  
    
    %% Fig Marginal pdf of Dispersdal (D)
    figure(1)
    hold on
    histogram(D_post,250:50:6500,'Normalization','pdf')
    xlim([0,600]) %xlim([0,3000])
    xlabel('Mean distance dispersal ($d$) in km ','interpreter','latex','FontSize',16)
    ylabel('Empirical distribution','interpreter','latex','FontSize',16)
    line([d_med,d_med],[0,y_d(choice+1)],'color',Color(1,:))  
    
end

% figure(2)
% clf
% vp_D = violinplot(D_POST,Reg_D,'GroupOrder',D_pos);
% ylabel('Mean distance dispersal ($d$)','interpreter','latex','FontSize',16)
% ylim([0,3000])


% for i = 1%1:3
% line([D_mean(i),D_mean(i)],[0,0.01],'color',Color(i,:))
% end





