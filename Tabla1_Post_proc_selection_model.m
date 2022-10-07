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

for choice = 0:2 % 0 = full random / 1 = informed disp - random search / 2 = full informed
    
    Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS','AtoBe'};
    %% Mu fixed in colonies sampled
    if (choice == 0)
        load('post_proc_EP_7pm_random')
        Dbar = 2*sum(wm.*lL);
        mu_mean = Mus_mean;
        pm_mean = sum(wm.*Pm);
        D_mean  = sum(wm.*D);
        Dthetabar = 2*LogL_real_7pm(mu_mean,D_mean,pm_mean,0,0);

    elseif (choice == 1)
        load('post_proc_EP_7pm_informed_jd0')
        Dbar = 2*sum(wm.*lL);
        mu_mean = Mus_mean;
        pm_mean = sum(wm.*Pm);
        D_mean  = sum(wm.*D);
        Dthetabar = 2*LogL_real_7pm_informed(mu_mean,D_mean,pm_mean,0,0);

    elseif (choice == 2)
        load('post_proc_EP_7pm_informed_jd1')
        Dbar = 2*sum(wm.*lL);
        mu_mean = Mus_mean;
        pm_mean = sum(wm.*Pm);
        D_mean  = sum(wm.*D);
        Dthetabar = 2*LogL_real_7pm_informed(mu_mean,D_mean,pm_mean,1,0);
        
    end
    DIC1(choice+1) = 2*Dbar - Dthetabar;
    DIC2(choice+1) = Dbar +0.5*sum(wm.*(lL-Dbar*0.5).^2);
    IC(choice+1) = 3*Dbar - 2*Dthetabar;
    AIC(choice+1) = 2*lL_min; %2*min(lL);
end
SEL = [AIC;DIC1;DIC2;IC];
DAIC = (AIC-min(AIC));
score = exp(-DAIC)./sum(exp(-DAIC));

labelset = {'R', 'SI' , 'I'};% {'','',''}
figure(1)
clf
hold on
for i = 1:4
plot(SEL(i,:),'-o','MarkerSize',7,'MarkerFaceColor',Color(i,:),'color',Color(i,:))
end
legend('AIC','DIC1','DIC2','IC','Location','southeast')
xticks([1 2 3])
xticklabels(labelset)
xlabel('Dispersal behavior','Interpreter','latex','FontSize',16)


% figure(2)
% clf
% psi = PM{2};
% pi  = PM{3};
% for i = 1:7
%     subplot(4,3,i)
%     hold on
%     histogram(psi(:,i),0:0.05:1,'normalization','probability')
%     histogram(pi(:,i),0:0.05:1,'normalization','probability')
%     xlabel(Col{i})
% %     plot(p(:,i),'o')
% %     pause
% end
% %% Fig Marginal pdf of and Migration rate (Pm)
% if (choice == 0)
%     %%%% Random dispersal
%     figure(choice+1)
%     clf
%     hold on
%     histogram(PM{1})
% %     boxplot(Pm_post,'Labels', Col)
%     %%% Remplacer boxplot par violin plot
%     % violin(Pm_post,'Labels', Col)
%     % distributionPlot(Pm_post,'Labels', Col)
%     %%%
%     %boxplot(Pm,'Labels', Col)
%     hold off
%     % xlim([0,1])
%     xlabel('Mean migration rate')
% else
%     %%%% Informed search
%     figure(choice+1)
%     clf
%     hold on
%     boxplot(Em_post,'Labels', Col,'PlotStyle','compact')
%     %%% Remplacer boxplot par violin plot
%     % violin(Em_post(:,[3,4,5,7]),'Labels', Col)
%     % distributionPlot(Pm_post,'Labels', Col)
%     %%%
%     %boxplot(Pm,'Labels', Col)
%     hold off
%     % xlim([0,1])
%     xlabel('Mean migration rate')
% end

