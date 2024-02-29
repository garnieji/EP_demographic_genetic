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
[nDr,mDr]=size(D_r);
I_Mu = (1:nDr)'*ones(1,mDr);
Muu = Mu(I_Mu,:);
%% Likelihood with respect to the differnet parameters
%% Dispersal d
lL_min = lL;
I_min = 1:length(lL);
D_min = D(I_min);

fl = @(x) 1./(1+(x-min(lL))/20);
fl_min = fl(lL_min);

figure(2)
clf
plot(D_min,fl_min,'+')
xlim([250,6500])
xlabel('Mean distance dispersal ($d$) in km ','interpreter','latex','FontSize',16)
ylabel('$f(\Theta)$','interpreter','latex','FontSize',16)

%% Mu
[lL_min,I_min] = mink(lL,1e2);
fl = @(x) 1./(1+(x-min(lL))/100);

fl_min = fl(lL_min);

Mu_min = Muu(I_min,:);

figure(20)
clf
scatter3(Mu_min(:,1),Mu_min(:,2),Mu_min(:,3),20*fl_min,fl_min,'filled')
hold on
plot3(Mu_min(1,1),Mu_min(1,2),Mu_min(1,3),'o')


%% Pm
lL_min = lL;
I_min = 1:length(lL);
Pm_min = Pm(I_min,:);

fl = @(x) 1./(1+(x-min(lL))/20);
fl_min = fl(lL_min);
for i = 1:7
figure(30+i)
clf
hold on
plot(Pm_min(:,i),fl_min,'+')
xlabel(['emigration parameter ($p_{m}$) in ',Col(i)],'interpreter','latex','FontSize',16)
ylabel('$f(\Theta)$','interpreter','latex','FontSize',16)
end






