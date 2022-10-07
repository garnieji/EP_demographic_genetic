clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ir = [7,3,6];     % Number of iteration process
Ir_rand = 2; % Number of iteration process parallel
ir = 0;
choice = 1; % 0 = full random / 1 = informed disp - random search / 2 = full informed

Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS','AtoBe'};

D  = [];  % Mean distance dispersal
Pm = [];  % Migration proportion
Mu = [];  % Initial proportion of subgroup
lL = [];  % Log Likelyhood
%% Random dispersal
D_r  = [];  % Mean distance dispersal
Pm_r = [];  % Migration proportion
Mu_r = [];  % Initial proportion of subgroup
lL_r = [];  % Log Likelyhood
%% Mu fixed in colonies sampled
for (ir= 0:Ir(choice+1))
    if (choice == 0)
        %     load(['par_main_EP_gene_7pm_pop_real_jd0_L_',num2str(ir),'_rand0.mat'])
        load(['par_main_EP_gene_7pm_pop_real_new_jd0_L_',num2str(ir),'_rand0.mat'])
        %     load(['par_main_EP_gene_7pm_pop_real_old_jd0_L_',num2str(ir),'_rand0.mat'])
    elseif (choice == 1)
        load(['par_main_EP_gene_7pm_pop_real_informed_jd0_L_',num2str(ir),'_rand0.mat'])
    elseif (choice == 2)
        load(['par_main_EP_gene_7pm_pop_real_informed_jd1_L_',num2str(ir),'_rand0.mat'])
    end
    D_r = [D_r;DDD];
    Pm_r = [Pm_r;PPm];
    Mu_r = [Mu_r;MU];
    lL_r = [lL_r;llL];
end
[nDr,mDr]=size(D_r);
D  = D_r(:);
Pm = zeros(length(D),7);
for i=1:7
    Pmi = Pm_r(i:7:end,:)';
    Pm(:,i) = Pmi(:);
end
Mu = Mu_r;
I_Mu = (1:nDr)'*ones(1,mDr);
lL = lL_r(:);

%% Minimizer of log(L)
[lL_min,Imin] = min(lL);
D_est  = D(Imin);
Pm_est = Pm(Imin,:);
Mu_est = Mu(I_Mu(Imin),:);

%% Important sampling
wm = exp(-lL)./(sum(exp(-lL)));
% [Dsort,Isort] = sort(D);
% wsort = wm(Isort);
% [Dsort,Isort] = unique(Dsort);
% wsort = wsort(Isort);
% dD = Dsort(2:end)-Dsort(1:end-1);
% dD = [dD;mean(dD)];
% wsort = wsort./sum(wsort.*dD);
% dd = 0.1;
% DD = 300:dd:6400;
% Wm = interp1(Dsort,wsort,DD);
% I_post = randsample(length(Wm),1e3,true,Wm);
% D_post = DD(I_post);

n_post = length(D);
N_data = length(D);
I_post = randsample(N_data,n_post,true,wm);
% [ll,I_post] = mink(lL,10);

D_post = D(I_post);
Pm_post = Pm(I_post,:);
Mu_post = Mu(I_Mu(I_post),:);


%%%% Emigration rate with informed dispersal
% if (choice ~= 0)
%     load('Rmedian.mat')
%     r = Rmedian(:,:,1);
%     i_subgroup = [1,5,8,21,25,39,46,55];
%     n_subgroup = length(i_subgroup)-1;
%     nsubpop = 4;
%     rM = 0.25;      % maximal death rate
%     n_post = 1e3;
%     i_post = randsample(N_data,n_post,true,wm);
%     d_post = D(i_post);
%     pm_post = Pm(i_post,:);
%     mu_post = Mu(I_Mu(i_post),:);
% 
%     Em_post = []; EE_post = [];
%     Tf = 2014;
%     It = Tf-2009;
%     parfor i = 1:n_post
%         pm = pm_post(i,:);
%         Em = [];EEm = [];
%         for s = 1:n_subgroup
%             Is = i_subgroup(s):(i_subgroup(s+1)-1);
%             rm = rM*(1-pm(s))*(pm(s)<1);
%             rs = r(Is,1:It);
%             em = ((rs<-rm) - (rs./(rm+(pm(s)==1))).*(-rm<=rs).*(rs<0))*(pm(s)<1) + (rs<0)*(pm(s)==1);
%             EEm{s} = em(:);
% %             EEm = [EEm,em(:)];
%             Em = [Em,median(em,'all')];
%         end
%         Em_post = [Em_post;Em];
%         EE_post = [EE_post;EEm];
%     end
% end

%% Figure Mu genetic repartition
Mus_mean = [];
figure(1)
clf
nMu = [];nnMu = []; W1MU_equ = [];W1MU_AMPG = []; W1MU_MAWS = []; W1MU_WEDD = [];
d_MU_pi = [];
pi_Clusters = [eye(4);1/4*ones(1,4)];
%%% Post proc
for i = 1:7
    Mus_mean = [Mus_mean,mean(Mu_post(:,4*i-3:4*i))'];
    for j = 1:5
    d_MU_pi(:,j) = sum((Mu_post(:,4*i-3:4*i) - pi_Clusters(j,:)).^2,2);
    end
     
%     
%     W1MU_equ = [W1MU_equ, sum(abs(cumsum(Mu_post(:,4*i-3:4*i),2)-(1:4)/4),2) ];
% 
%     W1MU_AMPG = [W1MU_AMPG, sum(abs(cumsum(Mu_post(:,4*i-3:4*i)-pi_AMPG,2)),2) ];
%     W1MU_MAWS = [W1MU_MAWS, sum(abs(cumsum(Mu_post(:,4*i-3:4*i)-pi_MAWS,2)),2) ];
%     W1MU_WEDD = [W1MU_WEDD, sum(abs(cumsum(Mu_post(:,4*i-3:4*i)-pi_WEDD,2)),2) ];
% %     nnMu = [nnMu,max((Mu_post(:,4*i-3:4*i) - 1/4),[],2)];
%     nnMu = [nnMu,(Mu_post(:,4*i-3:4*i)>1/4)];
    subplot(2,4,i)
%     hold on
    pie(Mus_mean(:,i))
    legend('AMPG','MAWS','ROSS','WEDD')
%     boxplot(Mu_post(:,4*i-3:4*i),'Labels',{'AMPG','MAWS','ROSS','WEDD'},'medianstyle','target','BoxStyle','filled')
%     boxplot(d_MU_pi,'Labels',{'AMPG','MAWS','ROSS','WEDD','EQU'},'medianstyle','target')
    hold off
%     xlabel(Col(i))
%     axis([0,5,0,1])
end

subplot(2,4,i+1)
% bar(categorical(Col,Col),Mus_mean','stacked')
pie(Mus_mean(:,1))
xlabel('Mean repartition of sub--population')
ylim([0 1])

figure
clf
bar(categorical(Col,Col),Mus_mean','stacked')
legend('AMPG','MAWS','ROSS','WEDD')
xlabel('Mean repartition of cluster')
ylim([0 1])

%% Fig Marginal pdf of Dispersdal (D)
figure(2)
clf
hold on
histogram(D_post,250:50:6500,'Normalization','probability')
D_mean = median(D_post);
line([D_mean,D_mean],[0,1],'color','r')
xlim([0,6500])
xlabel('Mean distance dispersal')

%% Fig Marginal pdf of and Migration rate (Pm)
if (choice == 0)
    %%%% Random dispersal
    figure(32)
    clf
    hold on
    boxplot(Pm_post,'Labels', Col)
    %%% Remplacer boxplot par violin plot
    % violin(Pm_post,'Labels', Col)
    % distributionPlot(Pm_post,'Labels', Col)
    %%%
    %boxplot(Pm,'Labels', Col)
    hold off
    % xlim([0,1])
    xlabel('Mean migration rate')
    
    figure(33)
    clf
    hold on
    for i = 1:7
        subplot(3,3,i)
%         histogram2(D,Pm(:,i),30,'DisplayStyle','tile','Normalization','probability')
        histogram(Pm_post(:,i),30,'Normalization','probability')
        xlabel(Col{i})
        drawnow
    end
    figure(34)
    clf
    hold on
    for i = 1:7
        subplot(3,3,i)
        histogram2(D_post,Pm_post(:,i),30,'DisplayStyle','tile','Normalization','probability')
        xlabel(Col{i})
        drawnow
    end
    
    figure
    clf
    hold on
    histogram(Pm_post(:,2),0:0.1:1,'Normalization','probability')

else
    %%%% Informed search
    figure(32)
    clf
    hold on
    boxplot(Em_post,'Labels', Col)
    %%% Remplacer boxplot par violin plot
    % violin(Em_post(:,[3,4,5,7]),'Labels', Col)
    % distributionPlot(Pm_post,'Labels', Col)
    %%%
    %boxplot(Pm,'Labels', Col)
    hold off
    % xlim([0,1])
    xlabel('Mean migration rate')
end

%% Compute the connectivity between colonies
%% Random strategy
d_r  = median(D_post);
if (choice ==0)
    pm_r = median(Pm_post,1);
else
    pm_r = median(Em_post,1);
end

i_subgroup = [1,5,8,21,25,39,46,55];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-20 Dawson to Kloa Point
% 21-24 Fold Island to Cape Darnley
% 25-38 Amanda Bay Point Geologie Davis Bay
% 39-45 Ross Sea
% 46-54 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;
lat= xlsread('COL_EP.xlsx','C2:C55'); % count for each colony in 2009 from Fretwell et al.
long =xlsread('COL_EP.xlsx','D2:D55');
ncol=length(lat);

d    = zeros(1,ncol);
lat  = lat*pi/180;
long = long*pi/180;
for c = 1:ncol-1
    d(c) = 6374.8925* acos(...
        sin(lat(c))*sin(lat(c+1))+...
        cos(lat(c))*cos(lat(c+1))*cos(long(c+1)-long(c)));
end
d(end) = 6374.8925* acos(sin(lat(end))*sin(lat(1))+...
    cos(lat(end))*cos(lat(1))*cos(long(1)-long(end)));
dm = dispersion(d,ncol);


Tf = 2014;
Pm_r = [];
for s = 1:n_subgroup
    Is = i_subgroup(s+1)-i_subgroup(s);
    Pm_r = [Pm_r;ones(Is,1)*pm_r(s)];
end
if (choice == 0)
    [a_subgroup,a_col,N] = EP_project_random_7col_migration(Tf,2009,d_r,Pm_r);
elseif (choice==1)
    jd = 0;
    [a_subgroup,a_col,N] = EP_project_informed_7col_migration(Tf,2009,d_r,Pm_r,jd);
elseif ( choice ==2)
    jd = 1;
    [a_subgroup,a_col,N] = EP_project_informed_7col_migration(Tf,2009,d_r,Pm_r,jd);
end


load('post_proc_EP_7pm_informed_jd0_connectivity')
choice = 1;
na = length(D_post);
A_col  = {}; %zeros(ncol,ncol,na);
A_subgroup = {}; %zeros(7,7,na);

parfor ik = 1:na
    Pm_r = [];
    pm_r = Pm_post(ik,:);
    for s = 1:n_subgroup
        Is = i_subgroup(s+1)-i_subgroup(s);
        Pm_r = [Pm_r;ones(Is,1)*pm_r(s)];
    end
    if (choice == 0)
        [a_subgroup,a_col,N] = EP_project_random_7col_migration(Tf,2009,D_post(ik),Pm_r);
    elseif (choice==1)
        jd = 0;
        [a_subgroup,a_col,N] = EP_project_informed_7col_migration(Tf,2009,D_post(ik),Pm_r,jd);
    elseif ( choice ==2)
        jd = 1;
        [a_subgroup,a_col,N] = EP_project_informed_7col_migration(Tf,2009,D_post(ik),Pm_r,jd);
    end
%     a_subgroup_med = median(a_subgroup,3);
%     a_col_med      = median(a_col,3);
%     A_col(:,:,ik)  = a_col_med;
%     A_subgroup(:,:,ik) = a_subgroup_med;
    A_subgroup{ik} = a_subgroup;
    A_col{ik}      = a_col;
end
A_col_med = median(A_col,3);
A_subgroup_med = median(A_subgroup,3);


%% Figure of mean connectivity between colonies
% tt = find(time == Tt)-1;
na =length(D_post);
a_Col = zeros(ncol,ncol,5*na);
for ik = 1:na
    a_Col(:,:,5*ik-4:5*ik) = A_col{ik};
end
A_Col_med = median(a_Col,3);


[s,t] = meshgrid((1:ncol),(1:ncol));
Ga_col_med = digraph(100*A_Col_med');
LWidths_med = 5*Ga_col_med.Edges.Weight/max(Ga_col_med.Edges.Weight);

% Plot of the region of Antarctica
R = 6371;
latS = -pi/2;
longS = 0;
k = 2*R./(1+sin(latS)*sin(lat)+cos(latS)*cos(lat).*cos(long-longS));
x = k.*cos(lat).*sin(long-longS);
y = k.*(cos(latS)*sin(lat) - sin(latS)*cos(lat).*cos(long-longS)) ;
Xi = {};
Yi = {};
X = [];
Y = [];
for i = 1:7
    xi = x(i_subgroup(i):i_subgroup(i+1)-1);
    yi = y(i_subgroup(i):i_subgroup(i+1)-1);
    Xi = [Xi,xi];
    Yi = [Yi,yi];
    X = [X,median(xi)];
    Y = [Y,median(yi)];
end

figure
clf
hold on
plot([x;x(1)],[y;y(1)],'k-')
for i=1:7
    plot(Xi{1,i},Yi{1,i},'-o','linewidth',3)
end
plot(Ga_col_med,'XData',x,'YData',y,'EdgeLabel',Ga_col_med.Edges.Weight,'LineWidth',LWidths_med,'ArrowSize',15)
% plot(Ga_mean,'XData',x,'YData',y,'EdgeLabel',Ga_mean.Edges.Weight,'LineWidth',LWidths_mean,'ArrowSize',15)

%% Figure of the mean connectivity between regions
Ga_subgroup = digraph(100*A_subgroup_med);%.*(AA>0.0005));
LaaWidths = 5*Ga_subgroup.Edges.Weight/max(Ga_subgroup.Edges.Weight);

figure
clf
hold on
plot([x;x(1)],[y;y(1)],'k-')
for i=1:7
    plot(Xi{1,i},Yi{1,i},'-o','linewidth',3)
end
plot(Ga_subgroup,'XData',X,'YData',Y,'NodeLabel',Col,'EdgeLabel',Ga_subgroup.Edges.Weight,'LineWidth',LaaWidths,'ArrowSize',15)
for i =1:7
    plot(X(i),Y(i),'o')
end
drawnow


% 
% if (choice == 0)
%     save('post_proc_EP_7pm_random');
% else
%     save(['post_proc_EP_7pm_informed_jd',num2str(jd)]);
% end















% 
% A  = zeros(ncol,ncol,length(D_post));
% AA = zeros(7,7);
% Tt = 2010;
% for ik = 1:length(D_post)
%     Dd = D_post(ik);
%     if (choice == 0)
%         pm = Pm_post(ik,:);
%     else
%         pm = Em_post(ik,:);
%     end
%     ppm = [];
%     for s = 1:n_subgroup
%         Is = i_subgroup(s+1)-i_subgroup(s);
%         ppm = [ppm;ones(Is,1)*pm(s)];
%     end
%     
%     % Dispersal kernel       % Rem: k(0)=0, the probability of staying where we are, is zeroa
%     k = @(x) (x<Dd)/Dd .*(x>0);  % Uniform kernel
%     
%     % Matrix of connexion without renormalisation
%     %     dm       = dispersion(d,ncol); % Distance between colonies
%     dmat = k(dm);    % Unif kernel
%     
%     % Renormalization of the connexion matrix %
%     den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
%     Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
%     Dmat = dmat*diag(1./Den);
%     
%     
%     disp_sel = Dmat;
%     den_disp = sum(disp_sel,1)';
%     Disp_sel = disp_sel-diag(den_disp); % pour assurer que la somme des colonnes vaut 0
%     Disp = Disp_sel - diag(diag(Disp_sel));
%     a = (Disp*diag(ppm));
%     A(:,:,ik) = a;
% end



% %% Figure of the mean connectivity between regions
% for i = 1:7
%     for j=1:7
%         %         Aij = A( i_subgroup(j):i_subgroup(j+1)-1,i_subgroup(i):i_subgroup(i+1)-1,:);
%         Aij = A_med(i_subgroup(j):i_subgroup(j+1)-1,i_subgroup(i):i_subgroup(i+1)-1);
%         Aij = Aij(:);
%         if(sum(Aij)>0)
%             aij = median(Aij(Aij>0));
%             %             aij = median(Aij);
%         else
%             aij=0;
%         end
%         AA(i,j) = aij;
%     end
% end
% 
% % Gaa = digraph(100*AA);
% Gaa = digraph(100*AA);%.*(AA>0.0005));
% LaaWidths = 5*Gaa.Edges.Weight/max(Gaa.Edges.Weight);
% 
% figure
% clf
% hold on
% plot([x;x(1)],[y;y(1)],'k-')
% for i=1:7
%     plot(Xi{1,i},Yi{1,i},'-o','linewidth',3)
% end
% plot(Gaa,'XData',X,'YData',Y,'NodeLabel',Col,'EdgeLabel',Gaa.Edges.Weight,'LineWidth',LaaWidths,'ArrowSize',15)
% for i =1:7
%     plot(X(i),Y(i),'o')
% end
% drawnow

%% Fig boxplot migration rates
% figure
% clf
% boxplot(aa_pop)
%
% figure
% clf
% Aa_pop = permute(AA_pop,[3,1,2]);
% Aa_pop = Aa_pop(:,:);
% boxplot(Aa_pop)
% hold on
% plot(median(Aa_pop))
% plot(quantile(Aa_pop,0.01))
% plot(quantile(Aa_pop,0.99))

