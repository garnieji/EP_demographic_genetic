clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = get(gca,'colororder');
Color(end+1,:) = [0.5,0.5,0.5];
Marker = ['o','*','d','^'];
Ir = [7,3,2];     % Number of iteration process
Ir_rand = 2; % Number of iteration process parallel
ir = 0;
Col = {'StoS ','WEDD ','StoKP','MAWS ','AMPG ','ROSS ','A-B seas'};

%% Choice = 1 SI
choice = 1;
if (choice==0)
    load('post_proc_EP_7pm_random')
    EM_post = Pm_post;  
    Em_post = EM_post;
else
    if (choice == 1)
        load('post_proc_EP_7pm_informed_jd0')
    else
        load('post_proc_EP_7pm_informed_jd1')
    end
    %% Input
    % %%%%%%%%%%%%%%%%%%%%%%
    Tf = 2014;                 % Final time
    T0 = 2009;                 % Initial time
    % mu = randfixedsum(n_subpop,1,1,0,1)';
    % Mu = ones(ncol,1)*mu;  % Initial repartition of the sub group inside each colony
    % D  = 1000;                 % Mean distance dispersal
    % pm = 0.1*ones(ncol,1);     % mean proportion of migrant in each colony
    % % over time
    % jd = 1;                    % Choose 1) jd=0 Random strategy; 2) jd=1 Oriented strategy
    % nsubpop = 4;               % Nb of subpop
    
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
    
    %% Initialization %
    %%%%%%%%%%%%%%%%%%
    ncol = 54; % Nb of colonies
    % nt   = 92; % Nb of time iteration from 2009 to 2100
    time = T0:Tf;
    %%%% calcul la distance cotiere entre chaque colonnie en km
    lat  = xlsread('COL_EP.xlsx','C2:C55'); % count for each colony in 2009 from Fretwell et al.
    long = xlsread('COL_EP.xlsx','D2:D55');
    ncol = length(lat);
    
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
    dm = dispersion(d,ncol);   % Distance between colonies
    
    %% Growth Function %
    %%%%%%%%%%%%%%%%%%%
    % With Density dependance K
    % g: Individual growth rate f(N)/N
    g = @(N,r,K) exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*(r<=0);
    
    %*********************************
    jd = 0; % Choose 1) jd=0 Random strategy; 2) jd=1 Oriented strategy
    
    %%%% Growth rate per colony / Carrying capacity and initial condition
    load('Rmedian.mat')
    r = Rmedian(:,:,1);
    
    BE = xlsread('COL_EP.xlsx','F2:F55'); % count for each colony in 2009 from Fretwell et al.
    K  = 2*BE; % carrying capacity
    
    
    %% Computation **********
    nt = length(T0:Tf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialization of size sub-population inside eachcolonies %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Em_post = [];
    Ncol = [];
    Ntot = [];
    parfor i =1:length(D_post)
        em_post = [];
        pm = [];
        for j = 1:n_subgroup
            Is = i_subgroup(j+1)-i_subgroup(j);
            pm = [pm;ones(Is,1)*Pm_post(i,j)];
        end
        %% Movement scenario %
        %%%%%%%%%%%%%%%%%%%%%
        rM = 0.25;      % maximal death rate
        rm = rM*(1-pm);
        if (pm==1)
            movement=@(r)  (r<0); %HIGH
        else
            movement=@(r)  (r<-rm) - (r./rm).*(-rm<=r).*(r<0); %HIGH
        end
        
        %% Matrix of connexion %
        %%%%%%%%%%%%%%%%%%%%%%%
        D = D_post(i); % Mean distance dispersal
        % D = 1000;               % Mean distance dispersal
        % Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zeroa
        k = @(x) (x<D)/D .*(x>0); % Uniform kernel
        
        % Matrix of connexion without renormalisation
        dmat = k(dm);    % Unif kernel
        % Renormalization of the connexion matrix %
        den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
        Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
        Dmat = dmat*diag(1./Den);
        
        %% Time loop forward
        % Inidividuals
        t = T0;
        it = 0;
        Nnew = BE;  % N(i)  = # individuals in each colony
        N = Nnew;
        while (t<Tf)
            Nold = Nnew;
            
            %%% Reproduction at time t+1
            rstar = g(Nold,r(:,it+1),K)-1; % Effective growth rate at t ime t
            
            if (choice==1)
                % Random Dispersal
                disp_sel = Dmat;
            elseif (choice==2)
                %%%% Construction of informed dispersal matrix
                mat_r = meshgrid(rstar)';
                [rmax_vec, ind_col_sel] = max(mat_r.*(dm<D)+((dm<D)-1)); % ((dm<D)-1) to eliminate zeros and the max provieds the largest negative growth rate
                Mat_select = zeros(ncol,ncol);
                for j=1:ncol
                    Mat_select(ind_col_sel(j),j) = 1;
                end
                rmax_mat = meshgrid(rmax_vec);
                % Oriented dispersal
                disp_sel = Mat_select;
            end
            
            den_disp = sum(disp_sel,1)';
            Disp_sel = disp_sel-diag(den_disp); % pour assurer que la somme des colonnes vaut 0
            
            
            %%% Exploration
            lambda = movement(rstar); % mean proportion of indivudu that lives colonies
            fn     = (1+rstar).*Nold;
            Nnew  = Disp_sel*(lambda.*fn) + fn;
            
            em_post = [em_post,lambda];
            
            N = [N,Nnew];
            t = t+1;
            it = it+1;
        end
        Ncol(i,:,:) = N;
        Ntot(i,:)   = sum(N);
        Em_post(:,i,:) = em_post;
    end
    
    EM_post = [];
    for i = 1:ncol
        emi = Em_post(i,:,:);
            EM_post(i,:) = emi(:);
%         EM_post(i,:) = median(emi,2);
    end
end
% Em_post = sum(A_col,1);
% Em_post = permute(Em_post,[2,3,1]);
Em_mean = mean(Em_post,'all');
% Em_q    = quantile(Em_post(:),0.75);
Em_IPCC = quantile(EM_post(:),0.66);
% i_subgroup = [1,5,8,21,25,39,46,55];
% % 1-4 Snowhill to Smith
% % 5-7 Gould Bay to Halley Bay
% % 8-20 Dawson to Kloa Point
% % 21-24 Fold Island to Cape Darnley
% % 25-38 Amanda Bay Point Geologie Davis Bay
% % 39-45 Ross Sea
% % 46-54 Amundsen Bellington
% n_subgroup = length(i_subgroup)-1;
em_post    = [];
Origin = {};
EM = [];
ORIGIN = [];
Col = {'StoS ','WEDD ','StoKP','MAWS ','AMPG ','ROSS ','A-B seas'};
for s = 1:n_subgroup
    ems = EM_post(i_subgroup(s):i_subgroup(s+1)-1,:);
    em_post{s,1} = ems(:);
    em_post{s,2} = [quantile(ems(:),[0.5,0.66,0.75]),mean(ems(:))];
    EM = [EM;ems(:)];
    origin = cell(length(em_post{s}),1);
    origin(:,1) = Col(s);
    Origin{s} = origin;
    ORIGIN = [ORIGIN;origin];
end

figure(1)
clf
% histogram(Em_post(:),0:0.02:1,'Normalization','pdf')
histogram(Em_post(:),0:0.05:1,'Normalization','pdf','FaceColor',Color(end,:))
hold on
line([Em_mean,Em_mean],[0,15],'color','k','linewidth',1)
% text(Em_mean+0.01,2,{'mean','$0.15$'},'Interpreter','latex','FontSize',16)

% line([Em_IPCC,Em_IPCC],[0,35],'color','k','linewidth',1)
% text(Em_IPCC+0.01,15,'$P_{66\%}$','Interpreter','latex','FontSize',16)
xlim([0,1])
ylabel({'Distribution of emigration rate','per colony over Antarctica'},'interpreter','latex','FontSize',16)
xlabel('Emigration rate','interpreter','latex','FontSize',16)

%% location of subpart on figure
xstart = .4;
xend   = .85;
ystart = .4;
yend   = .85;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
%% define range of subpart
vp_em = violinplot(Em_post(:),'Antartica');
vp_em(1,1).ViolinColor = Color(end,:);
vp_em.ShowData = 0;
line([0.5,1.5],[Em_mean,Em_mean],'color','k','linewidth',1)
% line([0.5,1.5],[Em_IPCC,Em_IPCC],'color','k','linewidth',1)
% text(0.53,Em_IPCC+0.08,'$0.15$','Interpreter','latex','FontSize',16)
ylabel('Emigration rate per colony','interpreter','latex','fontsize',14)
xlabel('Antarctica','interpreter','latex','fontsize',14)
ylim([0,0.3])



set(gca,'fontsize',10)

figure(2)
clf
vp = violinplot(EM,ORIGIN,'GroupOrder',Col);
line([0,8],[Em_mean,Em_mean],'color','k','linewidth',2)
% line([0,8],[Em_IPCC,Em_IPCC],'color','k','linewidth',2)
for i = 1:length(Col)
    vp(1,i).ViolinColor = Color(i,:);
end
ylabel('Emigration rate per colony','interpreter','latex','FontSize',16)
xlabel('Geographical regions of Antarctica','interpreter','latex','FontSize',16)
% set(gca,'fontsize',10)


for s = 1:n_subgroup
    figure(10+s)
    clf
    Is = i_subgroup(s):i_subgroup(s+1)-1;
    ems = EM_post(Is,:);
    origin = Is'*ones(1,size(EM_post,2));
    origin = {origin(:)};
    vpts = violinplot(ems(:),origin{:});
    for j = 1:length(Is)
    vpts(1,j).ViolinColor = Color(s,:);
    end
        line([0,length(Is)+1],[Em_mean,Em_mean],'color','k','linewidth',2,'linestyle','--')
%     line([0,length(Is)+1],[Em_IPCC,Em_IPCC],'color','k','linewidth',2,'linestyle','--')
    ylabel('Emigration rate per colony','interpreter','latex','FontSize',16)
    xlabel(Col(s),'interpreter','latex','fontsize',14)
end



% figure(3)
% clf
% vpt = violinplot(Em_post');
% for i = 1:ncol
%     s = sum(i>=i_subgroup);
%     vpt(1,i).ViolinColor = Color(s,:);
% end
% xlabel('EP colonies', 'Interpreter','latex')
% ylabel('Emigration rate per colony','interpreter','latex','FontSize',16)
% xlim([0,55])

