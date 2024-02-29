clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];
Ir = [7,3,6];     % Number of iteration process
Ir_rand = 2; % Number of iteration process parallel
ir = 0;
choice = 1; % 0 = full random / 1 = informed disp - random search / 2 = full informed

Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS','A-B seas'};

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

n_post = length(D);
N_data = length(D);
I_post = randsample(N_data,n_post,true,wm);

D_post = D(I_post);
Pm_post = Pm(I_post,:);
Mu_post = Mu(I_Mu(I_post),:);

D_med = median(D_post);
Pm_med = median(Pm_post);
Mu_mean = mean(Mu_post);
%% Profile likelihood function

%% theta = D
dD = (Dmax-Dmin)/500;
dd = Dmin:dD:Dmax;
nd = length(dd);
h_D  = zeros(1,nd);
parfor i = 1:nd
    Id = logical((D<=dd(i)+dD).*(D>dd(i)-dD));
    hi = min(lL(Id));
    if (isempty(hi)==0)
        h_D(i) = hi;
    end
end
Ih = 2*(h_D - lL_min)<=3.84;
dh_D = dd(Ih);

% %% theta = Pm
n_subgroup = 7;
n_subpop   = 4; 
np = 101;
ppm = linspace(0,1,np);
h_Pm  = zeros(n_subgroup,np-1);
dh_pm = cell(1,n_subgroup);
for j = 1:n_subgroup
    parfor i = 1:nd-1
        Id = logical((Pm(:,j)<=ppm(i+1)).*(Pm(:,j)>ppm(i)));
        hi = min(lL(Id));
        if (isempty(hi)==0)
            h_Pm(j,i) = hi;
        end
    end
    Ih = 2*(h_Pm(j,:) - lL_min)<=3.84;
    dh_pm{j} = ppm(Ih);
end

%% theta = Mum
nm = 101;
mum = linspace(0,1,nm);
h_Mum  = zeros(n_subgroup*n_subpop,nm-1);
dh_mu = cell(1,n_subgroup*n_subpop);
for j = 1:n_subgroup*n_subpop
    parfor i = 1:nm-1
        Id = logical((Mu(:,j)<=mum(i+1)).*(Mu(:,j)>mum(i)));
        hi = min(lL(Id));
        if (isempty(hi)==0)
            h_Mum(j,i) = hi;
        end
    end
    Ih = 2*(h_Mum(j,:) - lL_min)<=3.84;
    dh_mu{j} = mum(Ih);
end

ALPHA = (mean(Mu_post).*(1+mean(Mu_post))./var(Mu_post)- 1 ).*mean(Mu_post);

%% Goodness of fit
load genotype_individuals.mat

% Genotype of each individuals
% load Individuals_genotypes.mat
M_CG = (M>1).*(M<4);         % 0 if A-T and 1 if C-G and 0 if missing
M_missing = (M<0);           % 0 if A-T and 0 if C-G and 1 if missing (-9)
M_AT = 1 - M_CG - M_missing; % 1 if A-T and 0 if C-G and 0 if missing
K_AT = M_AT(1:2:end,:).*M_AT(2:2:end,:);
K_CG = M_CG(1:2:end,:).*M_CG(2:2:end,:);
K    = K_AT +K_CG;
k = sum(K,2);
n_obs = size(M_CG,1)/2;
% SNP_ID with high variance
Var_AT = var(Freq_subpop_AT);  % Inter subpop variance of each SNP_ID
[Var_AT_s,Is] = sort(Var_AT);  % Sort the variance
I_SNP = Is(end-20:end);        % Take the 20 SNP_ID that maximize the variance
ki = sum(K(:,I_SNP),2);

nsubpop = 4;
nSNP = 21;
M = zeros(2^nSNP,nSNP);
for i=1:nSNP
    x = [zeros(2^(i-1),1);ones(2^(i-1),1)];
    M(:,i) = repmat(x,2^(nSNP-i),1);
end
nx = 1e3;
px = nx;
x = lhsdesign(nx,px);
y = lhsdesign(nx,px);
X = ceil(nSNP*x);
Y = ceil(nSNP*y);

%% Parametre du modèle
%%%%%%%%%%%%%%%%%%%%%%
T_sampled = [1992,1993,2010,2012,2013];
nt = length(T_sampled);
Tf = 2014; % Final time
T0 = 1991; % Final initial

%%%% mu : Initial repartition of the sub group inside each colony
% IL FAUT CLASSER LES COLONIES EN n_subgroup GROUP
i_subgroup = [1,5,8,21,25,39,46,55];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-20 Dawson to Kloa Point
% 21-24 Fold Island to Cape Darnley
% 25-38 Amanda Bay Point Geologie Davis Bay
% 39-45 Ross Sea
% 46-54 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;
nsubpop = 4;
choice_m = 'MED';
if (choice_m == 'MEL')
    %% More Likely Estimate
    D_GF = D_est;
    Mu_GF = [];
    Pm_GF = [];
    for s = 1:n_subgroup
        Is = i_subgroup(s+1)-i_subgroup(s);
        mus = Mu_est( (s-1)*nsubpop+1 : s*nsubpop );
        Mus = mus;
        Mu_GF = [Mu_GF;ones(Is,1)*Mus];
        Pm_GF = [Pm_GF;ones(Is,1)*Pm_est(s)];
    end
elseif (choice_m == 'MED')
    %% Median of estimations
    D_GF = D_med;
    Mu_GF = [];
    Pm_GF = [];
    for s = 1:n_subgroup
        Is = i_subgroup(s+1)-i_subgroup(s);
        mus = Mu_mean( (s-1)*nsubpop+1 : s*nsubpop );
        Mus = mus;
        Mu_GF = [Mu_GF;ones(Is,1)*Mus];
        Pm_GF = [Pm_GF;ones(Is,1)*Pm_med(s)];
    end
end
%% Compute outcome of the model for the defined parameters
if (choice == 0)
    [NNR,NNI,N] = EP_project_random(Tf,T0,Mu_GF,D_GF,Pm_GF,nsubpop);
elseif (choice==1)
    jd = 0;
    [NNR,NNI,N] = EP_project_informed_old(Tf,T0,Mu_GF,D_GF,Pm_GF,nsubpop,jd);
elseif (choice == 2)
    jd = 1;
    [NNR,NNI,N] = EP_project_informed_old(Tf,T0,Mu_GF,D_GF,Pm_GF,nsubpop,jd);
end


PPr = zeros(nx*px,nsubpop);
PPr_obs = zeros(n_obs,nsubpop);
for i = 1:nsubpop
    %% Possible outcomes
    k   = sum(M(X(:),:).*M(Y(:),:) + (1-M(X(:),:)).*(1-M(Y(:),:)),2);
    pr1 = M(X(:),:).*Freq_subpop_CG(i,I_SNP) + (1-M(X(:),:)).*Freq_subpop_AT(i,I_SNP);
    pr2 = M(Y(:),:).*Freq_subpop_CG(i,I_SNP) + (1-M(Y(:),:)).*Freq_subpop_AT(i,I_SNP);
    Pr = pr1.*pr2;%.*(2-K(:,I_SNP));
    PPr(:,i) = prod(Pr,2).*2.^k;
    %% Observations
    pr_obs = M_CG(:,I_SNP).*Freq_subpop_CG(i,I_SNP) + M_AT(:,I_SNP).*Freq_subpop_AT(i,I_SNP) + M_missing(:,I_SNP);
    Pr_obs = pr_obs(1:2:end,:).*pr_obs(2:2:end,:);%.*(2-K(:,I_SNP));
    PPr_obs(:,i) = prod(Pr_obs,2).*2.^ki;
end
% Index
nind_sampled = [16,15,16,16,10,11,13,13];
Is = [0,cumsum(nind_sampled)];
iPPr_obs = zeros(n_obs,1);
for i = 1:length(icol_sampled)
    iPPr_obs(Is(i)+1:Is(i+1)) = icol_sampled(i);
end

%% Sampling blood on alive penguins
PL = cell(1,nt);
PL_obs = cell(1,nt);
for it = 1:nt  % Sampling time
    iT_sampled = T_sampled(it)-T0+1;
    tau = icol_sampled(t_sampled==T_sampled(it));
    % P_i_to_h = NNI(tau,:,iT_sampled)./N(tau,iT_sampled);
    Mur = NNR(:,:,iT_sampled)./(sum(NNR(:,:,iT_sampled),2)+(sum(NNR(:,:,iT_sampled),2)<=0));
    %% Possible outcomes
    P_Gitau_in_h = Mur(tau,:,:)*PPr';        % Allele frequency without error
    PL{it} = P_Gitau_in_h;
    %% Observations
    for j = 1:length(tau)
        It = (iPPr_obs==tau(j));
        P_Gitau_in_h_obsj = Mur(tau(j),:,:)*PPr_obs(It,:)';
        PL_obs{it}{j} = P_Gitau_in_h_obsj;
    end
end

%% Figure 95% regions
n_95 = round(0.05*(nx*px));
Per_95 = [];
p_95 = 0;
for it = 1:nt
    It = T_sampled(it)-T0+1;
    tau = icol_sampled(t_sampled==T_sampled(it));
    per_95 = cell(1,length(tau));
    for j= 1:length(tau)
        PL_sort = sort(PL{it}(j,:));
        PL_95 = PL_sort(n_95);
        np_95 = sum(PL_obs{1,it}{j}>=PL_95);
        per_95{1,j} = np_95./length(PL_obs{1,it}{j})*100;
        p_95 = p_95 + np_95;
    end
    Per_95 = [Per_95,per_95{:}];
end
P_95 = p_95/n_obs*100;






