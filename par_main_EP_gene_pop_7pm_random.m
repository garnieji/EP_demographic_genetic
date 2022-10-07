% clear all
% close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metapopulation model with density dependence and behaviors to move/ select
% colonies according to lambda
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = par_main_EP_gene_pop_7pm_random(Irmax,nD,irand)
format longE;
%%% Radom estimate 
jd = 0;
%% Input
Irmax = 1;% 1e4;  % # of iteration in mu
nD    = 1;    % Number of the output data
jd    = 0;    % Choice dispersal strategy 1) jd=0 Random strategy; 2) jd=1 Oriented strategy
irand = 0;    % Total random irand = 0 ou random outside the region sampled irand = 0

%% Parameters
ic    = 0;     % Choice of clustering
% if (jd == 1)
%     Nmax  = 1e3;   % # of iteration for (D,pm)
% elseif (jd==0)
    Nmax = 1e3;
% end
Dmax  = 6500;  % Upper bound on D
Dmin  = 250;   % Lower bound on D
% M     = 12;    % Number of workers (core used)

rand('state',sum(100.*clock)) % reset of the rand function

% %% Initialization %
% %%%%%%%%%%%%%%%%%%%
% load('allele_frequencies.mat')
% init_garnier
% BE(46)   = 0;                  % Colony 46 is extincted
% r(end,:) = -0.25;              % growth rate unkown
% dm       = dispersion(d,ncol); % Distance between colonies
% n_subpop = 4;                  % Nb of subpop
%
% icol_sampled = [24 30 20 22 33 35 4 7]; % Location of sampled colony
% % icol_sampled = [4 7 20 22 24 30 33 35]; % Location of sampled colony
%
% ncol_sampled = length(icol_sampled);    % Number of colonies sampled
% namecol_sampled = namecol(icol_sampled); % Sampled colonies % # sampled individuals = 220
%
% % Initialization of the parametre du modèle %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % mu = rand(ncol,n_subpop); % Initial repartition of the sub group inside each colony mu_ij\in(0,1)
% i_subgroup = [1,4,8,20,24,33,41,47];
% % 1-3 Snowhill to Smith
% % 4-7 Gould Bay to Halley Bay
% % 8-19 Stancomb to Kloa Point
% % 20-23 Fold Island to Cape Darnley
% % 24-32 Amanda Bay Pointe Geologie
% % 33-40 Ross Sea
% % 41-46 Amundsen Bellington
% n_subgroup = length(i_subgroup)-1;
% mu = rand(1,n_subgroup*n_subpop);% Initial repartition of the sub_group inside each colony

% ir = 0;
DDD = [];
PPm = [];
MU  = [];
llL = [];

tic
% parpool
for ir=1:Irmax
    % mu : Proportion in each sub region
    if (irand == 1)
        %% Total random repartition of subpop in colonies
        mu_StoS  = randfixedsum(4,1,1,0,1)';
        mu_WEDD  = randfixedsum(4,1,1,0,1)';
        mu_StoKP = randfixedsum(4,1,1,0,1)';
        mu_MAWS  = randfixedsum(4,1,1,0,1)';
        mu_AMPG  = randfixedsum(4,1,1,0,1)';
        mu_ROSS  = randfixedsum(4,1,1,0,1)';
        mu_AtoBe = randfixedsum(4,1,1,0,1)';
    else
        %% Clustering Laurent&co: Assumption = subpop per location sampled defined by the clustering algorithm
        mu_StoS = randfixedsum(4,1,1,0,1)';
        mu_WEDD  = [0,0,0,1];
        mu_StoKP = randfixedsum(4,1,1,0,1)';
        mu_MAWS = [0 1 0 0];
        mu_AMPG = [1 0 0 0];
        mu_ROSS = [0 0 1 0];
        mu_AtoBe = randfixedsum(4,1,1,0,1)';
    end
    mu = [mu_StoS mu_WEDD mu_StoKP mu_MAWS mu_AMPG mu_ROSS mu_AtoBe];
    
    DD = [];
    Pm = [];
    lL = [];
    % i = 0;
    %%% Minimizing lL with respect to d and pm
    %%%% Computation of lL(d,pm)
    tic
%     if (jd ==1)
%         parfor (i = 0:Nmax)
%             D  = (Dmax-Dmin)*rand(1)+Dmin; % Mean distance dispersal D\in(Dmin,Dmax)
%             pm = rand(1);                  % proportion pm\in(0,1)
%             DD = [DD,D]; Pm = [Pm,pm];
%             %     ll = LogL(mu,D,pm);
%             ll = LogL_real(mu,D,pm,jd,ic);
%             lL = [lL,ll];
%             %     i=i+1;
%         end
%     elseif(jd == 0)
        parfor (i = 0:Nmax)
            D  = (Dmax-Dmin)*rand(1)+Dmin; % Mean distance dispersal D\in(Dmin,Dmax)
            pm = rand(7,1);                  % proportion pm\in(0,1)
            DD = [DD,D]; Pm = [Pm,pm];
            %     ll = LogL(mu,D,pm);
            ll = LogL_real_7pm(mu,D,pm,jd,ic);
            lL = [lL,ll];
            %     i=i+1;
        end
%     end
    %     [lL_min,I] = min(lL);
    [lL_min,I] = mink(lL,100);
    toc
    
%     %%%% Minimizing with fmincon
%     tic
%     fun = @(X) LogL_real_7pm(mu,X(1),X(2:8),jd,ic);
%     lb = [Dmin,zeros(1,7)];
%     ub = [Dmax,ones(1,7)];
%     A = []; b = []; Aeq = []; beq = []; nonlcon = [];
%     X0 = [(Dmax-Dmin)*rand(1)+Dmin,rand(1,7)];
%     [X_min,lL_min] = fmincon(@(X) fun(X),X0,A,b,Aeq,beq,lb,ub);
%     toc
    
    % figure(1)
    % clf
    % hold on
    % plot3(DD,Pm,lL,'*')
    % plot3(DD(I),Pm(I),lL_min,'r*')
    
    DDD = [DDD;DD(I)];
    PPm = [PPm;Pm(:,I)];
    MU  = [MU;mu];
    llL = [llL;lL(I)];
    
    if (mod(ir,50)==0)
        save(['par_main_EP_gene_7pm_pop_real_jd',num2str(jd),'_L_',num2str(nD),'_rand',num2str(irand)],'MU','DDD','PPm','llL','Dmin','Dmax','Irmax')
    end
    
%     save('Test_par','MU','DD','lL','Pm','DDD','PPm','llL','Dmin','Dmax','Irmax')
end
% toc
end
% plot3(DD(I),Pm(I),lL_min,'r*')

% % Minimization of the -log(L)(mu,D,pm)
% % Initialization
% X0 = [mu,D,pm];
% % Equality constraint for Mu
% Aeq = zeros(n_subgroup,length(X0));
% for i = 1:n_subgroup
%     Aeq(i,(i-1)*n_subpop+1 : i*n_subpop ) = 1;
% end
% beq = ones(n_subgroup,1);
%
% Mmu = X0*(Aeq'*Aeq);
% mmu = mu./Mmu(1:end-2);
% X0 = [mmu,D,pm];
%
% % Constraint inequality
% Mu_min = zeros(1,n_subgroup*(n_subpop));
% Mu_max = ones(1,n_subgroup*(n_subpop));
%
% LU = [Mu_min,1000,0];
% LM = [Mu_max,6000,1];
% ll = @(X) LogL(X(1:end-2),X(end-1),X(end));
%
% options = optimoptions('fmincon','MaxFunctionEvaluation',6000,'MaxIterations',2000);
% [X,Fval] = fmincon(@(X) LogL(X(1:end-2),X(end-1),X(end)),X0,[],[],Aeq,beq,LU,LM,[], ...
% options);

% % Minimization of the -log(L)(D,pm)
% % Initialization
% D  = 6250*rand(1)+1500;    % Mean distance dispersal D\in(250,6500)
% pm = rand(1);             % proportion pm\in(0,1)
% X0 = [mu,D,pm];
% % Equality constraint for Mu
% Aeq = zeros(n_subgroup,length(X0));
% for i = 1:n_subgroup
%     Aeq(i,(i-1)*n_subpop+1 : i*n_subpop ) = 1;
% end
% beq = ones(n_subgroup,1);
%
% Mmu = X0*(Aeq'*Aeq);
% mmu = mu./Mmu(1:end-2);
% X0 = [D,pm];
%
% % Constraint inequality
% Mu_min = zeros(1,n_subgroup*(n_subpop));
% Mu_max = ones(1,n_subgroup*(n_subpop));
%
% LU = [1500,0];
% LM = [6000,1];
% ll = @(X) LogL(X(1:end-2),X(end-1),X(end));
%
% options = optimoptions('fmincon','MaxFunctionEvaluation',6000,'MaxIterations',2000);
% [X,Fval] = fmincon(@(X) LogL(mu,X(end-1),X(end)),X0,[],[],[],[],LU,LM,[], ...
% options);
%
% plot3(X(1),X(2),Fval,'g*')


%% Figure test
% figure(1)
% clf
% plot(N,'r')
%
%
% figure(2)
% plot(cumsum(prop,2))
% hold on
% plot((r(:,Tf)>0),'*-')













