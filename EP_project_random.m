% clear all
% close all
% rand('state',sum(100.*clock)) % reset of the rand function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metapopulation model with density dependence and behaviors to move/ select
% colonies according to lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NNR,NNI,N] = EP_project_random(Tf,T0,Mu,D,pm,nsubpop)
%% Input 
% %%%%%%%%%%%%%%%%%%%%%%
% Tf = 2100;                 % Final time
% T0 = 1992;                 % Initial time
% ncol = 54; n_subpop = 4;
% mu = randfixedsum(n_subpop,1,1,0,1)';
% Mu = ones(ncol,1)*mu;  % Initial repartition of the sub group inside each colony
% D  = 1000;                 % Mean distance dispersal
% pm = 0.9*ones(ncol,1);     % mean proportion of migrant in each colony
% % over time
% jd = 0;                    % Choose 1) jd=0 Random strategy; 2) jd=1 Oriented strategy
% nsubpop = 4;               % Nb of subpop
%% Output
% N = sum(NN,2);              % n(i,j): # individuals in colony i at time t
% NNR = zeros(ncol,nsubpop,nt); % n(i,r,t)= # individuals of subgroup after reproduction
% NNI = zeros(ncol,ncol,nt);     % n(i,j,t)= # individuals that migrate
%%

%% Initialization %
%%%%%%%%%%%%%%%%%%
ncol = 54; % Nb of colonies
% nt   = 92; % Nb of time iteration from 2009 to 2100
time = 2009:2100;
%%%% calcul la distance cotiere entre chaque colonnie en km
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
dm = dispersion(d,ncol);   % Distance between colonies

%% Growth Function %
%%%%%%%%%%%%%%%%%%%
% With Density dependance K
% g: Individual growth rate f(N)/N
g = @(N,r,K) exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*(r<=0);

%*********************************
% jd = 0; % Choose 1) jd=0 Random strategy; 2) jd=1 Oriented strategy

%%%% Growth rate per colony / Carrying capacity and initial condition
load('Rmedian.mat')
% r = Rmedian(:,:,1);
it = sum((T0:Tf)<=2009);
r = [Rmedian(:,1,1)*ones(1,it-1),Rmedian(:,:,1)];

BE = xlsread('COL_EP.xlsx','F2:F55'); % count for each colony in 2009 from Fretwell et al. 
K  = 2*BE; % carrying capacity


%% Movement scenario %
%%%%%%%%%%%%%%%%%%%%%
% pm = 0.5;         % proportion pm
% movement = pm;    

%% Matrix of connexion %
%%%%%%%%%%%%%%%%%%%%%%%
% D = 1000;               % Mean distance dispersal
% Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zeroa
k = @(x) (x<D)/D .*(x>0); % Uniform kernel
 
% Matrix of connexion without renormalisation
dmat = k(dm);    % Unif kernel
% Renormalization of the connexion matrix %
den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
Dmat = dmat*diag(1./Den);
% Random Dispersal
disp_sel = Dmat;
den_disp = sum(disp_sel,1)';
Disp_sel = disp_sel-diag(den_disp); % pour assurer que la somme des colonnes vaut 0

%% Computation **********
%% Time
% Tf = 2014;
% it = 0;

% t = T0;
nt = length(T0:Tf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of size sub-population inside eachcolonies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN = zeros(ncol,nsubpop,nt);  % n(i,r,t)= # individuals of subgroup
NNR = zeros(ncol,nsubpop,nt); % n(i,r,t)= # individuals of subgroup after reproduction
NNd = zeros(ncol,nsubpop,nt); % n(i,r,t)= cumulated # dead individuals of subgroup
NNI = zeros(ncol,ncol,nt);     % n(i,j,t)= # individuals that migrate
% r\in(0,4) in colony i\in(1,46) at  ime t
% N = sum(NN,2);              % n(i,j): # individuals in colony i at time t

% mu = rand(ncol,n_subpop);
% Mu = mu./sum(mu,2);       % Mu(i,r): initial proportion of subpop r inside colony i
% NN0 = BE*ones(1,nsubpop).*Mu;   % Initial # individuals of each subpop in each colony
% 
% % Inidividuals
% NNnew = NN0;          % n(i,j,r)= # individuals of subgroup r\in(0,4) initially in colony j\in(1,46) inside colony i
% Nnew = sum(NNnew,2);  % N(i)  = # individuals in each colony


% NN(:,:,it) = NNnew;
Nnew = BE; 
%% Time loop backward
rold = r(:,1);
t = 2009;
while (t>T0)
    Nold = Nnew;
    Nnew = Nold./(1+rold);
    t = t-1;
    it = it-1;
end
N0 = Nnew;
% Gold = @(N) exp(log(1+rold).*(1-N./K)).*(rold>0) + (1+rold).*(rold<=0);
% Fold = @(N) Disp_sel*(pm.*Gold(N).*N)+Gold(N).*N;
% Gfold = @(N) (Disp_sel.*(pm')+diag(ones(1,ncol)))*diag(Gold(N));
% t = 2009;
% it = sum((T0:Tf)<=2009);
% nn = NN0;
% while (t>T0)
%     nnold = NNnew;
%     Nold = Nnew;
%     Disp = (Disp_sel.*(pm')+diag(ones(1,ncol)));
%     Nnew = Nold./(1+rold);
%     NNnew = (Disp*diag(Gold(Nnew)))\nnold;
%     NN(:,:,it-1) = NNnew;
%     N = [Nnew,N];
%     t = t-1;
%     it = it-1;
% end

%% Time loop forward
% Inidividuals
t = T0;
it = 0;
NN0 = N0*ones(1,nsubpop).*Mu;   % Initial # individuals of each subpop in each colony

% Inidividuals
NNnew = NN0;          % n(i,j,r)= # individuals of subgroup r\in(0,4) initially in colony j\in(1,46) inside colony i
Nnew = N0;  % N(i)  = # individuals in each colony

NN(:,:,it+1) = NNnew;
N = Nnew;
while (t<Tf)
    nnold = NNnew;
    Nold = Nnew;
    
    %%% Reproduction at time t+1
    rstar = g(Nold,r(:,it+1),K)-1; % Effective growth rate at t ime t

    %%% Exploration
%     lambda = movement; % mean proportion of indivudu that lives colonies           
    fn     = (1+rstar).*nnold;
    NNR(:,:,it+1) = fn;
    NNI(:,:,it+1) = Disp_sel*diag(pm.*sum(fn,2))+diag(sum(fn,2));
    
    NNnew  = Disp_sel*(pm.*fn) + fn;
    NN(:,:,it+1) = NNnew;
    Nnew = sum(NNnew,2);
%     Disp = (Disp_sel.*(pm')+diag(ones(1,ncol)));
%     Nnew = Disp*((1+rstar).*Nold);
    N = [N,Nnew];
    t = t+1;
    it = it+1;  
end
end


