
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computation of log(L) log Likelihood with incorporation of dead and
%  alive sampled penguins
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lL] = LogL_real_7pm_informed(mu,D,pm,jd,ic)
% % Parametre du probleme
%%%%%%%%%%%%%%%%%%%%%%
% ncol = 54; nsubpop = 4;
% mu = ones(7,1)*randfixedsum(nsubpop,1,1,0,1)'; % Initial repartition of the sub group inside each colony
% D  = 1000;                 % Mean distance dispersal
% pm = 0.5*ones(7,1);                  % proportion pm
% jd = 0;                    % Choose 1) jd=0 Random strategy; 2) jd=1 Oriented strategy
% ic = 0;                    % Choose 1) ic = 0 pop clustering (Younger, 2017); 
%                                     2) ic = 1 indiv clustering Laurent
% Initialization %
%%%%%%%%%%%%%%%%%%
if (ic == 0)
    %% Clustering Younger (2017)
    load('allele_frequencies_new.mat')
    load('genotype_individuals_new')
else
    %% Clustering Laurent &co
    load genotype_individuals_only.mat
    load allele_freq_subpop.mat
    load proba_genotype_PPr.mat
end
%%
% Parametre du modèle
%%%%%%%%%%%%%%%%%%%%%%
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
% mu = rand(1,n_subgroup*n_subpop); % Initial repartition of the sub_group inside each colony 
Mu = [];
Pm = [];
for s = 1:n_subgroup
    Is = i_subgroup(s+1)-i_subgroup(s);
    mus = mu( (s-1)*nsubpop+1 : s*nsubpop );
    Mus = mus; % mus/sum(mus);
    Mu = [Mu;ones(Is,1)*Mus];
    Pm = [Pm;ones(Is,1)*pm(s)];
end

[NNR,NNI,N] = EP_project_informed(Tf,T0,Mu,D,Pm,nsubpop,jd);

lL = 0;
%% Sampling blood on alive penguins
for T_sampled = [1992,1993,2010,2012,2013]% [2010,2012,2013]    % Sampling time
    iT_sampled = T_sampled-T0+1;
    tau = icol_sampled(t_sampled==T_sampled); 
    P_i_to_h = NNI(tau,:,iT_sampled)./N(tau,iT_sampled);
    Mur = NNR(:,:,iT_sampled)./(sum(NNR(:,:,iT_sampled),2)+(sum(NNR(:,:,iT_sampled),2)<=0));
    P_Gitau_in_h = Mur*PPr';        % Allele frequency without error
    lLt = 0;
    i_tau = (iPPr==tau);
    for j = 1:length(tau)
        PPtau = P_i_to_h(j,:)*P_Gitau_in_h(:,i_tau(:,j));
        lLt = lLt + sum(log(PPtau));
    end
    lL = lL-lLt;
end
% %% Sampling on dead penguins
% for T_sampled = [2010,2012]
%     iT_sampled = T_sampled-2009+1;
%     tau = icol_sampled(t_sampled==T_sampled);
%     lLt = 0;
%     %% Chicks die the year of capture
%     iT = iT_sampled;
%     P_i_to_h = NNI(tau,:,iT)./N(tau,iT);
%     Mur = NNR(:,:,iT)./(sum(NNR(:,:,iT),2)+(sum(NNR(:,:,iT),2)<=0));
%     P_Gitau_in_h = Mur*PPr';        % Allele frequency without error 
%     i_tau = iPPr==tau;
%     for j = 1:length(tau)         
%         %% Chicks may have died before the year of capture
% %         PPtau = 0;
% %         for iT = 2:iT_sampled
% %             P_i_to_h = NNI(tau,:,iT)./N(tau,iT);
% %             Mur = NNR(:,:,iT)./(sum(NNR(:,:,iT),2)+(sum(NNR(:,:,iT),2)<=0));
% %             P_Gitau_in_h = Mur*PPr';        % Allele frequency without error 
% %             i_tau = iPPr==tau;
% %             PPtau = PPtau + P_i_to_h(j,:)*P_Gitau_in_h(:,i_tau(:,j)).*(-rlo(tau(j),iT+1));
% %          end
%         PPtau = P_i_to_h(j,:)*P_Gitau_in_h(:,i_tau(:,j)).*(abs(rlo(tau(j),iT+1)));
%         lLt =lLt + sum(log(PPtau));
%     end
%     lL  = lL-lLt;
% end

% 
%     % Proportion of subgroup
%     np = sum(nn,2);
%     npp = sum(np,3);
%     prop_subpop = zeros(ncol,n_subpop);
%     for i=1:n_subpop
%         prop_subpop(:,i) = np(:,1,i)./npp;
%     end
%     % Proportion coming from initial colonies
%     prop = n./N;
%     Prop = prop(icol_sampled,:);
%     SmPPr = Prop*Mu*PPr';
%     % Computation of the log of Likelyhood function -log(L)
%     lL = -sum(sum(log(SmPPr))) ;
end
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













