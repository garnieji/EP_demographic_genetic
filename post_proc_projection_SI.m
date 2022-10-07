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

%% Mu fixed in colonies sampled
%% Choice = 1 SI
load('post_proc_EP_7pm_informed_jd0')

Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS',{'A-B', 'seas'}};
NN = [];
Tf = 2100;
T0 = 2009;
parfor i = 1:length(D_post)
    [A_subgroup,A_col,N] = EP_project_informed_7col_migration(Tf,T0,D_post(i),Pm_post(i),0);
    NN = [NN;sum(N)];
end
