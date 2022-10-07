clear all
close all
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing minimization Log Likelyhood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color  = get(gca,'colororder');
Colori = [Color(2,:);[0,0,0];Color(1,:);Color(3,:);Color(5,:)];
Marker = ['o','*','d','^'];
Ir = [7,3,2];     % Number of iteration process
Ir_rand = 2; % Number of iteration process parallel
ir = 0;

for choice = 1% 0:2 % 0 = full random / 1 = informed disp - random search / 2 = full informed

%% Mu fixed in colonies sampled
if (choice == 0)
    load('post_proc_EP_7pm_random')
elseif (choice == 1)
    load('post_proc_EP_7pm_informed_jd0_connectivity')
elseif (choice == 2)
    load('post_proc_EP_7pm_informed_jd1')
end
Col = {'StoS','WEDD','StoKP','MAWS','AMPG','ROSS',{'A-B', 'seas'}};
i_subgroup = [1,5,8,21,25,39,46,55];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-20 Dawson to Kloa Point
% 21-24 Fold Island to Cape Darnley
% 25-38 Amanda Bay Point Geologie Davis Bay
% 39-45 Ross Sea
% 46-54 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;
lat   = xlsread('COL_EP.xlsx','C2:C55'); % count for each colony in 2009 from Fretwell et al.
long  = xlsread('COL_EP.xlsx','D2:D55');

ncol=length(lat);
lat  = lat*pi/180;
long = long*pi/180;

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

    im = ceil(length(i_subgroup(i):i_subgroup(i+1)-1)/2);
    X = [X,xi(im)];
    Y = [Y,yi(im)];

end
%%% Relocation of text node
Xtext = X; Ytext = Y;
Xtext(1) = X(1)+100; Ytext(1) = Y(1)+200;
Xtext(2) = X(2)+200;
Xtext(3) = X(3)+100; Ytext(3) = Y(3)+200;
Xtext(4) = X(4)+500; Ytext(4) = Y(4)+300;
Xtext(5) = X(5)-400; Ytext(5) = Y(5)-1400;
Xtext(6) = X(6)-1200; Ytext(6) = Y(6)-200;
Xtext(7) = X(7)-1400; Ytext(7) = Y(7)+500;


%% Regarder les violinplot des A_subgroup!!
na =length(D_post);
a_sub = zeros(7,7,5*na);
for ik = 1:na
    a_sub(:,:,5*ik-4:5*ik) = A_subgroup{ik};
end
A_sub_med = median(a_sub,3);

Ga_subgroup = digraph(100*A_sub_med);
LaaWidths = 5*Ga_subgroup.Edges.Weight/max(Ga_subgroup.Edges.Weight);
Ga_subgroup.Edges.Weight = round(Ga_subgroup.Edges.Weight,2);
figure(choice+1)
clf
hold on
% plot([x;x(1)],[y;y(1)],'k-')
for i=1:7
     plot(Xi{1,i},Yi{1,i},'-','linewidth',1,'color',Color(i,:))

     plot(Ga_subgroup,'XData',X,'YData',Y,'NodeLabel',{},...
      'EdgeLabel',Ga_subgroup.Edges.Weight,'LineWidth',LaaWidths,'ArrowSize',10,...
      'Nodefontsize',12,'EdgeFontWeight','bold')
 for i =1:7
     plot(X(i),Y(i),'o','MarkerSize',7,'MarkerFaceColor',Color(C(i),:),'color',Color(i,:),'LineWidth',1)
 end
 axis([-3500,4500,-2500,3500])
% axis([-3500,4500,-3500,4500])
drawnow
end
