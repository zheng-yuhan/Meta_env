clear
%% Load gut samples from HMP, genus taxonomoc level, relative abundance, unassigned reads are excluded
%A=load('HMP_genus.txt');
%[f,p]=uigetfile({'*.xlsx'});
%%%[A,~,~]=xlsread([p,f],'L6');
%[A,~,~]=xlsread([p,f]);
MKSZ=3;
[f,p]=uigetfile({'*.txt'});
 T=readtable([p,f]);
% T=T(:,2:end);
% A=table2array(T);clear T;
%T=readtable('D:\Temp\jxsynwq\nwqv.txt');
% %T=readtable("F:\dzl\Temp\KO\2018ko.txt");T=T(:,8:13);
 % [filename,pathname]=uigetfile({'*.txt'});
 % T=readtable([pathname,filename]);figure (2)
c14_7=[  0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

c14_light_02=1-0.5*(1-c14_7);


hold on

% real
mean_BS=mean(BS,1);
Prctile_BS = prctile(BS,[10  90],1)'; 
h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ]);
h_area(1).FaceColor = 'none';
h_area(2).FaceColor = c14_light_02(6,:);
h_area(1).EdgeColor = 'none';
h_area(2).EdgeColor = 'none';

h_point=plot(Overlap(:),RootJSD(:),'.','color',[0.5 0.5 1],...
    'MarkerSize',MKSZ);

plot(xs,mean_BS','linewidth',1,'color',c14_7(1,:))
plot(xs,Prctile_BS(:,1),'linewidth',0.7,'color',c14_7(6,:))
plot(xs,Prctile_BS(:,2),'linewidth',0.7,'color',c14_7(6,:))

xc=detect_negative_slope_Lowess_03(xs,mean_BS');
plot(xc*[1 1],[0 1],'-g')


% xlim([0.5 1])
% ylim([0.05 0.6])


% box on
% axis square

% set(gca,'fontsize',18)
% set(gca,'xtick',0:0.1:1)
% set(gca,'ytick',0.1:0.1:0.5)

xlabel('Overlap','fontsize',12)
ylabel('Dissimilarity','fontsize',12)
%title({'HMP, stool, genus, '; [num2str(size(BS,1)) ' bootstrap realization(s)']},    'fontsize',12)


hold on

% null
mean_BS_null=mean(BS_null,1);
Prctile_BS_null = prctile(BS_null,[10  90],1)'; 
h_area_null=area(xs_null',[Prctile_BS_null(:,1), (Prctile_BS_null(:,2)-Prctile_BS_null(:,1)) ]);
h_area_null(1).FaceColor = 'none';
h_area_null(2).FaceColor = c14_light_02(3,:);
h_area_null(1).EdgeColor = 'none';
h_area_null(2).EdgeColor = 'none';

h_point_null=plot(Overlap_null(:),RootJSD_null(:),'.','color',[1 0.5 0.5],...
    'MarkerSize',MKSZ);


plot(xs_null,mean_BS_null','linewidth',1,'color',c14_7(2,:))
plot(xs_null,Prctile_BS_null(:,1),'linewidth',0.7,'color',c14_7(3,:))
plot(xs_null,Prctile_BS_null(:,2),'linewidth',0.7,'color',c14_7(3,:))

xsnull=detect_negative_slope_Lowess_03(xs_null,mean_BS_null');
Fns_null=sum(Overlap_null(:)>xsnull)/sum(~isnan(Overlap_null(:)));

plot(xsnull*[1 1],[0 1],'--','color',[0.9 0.5 0.2])

text(0.95,0.95,string())
%xlim([0.5 1])
%xlim([0.92 1])
ylim([0 1])
hold off
%box on
axis square
text(0.95,0.95,{'Fns=',string(Fns),'p=',string(P_value_slopes)});

%set(gca,'fontsize',18)
%set(gca,'xtick',0:0.1:1)
%set(gca,'ytick',0.1:0.1:0.5)

%xlabel('Overlap','fontsize',22)
%ylabel('Dissimilarity','fontsize',22)

%title({'Randomized samples'; [num2str(size(BS,1)) ' bootstrap realization(s)']},...
 %   'fontsize',22)


%pause(1)
%h_marker=h_point_null.MarkerHandle;  h_marker.EdgeColorData(4)=100;
%h_marker=h_point.MarkerHandle;  h_marker.EdgeColorData(4)=100;

alpha 0.3
%T=readtable("H:\TEMP\2018HB-v.txt");
% T=readtable('H:\TEMP\HB-Ko.txt');
 T=T(:,2:end);
A=table2array(T);clear T;
%while 1  
%% normalize to relative abundances
 [NumSpecies,NumSamples]=size(A);
 A=A./repmat(sum(A),NumSpecies,1);

%% make null model
 A_null=func_make_null(A,1);

%% get Overlap and Dissimilarity for the real samples
[Overlap,RootJSD]=func_Cal_Overlap_rJSD_from_relative_abundance(A);

%% get Overlap and Dissimilarity for the null model
 [Overlap_null,RootJSD_null]=func_Cal_Overlap_rJSD_from_relative_abundance(A_null);

%% make Bootstrap for the real samples
%N_times=1000;
N_times=1000;
[BS, xs]=func_cal_rlowess_bootstrap(Overlap,RootJSD,N_times);

%% make Bootstrap for the null model
[BS_null, xs_null]=func_cal_rlowess_bootstrap(Overlap_null,RootJSD_null,N_times);

%% plot figure
c14_7=[  0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

c14_light_02=1-0.5*(1-c14_7);


hold on

% real
mean_BS=mean(BS,1);
Prctile_BS = prctile(BS,[1  99],1)'; 
h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ]);
h_area(1).FaceColor = 'none';
h_area(2).FaceColor = c14_light_02(6,:);
h_area(1).EdgeColor = 'none';
h_area(2).EdgeColor = 'none';

h_point=plot(Overlap(:),RootJSD(:),'.','color',[0.5 0.5 1],...
    'MarkerSize',10);

plot(xs,mean_BS','linewidth',1,'color',c14_7(1,:))
plot(xs,Prctile_BS(:,1),'linewidth',0.7,'color',c14_7(6,:))
plot(xs,Prctile_BS(:,2),'linewidth',0.7,'color',c14_7(6,:))

xc=detect_negative_slope_Lowess_03(xs,mean_BS');
plot(xc*[1 1],[0 1],'-g')


% xlim([0.5 1])
% ylim([0.05 0.6])


% box on
% axis square

% set(gca,'fontsize',18)
% set(gca,'xtick',0:0.1:1)
% set(gca,'ytick',0.1:0.1:0.5)

xlabel('Overlap','fontsize',12)
ylabel('Dissimilarity','fontsize',12)
%title({'HMP, stool, genus, '; [num2str(size(BS,1)) ' bootstrap realization(s)']},    'fontsize',12)


hold on

% null
mean_BS_null=mean(BS_null,1);
Prctile_BS_null = prctile(BS_null,[10  90],1)'; 
h_area_null=area(xs_null',[Prctile_BS_null(:,1), (Prctile_BS_null(:,2)-Prctile_BS_null(:,1)) ]);
h_area_null(1).FaceColor = 'none';
h_area_null(2).FaceColor = c14_light_02(3,:);
h_area_null(1).EdgeColor = 'none';
h_area_null(2).EdgeColor = 'none';

h_point_null=plot(Overlap_null(:),RootJSD_null(:),'.','color',[1 0.5 0.5],...
    'MarkerSize',MKSZ);


plot(xs_null,mean_BS_null','linewidth',1,'color',c14_7(2,:))
plot(xs_null,Prctile_BS_null(:,1),'linewidth',0.7,'color',c14_7(3,:))
plot(xs_null,Prctile_BS_null(:,2),'linewidth',0.7,'color',c14_7(3,:))


%xlim([0.5 1])
ylim([0 1])
hold off
%box on
%axis square

%set(gca,'fontsize',18)
%set(gca,'xtick',0:0.1:1)
%set(gca,'ytick',0.1:0.1:0.5)

%xlabel('Overlap','fontsize',22)
%ylabel('Dissimilarity','fontsize',22)

%title({'Randomized samples'; [num2str(size(BS,1)) ' bootstrap realization(s)']},...
 %   'fontsize',22)


%pause(1)
%h_marker=h_point_null.MarkerHandle;  h_marker.EdgeColorData(4)=100;
%h_marker=h_point.MarkerHandle;  h_marker.EdgeColorData(4)=100;

alpha 0.3

%%
% Fraction of data with negative slope
Fns=sum(Overlap(:)>xc)/sum(~isnan(Overlap(:)));

% P-value for the slope to the right of the changing point
slope_BS_real=nan(N_times,1);
slope_BS_null=nan(N_times,1);
parfor k=1:N_times
    P_real=polyfit(xs(xs>xc),BS(k,xs>xc),1);
    slope_BS_real(k)=P_real(1);

    P_null=polyfit(xs_null(xs>xc),BS_null(k,xs>xc),1);
    slope_BS_null(k)=P_null(1);
end
[~,P_value_slopes]=ttest2(slope_BS_real,slope_BS_null);

% Display the universality scores
disp(['Fns=' num2str(Fns)])
disp(['p-value=' num2str(P_value_slopes)])

% if P_value_slopes<0.05
%     break
% end
% end

 figure (2)
c14_7=[  0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

c14_light_02=1-0.5*(1-c14_7);


hold on

% real
mean_BS=mean(BS,1);
Prctile_BS = prctile(BS,[10  90],1)'; 
h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ]);
h_area(1).FaceColor = 'none';
h_area(2).FaceColor = c14_light_02(6,:);
h_area(1).EdgeColor = 'none';
h_area(2).EdgeColor = 'none';

h_point=plot(Overlap(:),RootJSD(:),'.','color',[0.5 0.5 1],...
    'MarkerSize',MKSZ);

plot(xs,mean_BS','linewidth',1,'color',c14_7(1,:))
plot(xs,Prctile_BS(:,1),'linewidth',0.7,'color',c14_7(6,:))
plot(xs,Prctile_BS(:,2),'linewidth',0.7,'color',c14_7(6,:))

xc=detect_negative_slope_Lowess_03(xs,mean_BS');
plot(xc*[1 1],[0 1],'-g')


% xlim([0.5 1])
% ylim([0.05 0.6])


% box on
% axis square

% set(gca,'fontsize',18)
% set(gca,'xtick',0:0.1:1)
% set(gca,'ytick',0.1:0.1:0.5)

xlabel('Overlap','fontsize',12)
ylabel('Dissimilarity','fontsize',12)
%title({'HMP, stool, genus, '; [num2str(size(BS,1)) ' bootstrap realization(s)']},    'fontsize',12)


hold on

% null
mean_BS_null=mean(BS_null,1);
Prctile_BS_null = prctile(BS_null,[10  90],1)'; 
h_area_null=area(xs_null',[Prctile_BS_null(:,1), (Prctile_BS_null(:,2)-Prctile_BS_null(:,1)) ]);
h_area_null(1).FaceColor = 'none';
h_area_null(2).FaceColor = c14_light_02(3,:);
h_area_null(1).EdgeColor = 'none';
h_area_null(2).EdgeColor = 'none';

h_point_null=plot(Overlap_null(:),RootJSD_null(:),'.','color',[1 0.5 0.5],...
    'MarkerSize',MKSZ);


plot(xs_null,mean_BS_null','linewidth',1,'color',c14_7(2,:))
plot(xs_null,Prctile_BS_null(:,1),'linewidth',0.7,'color',c14_7(3,:))
plot(xs_null,Prctile_BS_null(:,2),'linewidth',0.7,'color',c14_7(3,:))

xsnull=detect_negative_slope_Lowess_03(xs_null,mean_BS_null');
Fns_null=sum(Overlap_null(:)>xsnull)/sum(~isnan(Overlap_null(:)));

plot(xsnull*[1 1],[0 1],'--','color',[0.9 0.5 0.2])

text(0.95,0.95,string())
%xlim([0.5 1])
%xlim([0.92 1])
ylim([0 1])
hold off
%box on
axis square
text(0.95,0.95,{'Fns=',string(Fns),'p=',string(P_value_slopes)});

%set(gca,'fontsize',18)
%set(gca,'xtick',0:0.1:1)
%set(gca,'ytick',0.1:0.1:0.5)

%xlabel('Overlap','fontsize',22)
%ylabel('Dissimilarity','fontsize',22)

%title({'Randomized samples'; [num2str(size(BS,1)) ' bootstrap realization(s)']},...
 %   'fontsize',22)


%pause(1)
%h_marker=h_point_null.MarkerHandle;  h_marker.EdgeColorData(4)=100;
%h_marker=h_point.MarkerHandle;  h_marker.EdgeColorData(4)=100;

alpha 0.3