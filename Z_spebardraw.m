%前10物种柱状图
%% 数据读取与处理
clear
[filename,pathname]=uigetfile({'\*.txt'});
 microb=readtable([pathname,filename]); 
 T1=microb.Properties.VariableNames;
 ta=microb(:,1);ta=table2cell(ta);
 T=microb(:,2:end);
spe=table2array(T);clear T;
%% 前十筛选
  s_spe=sum(spe);st_spe=sum(spe'); 
  m_spe=spe./s_spe;
[B,I]=sort(st_spe,'descend');       %排序
top10_L=I(1:10);oth_L=I(11:end);    %选择前10和归入其他
top10_L=sort(top10_L);              %取消按多少排序，改为按物种名称排序
top10_n=ta(top10_L);                %选择名称
oth_n=ta(oth_L);                    %选择其他

top10_a=spe(top10_L,:);             %选择前十数据
oth_a=sum(spe(oth_L,:));            %归入其他

num=[top10_a;oth_a]; nam=top10_n; nam(11)={'other'};
%%
figure
h=bar((num./s_spe)','stack','EdgeColor',"none");
%% 定义颜色
h(11).FaceColor=("#a6cee3");
h(10).FaceColor=("#2d82af");
h(9).FaceColor=("#98d277");
h(8).FaceColor=("#6f9e4c");
h(7).FaceColor=("#f16667");
h(6).FaceColor=("#f06c45");
h(5).FaceColor=("#fe982c");
h(4).FaceColor=("#d9a295");
h(3).FaceColor=("#7d54a5");
h(2).FaceColor=("#f0eb99");
h(1).FaceColor=("#b15928");

ax = gca;
ax.XTick=1:(length(T1)-1);
ax.XTickLabels = T1(2:end);
%% 根据样品数量X轴坐标方向
switch ceil(length(T1)/6)
    case 1
        ax.XTickLabelRotation = 0;
    case 2
        ax.XTickLabelRotation = 45;
    case 3
        ax.XTickLabelRotation = 45;
    otherwise
        ax.XTickLabelRotation = 90;
end
%% 在坐标轴外绘制
lgd=legend(nam,'Location','eastoutside');
title(lgd,'My Legend Title')
 legend('boxoff')
 