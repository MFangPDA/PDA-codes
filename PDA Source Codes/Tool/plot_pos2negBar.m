% 正负距平值百分比的条状图
function plot_pos2negBar(stat)

% data=load('./Result/AnaResult/FieldStats.mat');
% stat=data.FieldStats;

index=1:10;

b=bar(index,stat(1:2,:)',1);  % 1是width
b(1).FaceColor=[0.5451 0 0];    % 设置不同线条颜色
b(2).FaceColor=[0 0 0.5451];    % 设置不同线条颜色
axis([0 11 0 100]);   % 控制坐标范围
h=legend('Positive anomaly','Negative anomaly');       % 设置图例
set(h,'Box','off');
set(gca, 'XTicklabel',{'11^{th}','12^{th}','13^{th}','14^{th}','15^{th}'...
    '16^{th}','17^{th}','18^{th}','19^{th}','20^{th}'});   % X轴的记号
xlabel('Century'); % 设置X坐标标签
ylabel('Spatial percentage (%)');       % 设置Y坐标标签
set(gca,'fontname','Arial','fontsize',9,'linewidth',1,'Tickdir','out');%统一设置matlab图的字体，大小

end

