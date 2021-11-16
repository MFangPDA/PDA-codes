% ������ƽֵ�ٷֱȵ���״ͼ
function plot_pos2negBar(stat)

% data=load('./Result/AnaResult/FieldStats.mat');
% stat=data.FieldStats;

index=1:10;

b=bar(index,stat(1:2,:)',1);  % 1��width
b(1).FaceColor=[0.5451 0 0];    % ���ò�ͬ������ɫ
b(2).FaceColor=[0 0 0.5451];    % ���ò�ͬ������ɫ
axis([0 11 0 100]);   % �������귶Χ
h=legend('Positive anomaly','Negative anomaly');       % ����ͼ��
set(h,'Box','off');
set(gca, 'XTicklabel',{'11^{th}','12^{th}','13^{th}','14^{th}','15^{th}'...
    '16^{th}','17^{th}','18^{th}','19^{th}','20^{th}'});   % X��ļǺ�
xlabel('Century'); % ����X�����ǩ
ylabel('Spatial percentage (%)');       % ����Y�����ǩ
set(gca,'fontname','Arial','fontsize',9,'linewidth',1,'Tickdir','out');%ͳһ����matlabͼ�����壬��С

end

