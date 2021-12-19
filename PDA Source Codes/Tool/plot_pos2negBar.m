function plot_pos2negBar(stat)

index=1:10;

b=bar(index,stat(1:2,:)',1);    % 1----width
b(1).FaceColor=[0.5451 0 0];    
b(2).FaceColor=[0 0 0.5451];    
axis([0 11 0 100]);   
h=legend('Positive anomaly','Negative anomaly');       
set(h,'Box','off');
set(gca, 'XTicklabel',{'11^{th}','12^{th}','13^{th}','14^{th}','15^{th}'...
    '16^{th}','17^{th}','18^{th}','19^{th}','20^{th}'});   
xlabel('Century'); 
ylabel('Spatial percentage (%)');       
set(gca,'fontname','Arial','fontsize',9,'linewidth',1,'Tickdir','out');
end

