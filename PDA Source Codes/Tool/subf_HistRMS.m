function subf_HistRMS(RMS_txt)
RMS_Matrix=RMS_txt;
RMS=subf_Matrix2Vector(RMS_Matrix);
X=0:0.05:2;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=0+(i-1)*0.05;
    tail=0+i*0.05;
    for j=1:length(RMS)
        if RMS(j)>=head && RMS(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,'linewidth',1);
set(gca,'XLim',[0 2]); set(gca,'XTick',0:0.5:2);
set(gca,'YLim',[0 500]); set(gca,'YTick',0:100:500);
xlabel('Root mean square error'); ylabel('Frequency');
text(0.1,450,strcat('Mean: ',num2str(mean(RMS),'%.2f')),'color','k','fontsize',20,'fontweight','bold'); 
% set(gcf,'color','none');   % 图形背景设为无色
set(gca,'color','none');     % 坐标轴背景为无色，这条很重要，通常图形背景的白色（实际为坐标轴默认背景色），设置为透明，就可以进行多幅图的叠加
set(gca,'linewidth',1);
set(gca,'fontname','cambria','fontsize',10);
box off;
hold on;
end
