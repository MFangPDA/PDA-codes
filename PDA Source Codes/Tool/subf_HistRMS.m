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
% set(gcf,'color','none');   
set(gca,'color','none');     
set(gca,'linewidth',1);
set(gca,'fontname','cambria','fontsize',10);
box off;
hold on;
end
