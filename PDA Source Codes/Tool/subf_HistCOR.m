function subf_HistCOR(Cor_txt)
Cor_Matrix=Cor_txt;
Cor=subf_Matrix2Vector(Cor_Matrix);
X=-1:0.05:1;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=-1+(i-1)*0.05;
    tail=-1+i*0.05;
    for j=1:length(Cor)
        if Cor(j)>=head && Cor(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,'linewidth',1);
set(gca,'XLim',[-1 1]); set(gca,'XTick',-1:0.5:1);
set(gca,'YLim',[0 500]); set(gca,'YTick',0:100:500);
xlabel('Correlation coefficient'); ylabel('Frequency');
text(-0.9,450,strcat('Mean: ',num2str(mean(Cor),'%.2f')),'color','k','fontsize',20,'fontweight','bold'); 
% set(gcf,'color','none');   % ͼ�α�����Ϊ��ɫ
set(gca,'color','none');     % �����ᱳ��Ϊ��ɫ����������Ҫ��ͨ��ͼ�α����İ�ɫ��ʵ��Ϊ������Ĭ�ϱ���ɫ��������Ϊ͸�����Ϳ��Խ��ж��ͼ�ĵ���
set(gca,'linewidth',1);
set(gca,'fontname','cambria','fontsize',10);
box off;
hold on;
end


