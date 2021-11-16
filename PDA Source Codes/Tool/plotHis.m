function plotHis()
data=load('vctAnomalys.mat');
vctValues=data.AnomalyV;
Cor=vctValues(:,1);
X=-1:0.05:0;
N=length(X); 
NegCount=zeros(N,1);
for i=1:N-1
    head=-1+(i-1)*0.05;
    tail=-1+i*0.05;
    for j=1:length(Cor)
        if Cor(j)>head && Cor(j)<=tail
            NegCount(i)=NegCount(i)+1;
        end
    end
end
A=NegCount;


X=0:0.05:1;
N=length(X); 
PosCount=zeros(N,1);
for i=1:N-1
    head=0+(i-1)*0.05;
    tail=0+i*0.05;
    for j=1:length(Cor)
        if Cor(j)>head && Cor(j)<=tail
            PosCount(i)=PosCount(i)+1;
        end
    end
end
B=PosCount;

range=-1:0.05:1;
sumN=length(range);
AA=zeros(sumN,1);  AA(1:21)=A;
BB=zeros(sumN,1);  BB(21:41)=B;

histogram(Cor(Cor<0),-1:0.05:1,'edgealpha',0,'facealpha',1,'facecolor','b');   % -1:0.05:1 X�᷶Χ
hold on
histogram(Cor(Cor>=0),-1:0.05:1,'edgealpha',0,'facealpha',1,'facecolor','r');
hold on
% stairs(X,Count,'linewidth',0.5);
% set(gca,'XLim',[-2 2]); set(gca,'XTick',-2:1:2);
% set(gca,'YLim',[0 2000]); set(gca,'YTick',0:500:2000);
% if (flag>=6)
%     xlabel('Anoamly values');
% end
% if (flag==1 || flag==6)
%     ylabel('Frequency');
% end
% % set(gcf,'color','none');   % ͼ�α�����Ϊ��ɫ
% % set(gca,'color','none');   % �����ᱳ��Ϊ��ɫ����������Ҫ��ͨ��ͼ�α����İ�ɫ��ʵ��Ϊ������Ĭ�ϱ���ɫ��������Ϊ͸�����Ϳ��Խ��ж��ͼ�ĵ���
% % set(gca,'linewidth',0.5,'color','k');
% set(gca,'fontsize',10,'fontname','cambria');
% box off;
hold on;
end

