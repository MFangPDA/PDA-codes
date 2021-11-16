function postf2_plotVerification()

% 读取伪代用资料实验的结果
data1=load('./Result/VerResult/corMPIvct.mat'); corMPIvct=data1.corMPIvct;
data2=load('./Result/VerResult/corPDAvct.mat'); corPDAvct=data2.corPDAvct;
data3=load('./Result/VerResult/ceMPIvct.mat');  ceMPIvct=data3.ceMPIvct;
data4=load('./Result/VerResult/cePDAvct.mat');  cePDAvct=data4.cePDAvct;

% 读取空间验证结果
data5=load('./Result/VerResult/corPDAFieldvct.mat');  corPDAFieldvct=data5.noNaN1;
data6=load('./Result/VerResult/cor20CRFieldvct.mat'); cor20CRFieldvct=data6.noNaN2;
data7=load('./Result/VerResult/corMPIFieldvct.mat');  corMPIFieldvct=data7.noNaN3;

data8=load('./Result/VerResult/cePDAFieldvct.mat');  cePDAFieldvct=data8.noNaN1;
data9=load('./Result/VerResult/ce20CRFieldvct.mat'); ce20CRFieldvct=data9.noNaN2;
data10=load('./Result/VerResult/ceMPIFieldvct.mat'); ceMPIFieldvct=data10.noNaN3;

% 作图，直方图
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
subplot(1,4,3);
plotHistCOR(corMPIvct,'Correlation coefficient','b'); hold on
plotHistCOR(corPDAvct,'Correlation coefficient','r'); hold on
subplot(1,4,4);
plotHistCE(ceMPIvct,'Coefficient of efficiency','b'); hold on
plotHistCE(cePDAvct,'Coefficient of efficiency','r'); hold on

subplot(1,4,1);
plotHistCOR(corPDAFieldvct,'Correlation coefficient','r');  hold on
plotHistCOR(cor20CRFieldvct,'Correlation coefficient','b'); hold on
plotHistCOR(corMPIFieldvct,'Correlation coefficient','g');  hold on
subplot(1,4,2);
plotHistCE(cePDAFieldvct,'Coefficient of coefficiency','r');  hold on
plotHistCE(ce20CRFieldvct,'Coefficient of coefficiency','b'); hold on
plotHistCE(ceMPIFieldvct,'Coefficient of coefficiency','g');  hold on
end


function plotHistCOR(vct,paraName,colorName)
X=-1:0.05:1;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=-1+(i-1)*0.05;
    tail=-1+i*0.05;
    for j=1:length(vct)
        if vct(j)>=head && vct(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,colorName,'linewidth',1.5);
set(gca,'XLim',[-1 1]); set(gca,'XTick',-1:0.5:1);
set(gca,'YLim',[0 100]); set(gca,'YTick',0:20:100);
xlabel(paraName); ylabel('Frequency');
text(-0.9,80,strcat('Mean: ',num2str(mean(vct),'%.2f')),'color',colorName,'fontsize',14,'fontweight','bold'); 
% set(gcf,'color','none');   % 图形背景设为无色
set(gca,'color','none');     % 坐标轴背景为无色，这条很重要，通常图形背景的白色（实际为坐标轴默认背景色），设置为透明，就可以进行多幅图的叠加
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end

function plotHistCE(vct,paraName,colorName)
X=-5:0.15:1;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=-5+(i-1)*0.15;
    tail=-5+i*0.15;
    for j=1:length(vct)
        if vct(j)>=head && vct(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,colorName,'linewidth',1.5);
set(gca,'XLim',[-5 1]); set(gca,'XTick',-5:1:1);
set(gca,'YLim',[0 200]); set(gca,'YTick',0:40:200);
xlabel(paraName); ylabel('Frequency');
text(-4.5,160,strcat('Mean: ',num2str(mean(vct),'%.2f')),'color',colorName,'fontsize',12,'fontweight','bold'); 
% set(gcf,'color','none');   % 图形背景设为无色
set(gca,'color','none');     % 坐标轴背景为无色，这条很重要，通常图形背景的白色（实际为坐标轴默认背景色），设置为透明，就可以进行多幅图的叠加
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end
