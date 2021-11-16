function postf5_plotAAForcing()

addpath(genpath(pwd));

[NUM,TXT]=xlsread('./Alldata.xlsx','alldata_original');

figure
GHG=NUM(16:986,2); AMO=NUM(16:986,7); AAI=NUM(16:986,10);
data=[GHG,AMO,AAI];

maxAAI=max(data(:,3));      minAAI=min(data(:,3));      meanAAI=mean(data(:,3));
MCAAAI=mean(data(1:86,3));  pre=mean(data(700:835,3));  ind=mean(data(836:end,3));

% 做31年平滑
data_smooth=zeros(size(data));
for i=1:3
    data_smooth(:,i)=smooth(data(:,i),31,'moving');
end

% 去除线性趋势(Remove linear trends)
data_detrend=zeros(size(data));
data_smooth_detrend=zeros(size(data_smooth));
for i=1:3
    data_detrend(:,i)=detrend(data(:,i));
    data_smooth_detrend(:,i)=detrend(data_smooth(:,i));
end

% 分成两个时间段
% preindustrial era  1849-1015+1=835   
[r1,p1]=partialcorr(data_detrend(1:835,:),'type','Pearson');

% industrial era  1849-1015+1=835   
[r2,p2]=partialcorr(data_detrend(836:end,:),'type','Pearson');

parameters=roundn([r1,p1;r2,p2],-3);
disp(parameters);

set(gcf,'unit','inches','position',[0.2,0.2,7.3,5]);
ha = tight_subplot(2,2,[.08 .08],[.08 .08],[.06 .02]);
t=1015:1985;
Y=zeros(length(t),1);
axes(ha(1)); plotForcingBar(t,data(:,2),Y); 
axes(ha(2)); plot(t,data(:,3),'k-','linewidth',1);  set(gca,'ylim',[1 2.5]);
str1=strcat('Preindustrial era:', '{\itr}','=',32,num2str(parameters(3,1)),', p < 0.001');
str2=strcat('Industrial era:', '{\itr}','=',32,num2str(parameters(6,1)),', p < 0.001');
text(1015,2,{str1;str2});

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', []);
set(gca,'color','none');
set(gca,'YColor','k');     % Y轴设为蓝色
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',9,'fontname','Arial');

axes(ha(3)); plotForcingBar(t,data(:,1),Y); 
axes(ha(4)); plot(t,data(:,3),'k-','linewidth',1);  set(gca,'ylim',[1 2.5]);
str3=strcat('Preindustrial era: ', '{\itr}','=',32,num2str(parameters(3,2)),', {\itp} < 0.001');
str4=strcat('Industrial era: ', '{\itr}','=',32,num2str(parameters(6,2)),', {\itp} < 0.001');
text(1015,2,{str3;str4});

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', []);
set(gca,'color','none');
set(gca,'YColor','k');     % Y轴设为蓝色
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',7,'fontname','Arial'); 

% 重新作图，把AMO和GHG做30年滑动
% 做31年平滑
figure
set(gcf,'unit','inches','position',[0.2,0.2,7.3,5]);
ha = tight_subplot(2,4,[.08 .08],[.08 .08],[.08 .08]);
t=1015:1985;
Y=zeros(length(t),1);
axes(ha(1)); plotForcingBar(t,data_smooth(:,2),Y); 
axes(ha(2)); plot(t,data(:,3),'k-','linewidth',1);   set(gca,'ylim',[1 2.5],'fontsize',7,'fontname','Arial'); 
str1=strcat('Preindustrial era:', '{\itr}','=',32,', {\itp} = ');
str2=strcat('Industrial era:', '{\itr}','=',32,', {\itp} = ');
text(1015,2,{str1;str2},'fontsize',7,'fontname','Arial');

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', [],'fontsize',7,'fontname','Arial');
set(gca,'color','none');
set(gca,'YColor','k');     % Y轴设为黑
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',7,'fontname','Arial');

axes(ha(5)); plotForcingBar(t,data_smooth(:,1),Y); 
axes(ha(6)); plot(t,data(:,3),'k-','linewidth',1);   set(gca,'ylim',[1 2.5],'fontsize',7,'fontname','Arial');
str3=strcat('Preindustrial era: ', '{\itr}','=',32,', p = ');
str4=strcat('Industrial era: ', '{\itr}','=',32,', p = ');
text(1015,2,{str3;str4});

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', [],'fontsize',7,'fontname','Arial');
set(gca,'color','none');
set(gca,'YColor','k');     % Y轴设为蓝色
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',7,'fontname','Arial');

% 回归方程
% AMO only:
% AAI = 1.79361314779342 + 0.66728626*AMO
predAAI=1.79361314779342 + 0.66728626*data_smooth(:,2);
YY=ones(length(t),1);
axes(ha(3)); plot(t,predAAI,'b-','linewidth',1);   set(gca,'ylim',[1 2.5],'YColor','b','fontsize',7,'fontname','Arial'); hold on
ylabel({'Predicted AA index'});
set(gca,'color','none','YColor','b','tickdir','out','linewidth',1,'box','off');

axes(ha(4)); plot(t,data(:,3),'k-','linewidth',1); set(gca,'ylim',[1 2.5]); 
str3=strcat('Preindustrial era: ',32, '{\itr}','=',32,', {\itp} = ');
str4=strcat('Industrial era: ',32, '{\itr}','=',32,', {\itp} = ');
text(1015,2,{str3;str4},'fontsize',7,'fontname','Arial');

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', []);
set(gca,'color','none','YColor','k','tickdir','out','linewidth',1,'box','off');
set(gca,'fontsize',7,'fontname','Arial');

% AMO and GHG:
% AAI = 1.8462233474151262 + 0.75945257*AMO - 0.32731893*GHG
predAAI=1.8462233474151262 + 0.75945257*data_smooth(:,2) - 0.32731893*data_smooth(:,1);
axes(ha(7)); plot(t,predAAI,'b-','linewidth',1);   set(gca,'ylim',[1 2.5],'YColor','b','fontsize',7,'fontname','Arial');  hold on
ylabel({'Predicted AA index'});
set(gca,'color','none','YColor','k','tickdir','out','linewidth',1,'box','off','fontsize',7,'fontname','Arial');

axes(ha(8)); plot(t,data(:,3),'b-','linewidth',1); set(gca,'ylim',[1 2.5]);  
str3=strcat('Preindustrial era: ',32, '{\itr}','=',32,', {\itp} = ');
str4=strcat('Industrial era: ',32, '{\itr}','=',32,', {\itp} = ');
text(1015,2,{str3;str4},'fontsize',7,'fontname','Arial');

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', [],'fontsize',7,'fontname','Arial');
set(gca,'color','none','YColor','k','tickdir','out','linewidth',1,'box','off');
set(gca,'fontsize',7,'fontname','Arial');

% 补充PDO的数据作为补充材料
PDO=NUM(16:986,8);
PDO_smooth=smooth(PDO,31,'moving');
figure
set(gcf,'unit','inches','position',[0.2,0.2,7.3,5]);
ha = tight_subplot(1,2,[.08 .08],[.08 .08],[.06 .02]);
axes(ha(1)); plotForcingBar(t,PDO_smooth,Y); 
axes(ha(2)); plot(t,data(:,3),'k-','linewidth',1);   set(gca,'ylim',[1 2.5]); 
str1=strcat('Preindustrial era:',32, '{\itr}','=',32,', {\itp} = ');
str2=strcat('Industrial era:',32, '{\itr}','=',32,', {\itp} = ');
text(1015,2,{str1;str2});

ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); set(gca,'XTick', []);
set(gca,'color','none');
set(gca,'YColor','k');     % Y轴设为黑
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',7,'fontname','Arial');

[R,P]=corrcoef([PDO_smooth, data(:,3)]);
disp([R,P]);
end


% 以下函数默认左侧为红色，传入的参数为右侧坐标颜色、标题
function myPlotyy(x,y1,y2,yLim,yTick,yTicklabel,r1,p1,r2,p2)
% 获取坐标轴、图像句柄
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(AX,'FontSize',12,'FontName','cambria'); % 设置x轴、左y轴、右y轴刻度字号和字型

% 左轴
set(AX(1),'XColor','k','YColor',[0.5451 0 0],'fontsize',9,'fontname','Arial','Linewidth',1,'Tickdir','out');
% 右轴
set(AX(2),'XColor','k','YColor',[0 0 0.5451],'fontsize',9,'fontname','Arial','Linewidth',1,'Tickdir','out');

HH1=get(AX(1),'Ylabel');
set(HH1,'color',[0.5451 0 0]);     % set(HH1,'String','Left Y-axis');
set(H1,'LineStyle','-');  set(H1,'Linewidth',1.5); set(H1,'color',[0.5451 0 0]);

HH2=get(AX(2),'Ylabel');
set(HH2,'color',[0 0 0.5451]);     % set(HH2,'String','Right Y-axis');
set(H2,'LineStyle','-');  set(H2,'Linewidth',1.5); set(H2,'color',[0 0 0.5451]);

set(AX(1),'ylim',[1 2.5],'yTick',[1:0.5:2.5],'YTickLabel',{'1','1.5','2','2.5'});
set(AX(2),'ylim',yLim,'yTick',yTick,'YTickLabel',yTicklabel);

text(1000,2.7,'A','fontsize',9,'fontweight','bold');
text(1060,2.7,'\color[rgb]{0.5451 0 0}{AA index} \color{black}{\itvs} \color[rgb]{0 0 0.5451}{Model}','fontsize',9);

text(1020,1.32,strcat('\itr ','\rm =',32,num2str(r1),', p-value =',32,num2str(p1)),'fontsize',8,'Color',[0.5451 0 0.5451]);
text(1020,1.12,strcat('\itr ','\rm =',32,num2str(r2),', p-value =',32,num2str(p2)),'fontsize',8,'Color',[0 0.5451 0.5451]);
end

function [Slope, Pvalue, trendValue]=trendAndSlope(data,len)

[Slope,Pvalue]=MKTrendDetection(data,'');

% 对数据进行线性拟合
tt=1:length(data);
[p,S]=polyfit(tt',data,1);  % 必须是列向量
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),len);
end

function plotForcingBar(t,forcing,Y)
    pos=forcing; pos(pos<0)=0;     neg=forcing; neg(pos>0)=0;
    fill([t,fliplr(t)],[pos',fliplr(Y')],[0.5451 0 0],'edgealpha',0); hold on
    fill([t,fliplr(t)],[Y',fliplr(neg')],[0 0 0.5451],'edgealpha',0); hold on
    ylabel('ylabel','fontsize',7);  xlabel('Year (AD)','fontsize',7); title('A','fontweight','bold');
    set(gca,'Tickdir','out','fontname','Arial','fontsize',7,'Linewidth',1);
    set(gca,'box','off');
end

function plotAAIBar(t,AAI,YY)
    fill([t,fliplr(t)],[AAI',fliplr(YY')],[0.5451 0 0.5451],'edgealpha',0); hold on
    set(gca,'ylim',[1 2.5]);
    % plot(t,AAI,'color','b','linewidth',1.5); hold on
    ylabel('AA index'); xlabel('Year (AD)');  title('A','fontweight','bold');
    set(gca,'YAxisLocation','left','box','on');
    % set(gca,'xTick',[],'xTicklabel',[]);
    set(gca,'color','none','Ycolor','k','Tickdir','out','fontname','Arial','fontsize',9,'Linewidth',1);
end

function plotScatter(forcing,AAI)
    s=scatter(forcing,AAI);
    s.LineWidth = 0.2;
    s.MarkerEdgeColor = 'none';
    s.MarkerFaceColor = 'b';
    set(gca,'box','on');
    xlabel('xlabel'); ylabel('AA index');  title('A','fontweight','bold');
    set(gca,'color','none','Tickdir','out','fontname','Arial','fontsize',9,'Linewidth',1);
    hold on
    
    % 拟合
    [p2,S]=polyfit(forcing,AAI,1);  
    predAAI=p2(1)*forcing+p2(2);
    plot(forcing,predAAI,'k-');
        
    % R^2=1-(S.normr/norm(y-mean(y)))^2
    R2=1-(S.normr/norm(AAI-mean(AAI)))^2;  
    
    disp(p2); disp(R2);
    % text(1,2,strcat('y = ',num2str(roundn(p2(1),-3)),'x + ',num2str(roundn(p2(2),-3))),'fontsize',8,'Color','k');
    % text(1,2,strcat('R^2 = ',num2str(roundn(R2,-3))),'fontsize',8,'Color','k');
end

% Text字符的颜色设置有下面两种方法。
%
% ①\color{颜色名}颜色名：颜色名有12种，分别为red、green、yellow、magenta、blue、black、white、 cyan、gray、barkGreen、orange和lightBlue。例如：\color{magenta}magenta。
%
% ②\color[rgb]{a b c}:设置字体颜色为RGB矩阵[a b c]所表示的颜色。 a、b和c都在[0 1] 范围内。例如：color[rgb]{0 .5 .5}。