% ****************************************************************************
% Instruction£º
%     To show the PDA-based reconstruction. All of the time series are anomalous 
%     relative to their individual millennial means
% ****************************************************************************
function postf3_showResult()
%
% To load subfunction in tool
addpath(genpath(pwd));
%%
%
% Global parameters
global startDA   endDA   lenDA;
global smoothSpan;

startDA=1000; endDA=2000;
lenDA=endDA-startDA+1;
smoothSpan=11;
t=startDA:endDA;
%
%
% To load PDA-based reconstruction
PDA=load('./Result/PDAResult/PDATas.mat');
PDATas=PDA.ANNTas;

% All of the time series are anomalous relative to their individual millennial means
Anomaly=subf_yrAnomaly(PDATas,startDA,endDA,startDA,endDA);

% 11-yrs smoothing with lowess method
Smooth=zeros(size(Anomaly));  A=size(Anomaly);
for i=1:A(1,1)
    for j=1:A(1,2)
        Smooth(i,j,:)=smooth(Anomaly(i,j,:),smoothSpan,'lowess');
    end
end
regSeries=zeros(1,lenDA);
for y=1:lenDA    % y: year
    regSeries(y)=mean2(Smooth(:,:,y));
end

%%
%
% 1. Century-by-century tas anomaly field
figure
set(gcf,'unit','inches','position',[0.2,0.2,7.3,3]);
ha = tight_subplot(2,5,[.02 .02],[.14 .02],[.08 .08]);
PosiPerc=zeros(1,10); NegaPerc=zeros(1,10);
MaxV=zeros(1,10);     MinV=zeros(1,10);
AnomalyVct=zeros(45*180,10);
for i=1:10
    axes(ha(i));
    head=1000+(i-1)*100+1; tail=1000+i*100;
    centruyMean=mean(Smooth(:,:,(head-startDA+1):(tail-startDA+1)),3);
    plot_Climatology(centruyMean,['\fontsize{9}',num2str(head),'-',num2str(tail),' CE']);
    hold on
    
    AnomalyVct(:,i)=reshape(centruyMean,45*180,1);
    
    D=size(centruyMean);
    A=find(centruyMean>0);  
    PosiPerc(i)=roundn((length(A)/(D(1,1)*D(1,2)))*100,-2);  NegaPerc(i)=100-PosiPerc(i);
    
    MaxV(i)=max(max(centruyMean));          MinV(i)=min(min(centruyMean));
end
fieldStats=[PosiPerc;NegaPerc;MaxV;MinV]; 

%%
% 2. Zonal variability of tas anomaly
% 2.1 raster figure of zonal variability
figure
set(gcf,'unit','inches','position',[0.2,0.2,7.087,3]);
ha = tight_subplot(1,2,[.1 .01],[.15 .15],[.05 .05]);

N=45;
zonalSeries=zeros(N,lenDA);  % with 2¡ãresolution, 45 zonal zones
for i=1:lenDA
    for j=1:N
        zonalSeries(j,i)=mean2(Smooth(j,:,i));
    end
end

% % 30-yrs mean
% zonalSeries30yrs=zeros(N,lenDA);
% % 1000-1400
% for k=1:N
%     zonalSeries30yrs(k,1:41)=mean(zonalSeries(k,1:41));
% end
% for i=1:32
%     head=(i-1)*30+1+41; tail=i*30+41;
%     for k=1:N
%         zonalSeries30yrs(k,head:tail)=mean(zonalSeries(k,head:tail));
%     end
% end

% Contour map of zonal variability of the NH tas anomaly
axes(ha(1)); imagesc(zonalSeries);  box off
xlabel('Year (AD)'); ylabel({'Latitude (\circN)'}); title('Zonal pattern of anomaly');
set(gca,'fontsize',7,'fontname','Arial','Tickdir','out','linewidth',1);
set(gca,'YLim',[1,45],'YTick',[1:5:45],'YTickLabel',...
    {'90','80','70','60','50','40','30','20','10','0'});
set(gca,'XLim',[0,1000],'XTick',[0:200:1000],'XTickLabel',{'1000','1200','1400','1600','1800','2000'});
set(gca,'fontsize',7,'fontname','Arial');

myColor=genColor();
colormap(myColor);
hcb=colorbar('Ticks',-1:0.5:1,'location','Westoutside','fontname','Arial','fontsize',7);
caxis([-1,1]);
ylabel(hcb,'Temperature anomaly (\circC)','Fontsize',7,'FontName','Arial');

% 2.2 Overlap the AA index onto the Contour map
[NUM,TXT]=xlsread('./Alldata.xlsx','alldata_original');
axes(ha(2));
plot(NUM(:,1),NUM(:,end),'k-','linewidth',1);  hold on;
ylabel({'AA index'});
set(gca,'YAxisLocation','right','XAxisLocation','top'); 
set(gca,'XTick', []);
set(gca,'color','none');
set(gca,'YColor','k');       
set(gca,'tickdir','out');
set(gca,'linewidth',1);
set(gca,'box','off');
set(gca,'fontsize',9,'fontname','Arial');
end

% Customed color
function myColor=genColor()
myColor=zeros(64,3);
% matlab colourmap is a marix with format of 64 rows and 3 colums
flag=0.5451;  % DarkBlue:[0 0 0.5451]   DarkRed:[0.5451 0 0]
% Darkblue
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=linspace(flag,1,32);
% Darkred
myColor(33:64,1)=flip(linspace(flag,1,32)); myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32));
end
