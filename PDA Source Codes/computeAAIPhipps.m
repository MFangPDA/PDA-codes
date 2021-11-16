% ****************************************************************************
% Description:
%     (1) Computing the AA index values based on different simulations from
%         CSIRO-MK-3L 1500yrs simulation with different forcings
%     (2) All of the time series are anomalous relative to their individual millennial means
%     Note that: the choice of the reference time span do not affect the AAI computation
% ****************************************************************************
function computeAAIPhipps()

% To call sub functions in subdirectory
addpath(genpath(pwd));

% Phipps et al., 2014, Journal of Climate
% 0001-2000 AD monthly

% 1) To load CSIRO simulation
% Three realizations with different initial conditions
% To compute the AA during 1000-2000

global AAIstart
AAIstart=1015; AAIend=1985;
lenAAI=AAIend-AAIstart+1;

yearAAI_O=zeros(lenAAI,3);      AAtrend_O=zeros(lenAAI,3);
yearAAI_OG=zeros(lenAAI,3);     AAtrend_OG=zeros(lenAAI,3);

% 1��������������AA index
for i=1:3
    % O
    [NHseries_O,ARCseries_O,tasField_O]=readPhipps(strcat('o',num2str(i)));
    % OG
    [NHseries_OG,ARCseries_OG,tasField_OG]=readPhipps(strcat('og',num2str(i)));
    
    % 2) To compute AA index
    [AAyearValue_O,trendValue_O]=computeAAI(NHseries_O((AAIstart-15):(AAIend+15)),ARCseries_O((AAIstart-15):(AAIend+15)));
    yearAAI_O(:,i)=AAyearValue_O;  AAtrend_O(:,i)=trendValue_O;
    
    [AAyearValue_OG,trendValue_OG]=computeAAI(NHseries_OG((AAIstart-15):(AAIend+15)),ARCseries_OG((AAIstart-15):(AAIend+15)));
    yearAAI_OG(:,i)=AAyearValue_OG;  AAtrend_OG(:,i)=trendValue_OG;
end

% 2�����㼯��ƽ�������
AAyrs_O_Mean=mean(yearAAI_O,2);  
AAyrs_OG_Mean=mean(yearAAI_OG,2);
% �������±߽���
data1=minmax(yearAAI_O);  AA_O_lower=data1(:,1);  AA_O_upper=data1(:,2); 
data2=minmax(yearAAI_OG); AA_OG_lower=data2(:,1); AA_OG_upper=data2(:,2); 

% 3������ƽ��ֵ����Ϊ��ҵ������ǰ�͹�ҵ����ʱ��
AA_O=mean(AAyrs_O_Mean);
AA_O_mean_preindus=mean(AAyrs_O_Mean(1:(1849-AAIstart+1)));    AA_O_mean_indus=mean(AAyrs_O_Mean((1850-AAIstart+1):end));

AA_OG=mean(AAyrs_OG_Mean);
AA_OG_mean_preindus=mean(AAyrs_OG_Mean(1:(1849-AAIstart+1)));  AA_OG_mean_indus=mean(AAyrs_OG_Mean((1850-AAIstart+1):end));
disp([AA_O AA_OG; AA_O_mean_preindus AA_O_mean_indus; AA_OG_mean_preindus AA_OG_mean_indus]);

% 4������1000���AAI����
figure
set(gcf,'unit','inches','position',[0.2,0.2,7.087,5]);
ha = tight_subplot(2,1,[.08 .06],[.08 .08],[.08 .08]);

% 4.1 Industrial era
axes(ha(1));
t=1850:AAIend;

fill([t,fliplr(t)],[AA_O_upper((1850-AAIstart+1):end)',fliplr(AA_O_lower((1850-AAIstart+1):end)')],[0.5451 0 0],'edgealpha',0,'facealpha',0.3); hold on
fill([t,fliplr(t)],[AA_OG_upper((1850-AAIstart+1):end)',fliplr(AA_OG_lower((1850-AAIstart+1):end)')],[0 0 0.5451],'edgealpha',0,'facealpha',0.3); hold on

AAyrs_O_ind=AAyrs_O_Mean((1850-AAIstart+1):end);
AAyrs_OG_ind=AAyrs_OG_Mean((1850-AAIstart+1):end);

h1=plot(t,AAyrs_O_ind,'color',[0.5451 0 0],'linewidth',1.5); hold on
h2=plot(t,AAyrs_OG_ind,'color',[0 0 0.5451],'linewidth',1.5); hold on
h3=legend([h1,h2],'Without GHGs forcing','With GHGs forcing','Location','NorthEast'); set(h3,'Box','off');

name1=strcat('Mean =',32,' '); text(1900,2.8,name1,'color',[0.5451 0 0],'fontsize',7,'fontname','Arial');
name2=strcat('Mean =',32,' '); text(1900,2.2,name2,'color',[0 0 0.5451],'fontsize',7,'fontname','Arial');

xlabel('Year (AD)');  ylabel('AA index');
set(gca,'Linewidth',1,'fontname','Arial','fontsize',7,'ylim',[1.5,3.5],'Tickdir','out');
set(gca,'XLim',[1850 2000],'XTick',1850:10:2000);

% 4.2 The millennium
axes(ha(2));
t=AAIstart:1985;

fill([t,fliplr(t)],[AA_O_upper',fliplr(AA_O_lower')],[0.5451 0 0],'edgealpha',0,'facealpha',0.3); hold on
fill([t,fliplr(t)],[AA_OG_upper',fliplr(AA_OG_lower')],[0 0 0.5451],'edgealpha',0,'facealpha',0.3); hold on

h1=plot(t,AAyrs_O_Mean,'color',[0.5451 0 0],'linewidth',1.5); hold on
h2=plot(t,AAyrs_OG_Mean,'color',[0 0 0.5451],'linewidth',1.5); hold on
h3=legend([h1,h2],'Without GHGs forcing','With GHGs forcing','Location','NorthEast'); set(h3,'Box','off');

name1=strcat('Mean =',32,''); text(1500,2.8,name1,'color',[0.5451 0 0],'fontsize',7,'fontname','Arial');
name2=strcat('Mean =',32,''); text(1500,2.2,name2,'color',[0 0 0.5451],'fontsize',7,'fontname','Arial');

xlabel('Year (AD)');  ylabel('AA index');
set(gca,'Linewidth',1,'fontname','Arial','fontsize',7,'ylim',[1.5,3.5],'Tickdir','out');

% text(1866,3.7,'AA index \color{blue}with GHG \color{black}{\itvs} \color{red}without GHG','fontsize',9);
% �����Լ���
% (1) ANOVA���鷽��
% ��ҵ����֮ǰ
% Y=[AAyrs_O_Mean(1:(1850-AAIstart+1)) AAyrs_OG_Mean(1:(1850-AAIstart+1))];
% group={'O','OG'};
% [p,tbl,stats] = anova1(Y,group,'on');
% % ��ҵ����ʱ��
% Y=[AAyrs_O_Mean((1850-AAIstart+1):end) AAyrs_OG_Mean((1850-AAIstart+1):end)];
% group={'O','OG'};
% [p,tbl,stats] = anova1(Y,group,'on');

% ��2��t���鷽��
% ��ҵ����֮ǰ
[h,p,ci,stats] = ttest2(AAyrs_O_Mean(1:(1850-AAIstart+1)),AAyrs_OG_Mean(1:(1850-AAIstart+1)),'Vartype','unequal')
% ��ҵ����ʱ��
[h,p,ci,stats] = ttest2(AAyrs_O_Mean((1850-AAIstart+1):end),AAyrs_OG_Mean((1850-AAIstart+1):end),'Vartype','unequal')

% h=0 ����������5%�����Ŷ��±����ܣ���x,y��ͳ���Ͽ��Կ�������ͬһ�ֲ������ݣ�
% h=1 �������豻�ܾ�����x,y��ͳ���ϱ���Ϊ�������ڲ�ͬ�ֲ������ݣ��������ֶȣ�
end

% To load the simulation of CSIROMk3L-1-2, time span is 000101-200012 AD
function [NHseries,Arcseries,tasField]=readPhipps(name)

tas_Amon=subf_readNC(strcat('./tsc_',name,'.nc'),'tsc');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% ����⣬CSIROMk3L-1-2���ݽṹ��PMIP3������һ��
% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end

field=geoFormat;
% tasField=field;

tasField=subf_yrAnomaly(field,1,2000,1961,1990);
imagesc(tasField(:,:,end));

% To compute the NH-averaged series and the Arctic-averaged series
A=size(tasField);
NHseries=zeros(A(1,3),1);
Arcseries=zeros(A(1,3),1);
for y=1:A(1,3)    % y: year
    NHseries(y)=mean2(tasField(1:28,:,y));
    Arcseries(y)=mean2(tasField(1:9,:,y));
end
end

% To compute AA index from NH-averaged series and Arctic-averaged series
% ֻ��ȡ1000-2000�������������
function [yrAAI,trendValue]=computeAAI(NHseries,ARCseries)

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-2000,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1985��
global AAIstart
lenAAI=(1985-AAIstart+1);
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+30;      % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail); % Arctic TAS series
    x=NHseries(head:tail);  % NH TAS series
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    yrAAI(i)=b(2);
end
%
% ��AAI���߽����������
tt=1:length(yrAAI);
[p,S]=polyfit(tt',yrAAI,1);  % ������������
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),length(yrAAI));
end

function [Slope,Pvalue]=MKTrendDetection(series,name)
X=series;
N=length(X);
U=zeros(N-1,1);
for t=2:N
    x=X(1:t);
    S=0;
    n=length(x);
    for k=1:(n-1)
        for j=(k+1):n
            S=S+sign(x(j)-x(k));
        end
    end
    VarS=n*(n-1)*(2*n+5)/18;
    if S>0
        Z=(S+1)/sqrt(VarS);
    elseif S==0
        Z=0;
    else
        Z=(S-1)/sqrt(VarS);
    end
    U(t-1)=Z;
end

% significan level
Alpha=1-normcdf(U(end),0,1);

% variation rate
Qi=zeros(N*(N-1)/2,1);
counter=1;
for k=1:(N-1)
    for j=(k+1):N
        Qi(counter)=(X(j)-X(k))/(j-k);
        counter=counter+1;
    end
end
% Sen's Slope:the trend
Q=median(Qi)*30;  % ע�⣺��Phipps�����ݼ������Ƶ�ʱ���������1��1��ֵ

fprintf('%12s \n',strcat('*******',32,'Mann-Kendall is performing on:',32,name,32,'*******'));
fprintf('%12s \n',strcat('>> AA index trend =',num2str(mean(Q),'%.3f'),'/30 yrs'));
fprintf('%12s \n',strcat('>> AA index trend is statistically significant at Alpha=',num2str(1-Alpha),32,'level'));

Slope=Q;
Pvalue=1-Alpha;

% sign(x)��Signum function
% x<0,sign(x)=-1
% x=0,sign(x)=0
% x>0,sign(x)=1
end

% name=strcat('Trend =',32,num2str(roundn(slope,-3)),'/30 yrs;','%12s \n','p =',32, num2str(roundn(pvalue,-3)));
