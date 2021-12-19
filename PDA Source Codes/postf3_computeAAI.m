% ****************************************************************************
% Description:
%     (1) Computing the AA index values based on different datasets,
%     including PDA, multi-models millennium simulations in PMIP3
%     (2) All of the time series are anomalous relative to their individual
%     means of 1961-1990
%     Note that: the choice of the reference time span do not affect the AAI computation
% ****************************************************************************
function postf4_computeAAI()

% To load subfunction in tool
addpath(genpath(pwd));
%%
%
% Defining global variables
global startDA    endDA     lenDA    startAA    endAA   intervalAA;

% Assigning values to global variables
startDA=1000;     endDA=2000;    lenDA=endDA-startDA+1;
t=startDA:endDA;
startAA=1015; endAA=1985;  intervalAA=30;

% 1. Loading PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;
anomalyPDA=subf_yrAnomaly(PDATas,startDA,endDA,1961,1990);
PDA=anomalyPDA;
%%
%
% 2.To compute AA index based on PDA reconstruction and PMIP3 simulations
% Defining a matrix to store all AA index series
AAindices=zeros(endAA-startAA+1,13);
AAtrends=zeros(endAA-startAA+1,13);

% 2.1 AA index derived from PDA reconstruction
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(PDA);  PDA_NHseries=NHseries;
AAindices(:,1)=yearAAI;    AAtrends(:,1)=trendValue;

% 2.2 AA index derived from Goosse et al.,2012 based on data assimilation
tasField=readGoosse2012();
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,2)=yearAAI;    AAtrends(:,2)=trendValue;

% 2.3 AA index derived from Steiger et al., 2018 based on data assimilation
[LMRNHseries,LMRArcNHseries]=readLMRv2();
[yearAAI,trendValue]=computeAAILMRv2(LMRNHseries,LMRArcNHseries);
AAindices(:,3)=yearAAI;    AAtrends(:,3)=trendValue;

% 2.4 AA index derived from PMIP3 simulations (1000-2000 AD)
CMIP5_NHseries=zeros(1,1001); CMIP5_ARCseries=zeros(1,1001);

tasField=readPMIP3('BCC');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,4)=yearAAI;    AAtrends(:,4)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readPMIP3('CCSM4');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,5)=yearAAI;    AAtrends(:,5)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readPMIP3('FGOALSg2');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,6)=yearAAI;    AAtrends(:,6)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readPMIP3('GISS');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,7)=yearAAI;    AAtrends(:,7)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readPMIP3('IPSL');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,8)=yearAAI;    AAtrends(:,8)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readPMIP3('MPI');
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,9)=yearAAI;    AAtrends(:,9)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

tasField=readCSIRO();
[NHseries,ARCseries,yearAAI,trendValue]=computeAAI(tasField);
AAindices(:,10)=yearAAI;    AAtrends(:,10)=trendValue;
CMIP5_NHseries=CMIP5_NHseries+NHseries; CMIP5_ARCseries=CMIP5_ARCseries+ARCseries;

% ����7��ģʽ����������ƽ��ֵ
CMIP5_NHseries=CMIP5_NHseries/7; CMIP5_ARCseries=CMIP5_ARCseries/7;

tasField=readFGOALSg1();  % ע�⣺FGOALSg1û��2000�������
[NHseries,ARCseries,yearAAI,trendValue]=computeAAIFGOALSg1(tasField);
AAindices(:,11)=yearAAI;    AAtrends(:,11)=trendValue;
% 1000-1999 AD
CMIP5_NHseries(1:1000)=(CMIP5_NHseries(1:1000)+NHseries)/2;
CMIP5_ARCseries(1:1000)=(CMIP5_ARCseries(1:1000)+ARCseries)/2;

tasField=readHadCM3();   % ע�⣺HadCM3û��1851-1859������ݣ����ʱ���������0������
[NHseries,ARCseries,yearAAI,trendValue]=computeAAIHadCM3(tasField);
AAindices(:,12)=yearAAI;    AAtrends(:,12)=trendValue;
% 1000-1850 AD
CMIP5_NHseries(1:851)=(CMIP5_NHseries(1:851)+NHseries(1:851))/2;
CMIP5_ARCseries(1:851)=(CMIP5_ARCseries(1:851)+ARCseries(1:851))/2;
% 1860-2000 AD
CMIP5_NHseries(861:end)=(CMIP5_NHseries(861:end)+NHseries(861:end))/2;
CMIP5_ARCseries(861:end)=(CMIP5_ARCseries(861:end)+ARCseries(861:end))/2;

% �����ģʽ����ƽ����AAI
[AAI_CMIP5,trendValue_CMIP5]=computeAAICMIP5(CMIP5_NHseries,CMIP5_ARCseries);
AAindices(:,13)=AAI_CMIP5;    AAtrends(:,13)=trendValue_CMIP5;

% 2.5 To call the Mann-Kendall trend test function to check the trend of AA index
% during the past millennium
dataName=cell(1,13);
dataName(1:13)={'PDA-based reconstruction','Goosse2012','LMRv2','BCC-CSM1','CCSM4',...
    'FGOALS-s2','GISS-E2-R','IPSL-CM5A-LR','MPI-ESM-P','CSIRO Mk3L-1-2','FGOALS-g1','HadCM3',...
    'Multi-model mean'};
Slopes=zeros(1,13);
pValues=zeros(1,13);
for i=1:13
    if i==11
        series=AAindices(1:(endAA-startAA),i);
    else
        series=AAindices(:,i);
    end
    [slope,pValue]=MKTrendDetection(series,dataName{i});
    Slopes(i)=slope;
    pValues(i)=pValue;
end

% % 2.6 Defining a cell matrix to save all AAI data except for the AA index of 20CR-V2c
AAMatrix=cell((endAA-startAA+1)+1,14);
AAMatrix(1,1)={'Year'}; AAMatrix(1,2:14)=dataName;

for i=1:(endAA-startAA+1)
    period={num2str(i)};
    AAMatrix(i+1,1)=period;
    for j=1:13
        AAMatrix(i+1,j+1)={num2str(roundn(AAindices(i,j),-3))};
    end
end
xlswrite('allAAIMatrix.xlsx', AAMatrix)

% 2.7 Poltting the millennial AA index series and its trend
figure
t=startAA:endAA;
set(gcf,'unit','inches','position',[0.2,0.2,7.087,5]);
ha = tight_subplot(2,2,[.1 .05],[.06 .04],[.08 .04]);

axes(ha(1));
plotAAI(t,AAindices(:,1),AAtrends(:,1),Slopes(1),pValues(1),'a',dataName{1}); hold on
axes(ha(2));
plotAAI(t,AAindices(:,2),AAtrends(:,2),Slopes(2),pValues(2),'b',dataName{2}); hold on
axes(ha(3));
plotAAI(t,AAindices(:,3),AAtrends(:,3),Slopes(3),pValues(3),'c',dataName{3}); hold on
xlabel('Year (AD)'); ylabel({'AA index'});
axes(ha(4));
plotAAI(t,AAindices(:,13),AAtrends(:,13),Slopes(13),pValues(13),'d',dataName{13}); hold on
xlabel('Year (AD)'); ylabel({'AA index'});
set(gca,'ylim',[0 3.5],'fontsize',7,'fontname','Arial','linewidth',1,'Tickdir','out');

figure
set(gcf,'unit','inches','position',[0.2,0.2,7.087,5]);
ha = tight_subplot(1,1,[.1 .05],[.06 .04],[.08 .04]);
axes(ha(1));
legendCell=cell(1,10);
for i=1:10
    if (i==10)  % multi-models mean
        plot(t,AAindices(:,i+3),'b-','linewidth',2); hold on
        plot(t,AAtrends(:,i+3),'b--','linewidth',2); hold on
    else
        plot(t,AAindices(:,i+3),'linewidth',0.5); hold on
    end
    xlabel('Year (AD)');
    ylabel({'AA index'});
    set(gca,'ylim',[0 3.5],'fontsize',9,'fontname','Arial','linewidth',1,'Tickdir','out'); 
    
    if pValues(i+2)<0.001  
        name=strcat('Trend =',32,num2str(roundn(Slopes(i+3),-2)),'/100 yrs,','p < 0.001');
    else
        name=strcat('Trend =',32,num2str(roundn(Slopes(i+3),-2)),'/100 yrs,','p =',32, num2str(roundn(pValues(i+2),-3)));
    end
    legendCell(1,i)=strcat(dataName(i+3),'(',name,')');
    hold on
end
legend(legendCell,'Box','off');
text(startAA,3.6,'A','fontsize',7,'fontweight','bold');
text(startAA+20,3.6,'CMIP5 multi-model','fontsize',7);

% 2.8 �����º�AAI����һ��ͼ���� (����)
figure
AAI_PDA=zeros(lenDA,1);
AAI_PDA(1:15)=NaN; AAI_PDA(16:986)=AAindices(:,1); AAI_PDA(987:1001)=NaN;

year=1000:2000;
[AX,H1,H2]=plotyy(year,PDA_NHseries,year, AAI_PDA,'plot');

% ��ȡ�����ᡢͼ����
set(AX,'FontSize',9,'FontName','Arial'); % ����x�ᡢ��y�ᡢ��y��̶��ֺź�����

% ����
set(AX(1),'XColor','k','YColor',[0.5451 0 0],'fontsize',9,'fontname','Arial','Linewidth',1,'Tickdir','out');
% ����
set(AX(2),'XColor','k','YColor',[0 0 0.5451],'fontsize',9,'fontname','Arial','Linewidth',1,'Tickdir','out');

HH1=get(AX(1),'Ylabel');
set(HH1,'color',[0.5451 0 0]);     set(HH1,'String','NH Temperature anomaly (\circC)');
set(H1,'LineStyle','-');  set(H1,'Linewidth',1.5); set(H1,'color',[0.5451 0 0]);

HH2=get(AX(2),'Ylabel');
set(HH2,'color',[0 0 0.5451]);     set(HH2,'String','AA index');
set(H2,'LineStyle','-');  set(H2,'Linewidth',1.5); set(H2,'color',[0 0 0.5451]);

set(AX(1),'ylim',[-1.5 0.5],'yTick',[-1.5:0.5:0.5]);
set(AX(2),'ylim',[1 2.5],'yTick',[1:0.5:2.5]);

xlabel('Year (AD)'); 
legend('NH temperautre anomaly','AA index');

% figure
% fill([t,fliplr(t)],[modelMean(:,2)',fliplr(modelMean(:,3)')],[0.52941 0.80784 0.98039],'edgealpha',0); hold on
% plot(t,modelMean(:,1),'color',[0 0 0.5451]);
% % fill�﷨��
% % fill([t,fliplr(t),[upperBoundary,fliplr(lowerBoundary)],[0.5451 0 0],'edgealpha',0)
% % t: ������
% % upperBoundary:������
% % lowerBoundary:������
end

%% *******************��ȡ��������***********************
% ��ȡBCC,CCSM4,FGOALg2,GISS,IPSL,MPI��ͨ�ó���
function tasField=readPMIP3(modelname)
% ��������ʱ�Σ�085001-200512
% ��ȡԭʼ���ݣ������пռ����ת��
tas_Amon=subf_readNC(strcat('./PMIP3_en/tas_Amon_2deg_',modelname,'_085001_200512_NH.nc'),'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% ֻ��ȡ1000-2000������
field=geoFormat(:,:,(1000-850+1):(2000-850+1));
tasField=subf_yrAnomaly(field,1000,2000,1961,1990);
end

%% ��ȡFGOALSg1�����ݣ�ע��2000�������ȱʧ
function tasField=readFGOALSg1()
% FGOALS-g1������ʱ�Σ�100001-199912
% ��ȡԭʼ���ݣ������пռ����ת��
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_FGOALSg1_100001_199912_NH.nc','tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% ֻ��1000-1999������
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,1999,1961,1990);
end

%% ��ȡHadCM3�����ݣ�ע��1851-1859�������ȱʧ
function tasField=readHadCM3()
% HadCM3���ݷ�Ϊ���Σ�185912-200512��085001-185012��ȱ��1851-1859������
% ��ȡԭʼ���ݣ������пռ����ת��
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_HadCM3_185912_200512_NH.nc','tas');
% ȥ��1859-12����µ����ݣ�����������ƽ��
tas_yrMean=subf_yrMean(tas_Amon(:,:,2:end))-273.5;
A=size(tas_yrMean);

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
field1=geoFormat(:,:,1:(2000-1860+1)); % �ǵø�Ϊ1860
%
% ������ȱ�ٵ�1851-1859����9���������ʱ���ܣ���0����
D=size(field1);
field1_add=zeros(D(1,1),D(1,2),D(1,3)+9);
field1_add(:,:,1:9)=0;
field1_add(:,:,10:end)=field1;
%
% ��������ȡǰ���
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_HadCM3_085001_185012_NH.nc','tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
B=size(tas_yrMean);
%
% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(B(1,2),B(1,1),B(1,3));
for i=1:B(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
field2=geoFormat(:,:,(1000-850+1):end);
%
% ǰ���������ݺϲ�
fieldHad=zeros(B(1,2),B(1,1),(2000-1000+1));
fieldHad(:,:,1:(1850-1000+1))=field2;
fieldHad(:,:,(1851-1000+1):(2000-1000+1))=field1_add;
tasField=subf_yrAnomaly(fieldHad,1000,2000,1961,1990);
end

%% ��ȡCSIROMk3L-1-2������
function tasField=readCSIRO()
% CSIROMk3L-1-2������ʱ�Σ�100001-200012
% ��ȡԭʼ���ݣ������пռ����ת��
tas_Amon=subf_readNC('./PMIP3_en/tsc_Amon_2deg_CSIRO_100001-200012_NH.nc','tsc');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));
% ����⣬���ݽṹ��PMIP3������һ��

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,2000,1961,1990);
end

%% ��ȡGoosse2012ͬ�������ݣ��Ѿ�����ֱ���
function tasField=readGoosse2012()
% Goosse2012������ʱ�Σ�1000-2000����ֱ��ʣ��Ѿ����˾�ƽ
% ��ȡԭʼ���ݣ������пռ����ת��
tas_yrMean=subf_readNC('./PMIP3_en/ts_Ayr_2deg_assim_loveclim_mann2009_1000_2000_NH.nc','ts');
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));
% ����⣬���ݽṹ��PMIP3������һ��

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,2000,1961,1990);
end

%%
% ******************������ֱ�ļ��طŴ�ָ��**************************
% ����ģʽ��AAI��������PDA,BCC,CCSM4,FGOALg2,GISS,IPSL,MPI,CSIRO��ģʽ
function [NHseries,ARCseries,yrAAI,trendValue]=computeAAI(tasField)

global endAA startAA intervalAA;

% ����������ֹʱ�䣺1000-2000��
A=size(tasField);
NHseries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    NHseries(y)=mean2(tasField(:,:,y));
end
% 60-90��N��30-60��N��0-30��N��0-90��N�ĸ�γ�ȴ�
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % ÿ30��һ��γ�ȴ����ܼ�3��γ�ȴ�
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
ARCseries=zonalSeries(1,:);

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-2000,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1985��
lenAAI=endAA-startAA+1;
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;       % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
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

%% ����ģʽ��AAI��������FGOALS-g1ģʽ
% FGOALS-g1ģʽģ����ʱ�䳤�ȣ�1000-1999��
function [NHseries,ARCseries,yrAAI,trendValue]=computeAAIFGOALSg1(tasField)

global endAA startAA intervalAA;

A=size(tasField);
NHseries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    NHseries(y)=mean2(tasField(:,:,y));
end
% 60-90��N��30-60��N��0-30��N��0-90��N�ĸ�γ�ȴ�
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % ÿ30��һ��γ�ȴ����ܼ�3��γ�ȴ�
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
ARCseries=zonalSeries(1,:);

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-1999,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1984��
lenAAI=(endAA-1-startAA+1);
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;       % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    yrAAI(i)=b(2);
end
%
% ��AAI���߽����������
tt=1:length(yrAAI);
[p,S]=polyfit(tt',yrAAI,1);  % ������������
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),length(yrAAI));

% !!! Ϊ�˱�֤���ݵĳ���һ�£�1985���ֵ��NaN����
yrAAI=[yrAAI;NaN];
trendValue=[trendValue';NaN];
end

%% ����ģʽ��AAI��������HadCM3ģʽ
function [NHseries,ARCseries,yrAAI,trendValue]=computeAAIHadCM3(tasField)

global endAA startAA intervalAA;

% HadCM3���ݷ�Ϊ���Σ�185912-200512��085001-185012��ȱ��1851-1859������
% ����������AAI��ʱ��ȱ�ٵ�1851-1859����9���������NaN����
A=size(tasField);
NHseries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    NHseries(y)=mean2(tasField(:,:,y));
end
NHseries((1851-1000+1):(1859-1000+1))=NaN;

% 60-90��N��30-60��N��0-30��N��0-90��N�ĸ�γ�ȴ�
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % ÿ30��һ��γ�ȴ����ܼ�3��γ�ȴ�
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
ARCseries=zonalSeries(1,:);
ARCseries((1851-1000+1):(1859-1000+1))=NaN;

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-2000,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1985��
lenAAI=(endAA-startAA+1);
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;       % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
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

% �����ģʽ����ƽ����AAI
function [yrAAI,trendValue]=computeAAICMIP5(NHseries,ARCseries)

global endAA startAA intervalAA;

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-2000,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1985��
lenAAI=(endAA-startAA+1);
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;  % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
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

% ����Steiger 2018ͬ�������AAI
function [yrAAI,trendValue]=computeAAILMRv2(NHseries,ARCseries)

global endAA startAA intervalAA;

% �û����ع�ķ�������AAI������������Ϊ31��
% 1000-2000,31��Ĵ��ڣ���һ��AAIֵΪ1015�꣬���һ��AAIֵΪ1985��
lenAAI=(endAA-startAA+1);
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;  % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
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


function plotAAI(t,AAI,trendValue,slope,pvalue,index,model)
plot(t,AAI,'color',[0.5451 0 0],'linewidth',1.5); hold on
plot(t,trendValue,'color',[0 0 0.5451],'linewidth',1.5);
% name=strcat('Trend =',32,num2str(roundn(slope,-3)),'/30 yrs;','%12s \n','p =',32, num2str(roundn(pvalue,-3)));

if pvalue<0.001  
    name=strcat('Trend =',32,num2str(roundn(slope,-2)),'/100 yrs,','p < 0.001');
else
    name=strcat('Trend =',32,num2str(roundn(slope,-2)),'/100 yrs,','p =',32, num2str(roundn(pvalue,-3)));
end
xlabel('Year (AD)');
ylabel({'AA index'});
set(gca,'ylim',[1 3],'fontsize',7,'fontname','Arial','linewidth',1,'Tickdir','out');
text(1015,3.1,index,'fontsize',7,'fontweight','bold');
text(1135,3.1,model,'fontsize',7);
text(1020,0.7,name,'fontsize',7,'fontname','Arial','color',[0 0 0.5451]);
hold on
end
