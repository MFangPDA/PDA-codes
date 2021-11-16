% ****************************************************************************
% Instructions:
%     (1) Verifications on the PDA reconstruction, 20CR-V2c and MPI-ESM-P
%     The reference data include observations, CRU TS4.01 and proxy-based
%     reconstruction
%     (2) Anomaly time span: 1961-1990 CE and the millennial mean
% ****************************************************************************

function postf1_verification()

% To load subfunction in tool
addpath(genpath(pwd));

%%
% 定义全局变量
global startDA   endDA   lenDA;
global startCR20 endCR20 pathCR20;
global startCRU  endCRU  pathCRU;

% 全局变量赋值
startDA=1000;   endDA=2000;   lenDA=endDA-startDA+1;
startCR20=1851; endCR20=2014; pathCR20='./Data/20CRV2c.air.2m.mon.mean.1851.01.2014.12_2deg_NH.nc';
startCRU=1901;  endCRU=2016;  pathCRU='./Data/cru_ts4.01.1901.2016.tmp.dat_2deg_NH.nc';

%%
% 1. Load the PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;

% Anomaly with respective to the mean of 1961-1990 
Anomaly=subf_yrAnomaly(PDATas,startDA,endDA,1961,1990);
fieldPDA=Anomaly;

% NH TAS series
seriesPDA=zeros(1,lenDA);
for y=1:lenDA    % y: year
    seriesPDA(y)=mean2(fieldPDA(:,:,y));
end

% Load the standard of PDA reconstruction
stdData=load('./Result/PDAResult/stdXa.mat');
stdXa=stdData.stdXa;
% 95% significance uncertainty range
upper=seriesPDA+1.96*stdXa;
lower=seriesPDA-1.96*stdXa;

% 2. Load MPI-ESM-P simulation during 1850-2000 CE
[seriesMPI,fieldMPI]=readMPI();

% 3. Load the 20CR-V2c data (1851-2000 CE)
[field20CR,series20CR]=readCR20();

% 4.Load multi-groups NH TAS observations
% 4.1 BEST
Obs1=importdata('./Data/Observations/Berkeley_tavg_anom.txt');
mon=Obs1(:,2);
BEST=mon2anl(mon);

% 4.2 GISTEMP
Obs2=importdata('./Data/Observations/GISS_250_T2m_SST_anom.txt');
mon=Obs2(:,2);
GISTEMP=mon2anl(mon);

% 4.3 HadCRUT4
Obs3=importdata('./Data/Observations/HadCRUT4_filled-in_T2m_SST.txt');
mon=Obs3(:,2);
HadCRUT4=mon2anl(mon);

% 4.4 NCDC
Obs4=importdata('./Data/Observations/NCDC_v3_SST_T2m_anom.txt');
mon=Obs4(:,2);
NCDC=mon2anl(mon);

% Save all TAS series into a matrix (1850-2000 CE)
tasMatrix=zeros((2000-1850+1),9);
tasMatrix(:,1)=seriesPDA((1850-startDA+1):(2000-startDA+1))'; % PDA reconstruction
tasMatrix(1,2)=NaN;tasMatrix(2:end,2)=series20CR;             % 20CR-V2c
tasMatrix(:,3)=seriesMPI;                      

tasMatrix(:,4)=BEST;                          
tasMatrix(1:30,5)=NaN;tasMatrix(31:end,5)=GISTEMP;   
tasMatrix(:,6)=HadCRUT4;                             
tasMatrix(1:30,7)=NaN;tasMatrix(31:end,7)=NCDC;  
tasMatrix(:,8)=lower((1850-startDA+1):(2000-startDA+1))';     % the lower bound of PDA 
tasMatrix(:,9)=upper((1850-startDA+1):(2000-startDA+1))';     % the upper bound of PDA

%% ************************ Part one **************************************
% Verifying PDA, simulation and reanalysis series against observations
t=1850:2000;
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
fill([t,fliplr(t)],[tasMatrix(:,8)',fliplr(tasMatrix(:,9)')],[0.8 0.8 0.8],'edgealpha',0);
legend('Uncertainty','Location','SouthWest');
hold on
p=plot(t,tasMatrix(:,1),'r-',t,tasMatrix(:,2),'b-',t,tasMatrix(:,3),'c-',t,tasMatrix(:,4),'k-',...
    t,tasMatrix(:,5),'g-',t,tasMatrix(:,6),'m-',t,tasMatrix(:,7),'y-');
% Defining the legends
seriesName=cell(1,7);
seriesName(1:7)={'PDA reconstruction',...
                 '20CR-V2c',...
                 'MPI-ESM-P',...
                 'BEST',...
                 'GISTEMP',...
                 'HadCRUT4',...
                 'NCDC'};
% Plotting parameters
xlabel('Year');  ylabel('Temperature anomaly (\circC)');
legend(p(1:3),seriesName(1:3),'Location','NorthWest');
set(gca,'ylim',[-1.5 1],'fontsize',12,'fontname','cambria','fontweight','bold');
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,p(4:7),seriesName(4:7),'Location','NorthEast');
set(gca,'ylim',[-1.5 1],'fontsize',12,'fontname','cambria','fontweight','bold');

% Correlations among all series in tasMatrix (during 1880-2000 CE)
[corMatrix,pValue]=corr(tasMatrix(31:end,1:7));
save('./Result/VerResult/corPDA2tasMatrix.mat','corMatrix');
save('./Result/VerResult/pValPDA2tasMatrix.mat','pValue');

% Coefficient of efficiency values through verifying PDA, simulation and
% reanalysis against observations (during 1880-2000 CE)
% 3: PDA, simulation, reanalysis; 4: the four observational TAS series
ceMatrix=zeros(4,3);  

for i=1:3
    for j=1:4
       rec=tasMatrix(31:end,i); ref=tasMatrix(31:end,j+3); 
       ceMatrix(j,i)= CE(rec,ref);
    end
end
save('./Result/VerResult/cePDA2tasMatrix.mat','ceMatrix');

%% ************************ Part two **************************************
% Verifying PDA TAS field, 20CR-V2c TAS field and MPI-ESM-P TAS field
% against CRU TS4.01 field (only land) (1901-2000 CE)

% Load CRU TS4.01 field
fieldCRU=readCRU();

% 2.1 Computing grid-to-grid correlations
corField1=subf_corField(fieldPDA(:,:,(1901-startDA+1):(2000-startDA+1)),fieldCRU);
noNaN1=corField1(~isnan(corField1));

corField2=subf_corField(field20CR(:,:,(1901-startCR20+1):end),fieldCRU);
noNaN2=corField2(~isnan(corField2));

corField3=subf_corField(fieldMPI(:,:,(1901-1850+1):end),fieldCRU);
noNaN3=corField3(~isnan(corField3));

% Plotting the spatial distributions of grid-to-grid correlations
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
ha = tight_subplot(1,4,[.03 .02],[.02 .02],[.01 .01]); 

axes(ha(1));
corField1(isnan(corField1))=0;
plot_corField(corField1,'COR:PDA Reconstruction');

axes(ha(2));
corField2(isnan(corField2))=0;
plot_corField(corField2,'COR:20CR-V2c');

axes(ha(3));
corField3(isnan(corField3))=0;
plot_corField(corField3,'COR:MPI-ESM-P');

% Plotting the histograms of grid-to-grid correlations
figure
plotHistCOR(noNaN1,'Correlation coefficient','r'); hold on
plotHistCOR(noNaN2,'Correlation coefficient','b'); hold on
plotHistCOR(noNaN3,'Correlation coefficient','g'); hold on

% 2.2 Computing grid-to-grid coefficient of efficiency values
ceField1=subf_ceField(fieldPDA(:,:,(1901-startDA+1):(2000-startDA+1)),fieldCRU);
noNaN1=ceField1(~isnan(ceField1));

ceField2=subf_ceField(field20CR(:,:,(1901-startCR20+1):end),fieldCRU);
noNaN2=ceField2(~isnan(ceField2));

ceField3=subf_ceField(fieldMPI(:,:,(1901-1850+1):end),fieldCRU);
noNaN3=ceField3(~isnan(ceField3));

% Plotting the spatial distributions of grid-to-grid coefficient of efficiency values
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
ha = tight_subplot(1,4,[.03 .02],[.02 .02],[.01 .01]); 

axes(ha(1));
ceField1(isnan(ceField1))=-5;
plot_ceField(ceField1,'CE:PDA Reconstruction');

axes(ha(2));
ceField2(isnan(ceField2))=-5;
plot_ceField(ceField2,'CE:20CR-V2c');

axes(ha(3));
ceField3(isnan(ceField3))=-5;
plot_ceField(ceField3,'CE:MPI-ESM-P');

% Plotting the histograms of grid-to-grid coefficient of coefficiency
figure
plotHistCE(noNaN1,'Coefficient of coefficiency','r'); hold on
plotHistCE(noNaN2,'Coefficient of coefficiency','b'); hold on
plotHistCE(noNaN3,'Coefficient of coefficiency','g'); hold on

%% ************************ Part three ************************************
% Correlation coefficients between PDA reconstruction and proxy-based
% reconstrions (11 proxy-based TAS reconstructios from previous studies)

t=startDA:endDA;

% All of the time series used in this part are anomalous relative to their 
% individual millennial means
Anomaly=subf_yrAnomaly(PDATas,startDA,endDA,startDA,endDA);

% performing 11-yrs lowess smoothing on the PDA reconsturction, because the
% proxy-based reconstructions had been orighinally processed by this
% operation
Smooth=zeros(size(Anomaly));  
A=size(Anomaly);
for i=1:A(1,1)
    for j=1:A(1,2)
        Smooth(i,j,:)=smooth(Anomaly(i,j,:),11,'lowess');
    end
end
seriesPDA=zeros(1,lenDA);
for y=1:lenDA    % y: year
    seriesPDA(y)=mean2(Smooth(:,:,y));
end

upper=seriesPDA+1.96*stdXa;
lower=seriesPDA-1.96*stdXa;

% Load NH proxy-based reconstructions
ProxyRec=readProxyRec();

% Plotting all of the 1000-yrs tas series
figure
% The first, plotting the PDA uncertainty range
fill([t,fliplr(t)],[lower,fliplr(upper)],[0.8 0.8 0.8],'edgealpha',0);     
legend('Uncertainty','Location','SouthWest');   
hold on

% The second, plotting all of the tas series
p=plot(t,ProxyRec(1,:),'c-',t,ProxyRec(2,:),'c--',...
    t,ProxyRec(3,:),'b-',t,ProxyRec(4,:),'b--',...
    t,ProxyRec(5,:),'m-',t,ProxyRec(6,:),'m--',t,ProxyRec(7,:),'m-.',...
    t,ProxyRec(8,:),'g-',t,ProxyRec(9,:),'g--',t,ProxyRec(10,:),'g-.',...
    t,ProxyRec(11,:),'y-.',t,seriesPDA,'r-','linewidth',1);  hold on

Y=zeros(1,length(t));  plot(t,Y,'k--','linewidth',1);

% Defining the legends
seriesName=cell(1,12);
seriesName(1:12)={'Crowly2000','DArrigo2006',...
                 'Esper2002','Jones1998',...
                 'Mann2008(EIV,Land)','Mann2008(CPS,Land)','Mann2008(EIV,Land&Ocean)',...
                 'Shi2013(PC10+AR2)','Shi2013(CPS)','Shi2013(EIV)',...
                 'Moberg2005','PDA reconstruction'};

% Plotting parameters
xlabel('Year');  ylabel('Temperature anomaly (\circC)');
legend(p(1:6),seriesName(1:6),'Location','NorthWest');
set(gca,'ylim',[-0.8 1],'fontsize',12,'fontname','cambria','fontweight','bold');
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,p(7:12),seriesName(7:12),'Location','NorthEast');
set(gca,'ylim',[-0.8 1],'fontsize',12,'fontname','cambria','fontweight','bold');

% Computing the correlations matrix (1000-1979, because Moberg2005 is end in 1979)
dataMatrix=zeros((1979-startDA+1),12);
dataMatrix(:,1:11)=ProxyRec(:,1:(1979-startDA+1))';
dataMatrix(:,12)=seriesPDA(1:(1979-startDA+1))';
[corMatrix,pValue]=corr(dataMatrix);

% Plotting the correlations matrix as a raster picture
% figure
XVarNames=seriesName;
matrixplot(corMatrix,'XVarNames',XVarNames,'YVarNames',XVarNames,...
    'TextColor',[0,0,1],'FigStyle','Tril','ColorBar','on');
colorbar;
colormap('hot');
caxis([0.4,1]);
set(gca,'fontsize',12,'fontname','cambria','fontweight','bold');
hold on

%% ************************ Part four ************************************
% Compute the century-mean TAS value
centuryMean=zeros(12+2,lenDA);
Data=[ProxyRec;seriesPDA;lower;upper];

for i=1:9
    head=(i-1)*100+1; tail=i*100;
    for j=1:(12+2)
        centuryMean(j,head:tail)=mean(Data(j,head:tail));
    end
end

for j=1:(12+2)
    head=(10-1)*100+1;
    lastCentury=Data(j,head:end);
    centuryMean(j,head:end)=mean(lastCentury(~isnan(lastCentury)));
end

xlswrite('./centuryMean.xlsx',centuryMean);

% Plotting the centruy-mean TAS series
figure
% The first, plotting the PDA uncertainty
fill([t,fliplr(t)],[centuryMean(13,:),fliplr(centuryMean(14,:))],[0.8 0.8 0.8],'edgealpha',0);
hold on
% The second, plotting the TAS series
for i=1:12
    if (i<12)
        lineW=1;
        plot(t,centuryMean(i,:),'k-','linewidth',lineW);  hold on
    else
        lineW=2;
        plot(t,centuryMean(i,:),'r-','linewidth',lineW);  hold on
        xlabel('Year');
        ylabel('Temperature anomaly (\circC)');
        set(gca,'ylim',[-0.8 0.8],'fontsize',12,'fontname','cambria','fontweight','bold');
        Y=zeros(1,length(t)); plot(t,Y,'k--','linewidth',0.5);
    end
end
end

%% 以月平均序列求年平均序列
function aul=mon2anl(mon)
aul=zeros(length(mon)/12,1);
for i=1:length(mon)/12
    head=(i-1)*12+1; tail=i*12;
    aul(i)=mean(mon(head:tail));
end
end

%% 读取CRU4.01数据，并转换成规则格式
function cruTas=readCRU()

global startCRU endCRU pathCRU;
startSel=1901; endSel=2000;

monMean=subf_readNC(pathCRU,'tmp');
yrMean=subf_yrMean(monMean);
A=size(yrMean);
% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(yrMean(:,:,i)');
end

% 把CRU4.01中默认空值（3.0e*10）赋值为NaN
% matlab中find函数只能用于2维矩阵，超过2维就要循环了
for k=1:A(1,3)
    [rows, cols]=find(reshape(geoFormat(:,:,k),A(1,2),A(1,1))>999);
    for g=1:length(rows)
        geoFormat(rows(g),cols(g),k)=NaN;
    end
end

% 求距平值
Anomaly=subf_yrAnomaly(geoFormat,startCRU,endCRU,1961,1990);

head=startSel-startCRU+1; tail=endSel-startCRU+1;
cruTas=Anomaly(:,:,head:tail);
end

%% 读取20CR-V2c的数据(AD1850-2000)
function [CRV2c,CRV2cSeries]=readCR20()

global startCR20 endCR20 pathCR20;

% 读取20CR-V2c数据，把其作为参照
CR20V2c_monMean=subf_readNC(pathCR20,'air');
CR20V2c_yrMean=subf_yrMean(CR20V2c_monMean);
% 求距平值
CR20V2c_Anomaly=subf_yrAnomaly(CR20V2c_yrMean,startCR20,endCR20,1961,1990);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
B=size(CR20V2c_Anomaly);
geoFormat=zeros(B(1,2),B(1,1),B(1,3));
for i=1:B(1,3)
    geoFormat(:,:,i)=flipud(CR20V2c_Anomaly(:,:,i)');
end
CR20Tas=geoFormat;

head=1851-startCR20+1; tail=2000-startCR20+1;
CRV2c=CR20Tas(:,:,head:tail);

CRV2cSeries=zeros(1,(2000-1851+1));
for i=1:(2000-1851+1)
    CRV2cSeries(i)=mean2(CRV2c(:,:,i));
end
end

%% 读取MPI-ESM-P的数据
function [series,field]=readMPI()
% 读MPI-ESM-P在185001-200512的数据
% 提取出1850-2000年的北半球气温序列与气温场
% 读取MPI数据，并进行空间规则化转化
tas_Amon=subf_readNC('./Data/tas_Amon_2deg_MPI_185001-200512.nc','tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
tas_yrAnomaly=subf_yrAnomaly(geoFormat,1850,2005,1961,1990);
field=tas_yrAnomaly(1:45,:,1:(2000-1850+1));
for i=1:(2000-1850+1)
    series(i)=mean2(field(:,:,i));
end
end

%% 读取代用资料重建的结果(AD1000-2000年之间的结果)
function ProxyRec=readProxyRec()
smoothSpan=11;

% 注意：所有数据都求出相对于（1000-终止年份）的距平值
% 有些序列的长度达不到2000年的话，就用NaN补充上
ProxyRec=zeros(11,(2000-1000+1));
%
%
% (1) Crowly(2000)Northern Hemisphere(AD1000-1993), Anomaly w.r.t.don't know
% 注意：Crowly(2000)数据从上往下是年份升序，正常方式，不用做上下反转
Crowly=importdata('./Data/ProxyRec/Crowly_2000_Northern Hemisphere_AD1000-1993.txt');
% 求距平值
anomalSeries=Crowly((1000-1000+1):end,2)-mean(Crowly((1000-1000+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(1,1:(1993-1000+1))=smoothSeries;
ProxyRec(1,(1994-1000+1):end)=NaN;
%
%
% (2) D'Arrigo(2006)Northern Hemisphere(AD713-1995), Anomaly w.r.t.1961-1990
% 注意：D'Arrigo(2006)数据从上往下是年份降序，不正常方式
DArrigo=importdata('./Data/ProxyRec/DArrigo_2006_Northern Hemisphere_AD713-1995.txt');
DArrigo=flipud(DArrigo);
% 求距平值
anomalSeries=DArrigo((1000-713+1):end,2)-mean(DArrigo((1000-713+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(2,1:(1995-1000+1))=smoothSeries;
ProxyRec(2,(1996-1000+1):end)=NaN;
%
%
% (3) Esper(2002)Northern Hemisphere(AD831-1992), Unitless Index Value
% 注意：Esper(2002)数据从上往下是年份降序，不正常方式
Esper=importdata('./Data/ProxyRec/Esper_2002_Northern Hemisphere_AD831-1992.txt');
Esper=flipud(Esper);
% 求距平值
anomalSeries=Esper((1000-831+1):end,2)-mean(Esper((1000-831+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(3,1:(1992-1000+1))=smoothSeries;
ProxyRec(3,(1993-1000+1):end)=NaN;
%
%
% (4) Jones(1998)Northern Hemisphere(AD1000-1991), Anomaly w.r.t.1961-1990
% 注意：Jones(1998)数据从上往下是年份降序，不正常方式
Jones=importdata('./Data/ProxyRec/Jones_1998_Northern Hemisphere_AD1000-1991.txt');
Jones=flipud(Jones);
% 求距平值
anomalSeries=Jones((1000-1000+1):end,2)-mean(Jones((1000-1000+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(4,1:(1991-1000+1))=smoothSeries;
ProxyRec(4,(1992-1000+1):end)=NaN;
%
%
% (5) Mann(2008)Northern Hemishpere(AD300-2006)Land Only_EIV, Anomaly w.r.t1961-1990
% 注意：Mann(2008)数据从上往下是年份降序，不正常方式
Mann=importdata('./Data/ProxyRec/Mann_2008_Northern Hemishpere_AD300-2006_Land Only_EIV.txt');
Mann=flipud(Mann);
% 求距平值
anomalSeries=Mann((1000-300+1):(2000-300+1),2)-mean(Mann((1000-300+1):(2000-300+1),2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(5,1:(2000-1000+1))=smoothSeries;
%
%
% (6) Mann(2008)Northern Hemisphere(AD200-1995)Land Only_CPS, Anomaly w.r.t1961-1990
% 注意：Mann(2008)数据从上往下是年份降序，不正常方式
Mann=importdata('./Data/ProxyRec/Mann_2008_Northern Hemisphere_AD200-1995_Land Only_CPS.txt');
Mann=flipud(Mann);
% 求距平值
anomalSeries=Mann((1000-200+1):end,2)-mean(Mann((1000-200+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(6,1:(1995-1000+1))=smoothSeries;
ProxyRec(6,(1996-1000+1):end)=NaN;
%
%
% (7) Mann(2008)Northern Hemisphere(AD300-2006)Land and Ocean_EIV, Anomaly w.r.t1961-1990
% 注意：Mann(2008)数据从上往下是年份降序，不正常方式
Mann=importdata('./Data/ProxyRec/Mann_2008_Northern Hemisphere_AD300-2006_Land and Ocean_EIV.txt');
Mann=flipud(Mann);
% 求距平值
anomalSeries=Mann((1000-300+1):(2000-300+1),2)-mean(Mann((1000-300+1):(2000-300+1),2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(7,1:(2000-1000+1))=smoothSeries;
%
%
% (8) Shi(2013)Northern Hemisphere(AD1000-1998), Anomaly w.r.t.1961-1990
% 注意：Shi(2013)数据从上往下是年份升序，正常方式
Shi=importdata('./Data/ProxyRec/Shi_2013_Northern Hemisphere_AD1000-1998.txt');
% 求距平值
anomalSeries=Shi((1000-1000+1):end,2)-mean(Shi((1000-1000+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(8,1:(1998-1000+1))=smoothSeries;
ProxyRec(8,(1999-1000+1):end)=NaN;
%
% 求距平值
anomalSeries=Shi((1000-1000+1):end,5)-mean(Shi((1000-1000+1):end,5));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(9,1:(1998-1000+1))=smoothSeries;
ProxyRec(9,(1999-1000+1):end)=NaN;
%
% 求距平值
anomalSeries=Shi((1000-1000+1):end,6)-mean(Shi((1000-1000+1):end,6));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(10,1:(1998-1000+1))=smoothSeries;
ProxyRec(10,(1999-1000+1):end)=NaN;
%
%
% (9) Moberg(2005)Northern Hemisphere(AD1-1979), Anomaly w.r.t1961-1990
% 注意：Mann(2008)数据从上往下是年份降序，不正常方式
Moberg=importdata('./Data/ProxyRec/Moberg_2005_Northern Hemisphere_AD1-1979.txt');
Moberg=flipud(Moberg);
% 求距平值
anomalSeries=Moberg((1000-1+1):end,2)-mean(Moberg((1000-1+1):end,2));
% 求平滑值
smoothSeries=smooth(anomalSeries,smoothSpan,'lowess');
% 序列赋值
ProxyRec(11,1:(1979-1000+1))=smoothSeries;
ProxyRec(11,(1980-1000+1):end)=NaN;
end

%% 求两个序列（rec,ref）之间的CE
function ceValue=CE(rec,ref)
    ceValue=1-(sum((ref-rec).*(ref-rec)))/(sum((ref-mean(ref)).*(ref-mean(ref))));
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
set(gca,'YLim',[0 500]); set(gca,'YTick',0:100:500);
xlabel(paraName); ylabel('Frequency');
text(-0.9,400,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',12,'fontweight','bold'); 
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
set(gca,'YLim',[0 1000]); set(gca,'YTick',0:200:1000);
xlabel(paraName); ylabel('Frequency');
text(-4.8,800,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',12,'fontweight','bold'); 
% set(gcf,'color','none');   % 图形背景设为无色
set(gca,'color','none');     % 坐标轴背景为无色，这条很重要，通常图形背景的白色（实际为坐标轴默认背景色），设置为透明，就可以进行多幅图的叠加
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end