% ****************************************************************************
% Description:
%     (1) Computing the AA index values based on different datasets,
%     including PDA, multi-models millennium simulations in PMIP3 
%     (2) All of the time series are anomalous relative to their individual millennial means
%     Note that: the choice of the reference time span do not affect the AAI computation
%   
%     注意：只有用于作图展示的数据，才可以进行平滑；而用于计算其他的次生变量的数据，在计算次生变量
%     之前，是不能进行平滑的，必须保持原始的数据
% ****************************************************************************
function postf4_othercomputeAAI()

% 调用子文件夹Tool中的m文件
addpath(genpath(pwd));
%%
%
% Defining global variables
global startDA    endDA     lenDA;
global startCR20  endCR20   pathCR20  lenCR20;

% Assigning values to global variables
startDA=1000;     endDA=2000;    lenDA=endDA-startDA+1;
startCR20=1851;   endCR20=2014;  lenCR20=endCR20-startCR20+1;
pathCR20='./Data/20CRV2c.air.2m.mon.mean.1851.01.2014.12_2deg_NH.nc';
t=startDA:endDA;

% 1. Loading PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;

anomalyPDA=subf_yrAnomaly(PDATas,startDA,endDA,1961,1990);

PDA=anomalyPDA;
%%
%
% 2.To compute AA index based on PDA reconstruction and PMIP3 simulations

% Defining a matrix to store all AA index series
allAAI=zeros(33,11);

% 2.1 AA index derived from PDA reconstruction
[AAI_PDA,trendValue_PDA,AAI30yr]=computeAAI(PDA);                allAAI(:,1)=AAI30yr;

AAIb1850=mean(AAI30yr(1:28));
AAIa1850=mean(AAI30yr(29:33));

% 2.2 AA index derived from Goosse et al., 2009 based on data assimilation
tasField=readGoosse2012();      
[AAI_Goosse,trendValue_Goosse,AAI30yr]=computeAAI(tasField);     allAAI(:,2)=AAI30yr;

% 2.3 AA index derived from PMIP3 simulations
tasField=readPMIP3('BCC');      
[AAI_BCC,trendValue_BCC,AAI30yr]=computeAAI(tasField);           allAAI(:,3)=AAI30yr;         
tasField=readPMIP3('CCSM4');    
[AAI_CCSM4,trendValue_CCSM4,AAI30yr]=computeAAI(tasField);       allAAI(:,4)=AAI30yr;  
tasField=readPMIP3('FGOALSg2'); 
[AAI_FGOALSg2,trendValue_FGOALSg2,AAI30yr]=computeAAI(tasField); allAAI(:,5)=AAI30yr;
tasField=readPMIP3('GISS');     
[AAI_GISS,trendValue_GISS,AAI30yr]=computeAAI(tasField);         allAAI(:,6)=AAI30yr;      
tasField=readPMIP3('IPSL');     
[AAI_IPSL,trendValue_IPSL,AAI30yr]=computeAAI(tasField);         allAAI(:,7)=AAI30yr;       
tasField=readPMIP3('MPI');      
[AAI_MPI,trendValue_MPI,AAI30yr]=computeAAI(tasField);           allAAI(:,8)=AAI30yr;     

tasField=readCSIRO();           
[AAI_CSIRO,trendValue_CSIRO,AAI30yr]=computeAAI(tasField);       allAAI(:,9)=AAI30yr;   

tasField=readFGOALSg1();
[AAI_FGOALSg1,trendValue_FGOALSg1,AAI30yr]=computeAAIFGOALSg1(tasField); allAAI(:,10)=AAI30yr;

tasField=readHadCM3();  
[AAI_HadCM3,trendValue_HadCM3,AAI30yr]=computeAAIHadCM3(tasField);       allAAI(:,11)=AAI30yr;   

% 2.4 AA index derived from 20CR-V2c reanalysis during 1850-2000 CE 
tasField=readCR20();            AAI_CR20=computeAAI20CR(tasField);

% 2.5 To call the Mann-Kendall trend test function to check the trend of AA index
% during the past millennium
dataName=cell(1,12);
dataName(1:12)={'PDA','Goosse2012','BCC-CSM1','CCSM4',...
               'FGOALS-s2','FGOALS-g1','GISS-E2-R','IPSL-CM5A-LR',...
               'MPI-ESM-P','HadCM3','CSIRO Mk3L-1-2','Multimodel mean'};           
Slopes=zeros(1,12);
pValues=zeros(1,12);
for i=1:11
    series=allAAI(:,i);
    [slope,pValue]=MKTrendDetection(series,dataName{i});
    Slopes(i)=slope;
    pValues(i)=pValue;
end

% 多模式平均
series=mean(allAAI,2);
[slope,pValue]=MKTrendDetection(series,dataName{12});
Slopes(12)=slope;
pValues(12)=pValue;

% % 2.6 Defining a cell matrix to save all AAI data except for the AA index of 20CR-V2c
AAMatrix=cell(34,12);
AAMatrix(1,1)={'Period'}; AAMatrix(1,2:12)=dataName;

for i=1:33
    if i==1
        period={'1000-1040'};
    else
        period={strcat(num2str((i-2)*30+1+1040),'-',num2str(1040+(i-1)*30))};
    end
    AAMatrix(i+1,1)=period;
    
    for j=1:11
       AAMatrix(i+1,j+1)={num2str(allAAI(i,j))};
    end
end
% save('./Result/AAIanalysis/AAMatrix.mat','AAMatrix');
xlswrite('./Result/AAIanalysis/AAMatrix_all.xlsx', AAMatrix)

% 2.7 Poltting the millennial AA index series and its trend
figure
subplot(3,4,1)
plotAAI(t,AAI_PDA,trendValue_PDA,'k-',Slopes(1),pValues(1)); 
subplot(3,4,2)
plotAAI(t,AAI_Goosse,trendValue_Goosse,'k-',Slopes(2),pValues(2)); 
subplot(3,4,3)
plotAAI(t,AAI_BCC,trendValue_BCC,'k-',Slopes(3),pValues(3));
subplot(3,4,4)
plotAAI(t,AAI_CCSM4,trendValue_CCSM4,'k-',Slopes(4),pValues(4));
subplot(3,4,5)
plotAAI(t,AAI_FGOALSg2,trendValue_FGOALSg2,'k-',Slopes(5),pValues(5));
subplot(3,4,6)
plotAAI(t,AAI_FGOALSg1,trendValue_FGOALSg1,'k-',Slopes(6),pValues(6));
subplot(3,4,7)
plotAAI(t,AAI_GISS,trendValue_GISS,'k-',Slopes(7),pValues(7));
subplot(3,4,8)
plotAAI(t,AAI_IPSL,trendValue_IPSL,'k-',Slopes(8),pValues(8));
subplot(3,4,9)
plotAAI(t,AAI_MPI,trendValue_MPI,'k-',Slopes(9),pValues(9));
subplot(3,4,10)
plotAAI(t,AAI_HadCM3,trendValue_HadCM3,'k-',Slopes(10),pValues(10));
subplot(3,4,11)
plotAAI(t,AAI_CSIRO,trendValue_CSIRO,'k-',Slopes(11),pValues(11));
subplot(3,4,12)
plotAAI(t,series,trendValue_CSIRO,'k-',Slopes(12),pValues(12));
end

%% *******************读取各类数据***********************
% 包括：再分析资料、PMIP3集合模拟
% 读取20CR-V2c的数据(AD1850-2000)
function fieldCRV2c=readCR20()

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
fieldCRV2c=CR20Tas(:,:,head:tail);
end

%% 读取BCC,CCSM4,FGOALg2,GISS,IPSL,MPI的通用程序
function tasField=readPMIP3(modelname)
% 以上数据时段：085001-200512
% 读取原始数据，并进行空间规则化转化
tas_Amon=subf_readNC(strcat('./PMIP3_en/tas_Amon_2deg_',modelname,'_085001_200512_NH.nc'),'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% 只提取1000-2000的数据
field=geoFormat(:,:,(1000-850+1):(2000-850+1));
tasField=subf_yrAnomaly(field,1000,2000,1961,1990);
end

%% 读取FGOALSg1的数据，注意2000年的数据缺失
function tasField=readFGOALSg1()
% FGOALS-g1的数据时段：100001-199912
% 读取原始数据，并进行空间规则化转化
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_FGOALSg1_100001_199912_NH.nc','tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% 只有1000-1999的数据
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,1999,1961,1990);
end

%% 读取HadCM3的数据，注意1851-1859年的数据缺失
function tasField=readHadCM3()
% HadCM3数据分为两段：185912-200512；085001-185012，缺少1851-1859的数据
% 读取原始数据，并进行空间规则化转化
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_HadCM3_185912_200512_NH.nc','tas');
% 去掉1859-12这个月的数据，才能再求年平均
tas_yrMean=subf_yrMean(tas_Amon(:,:,2:end))-273.5;
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
field1=geoFormat(:,:,1:(2000-1860+1)); % 记得改为1860
%
% ！！！缺少的1851-1859年这9年的数据暂时不管，用0代替
% ！！！计算AAI的时候，1851-1880年的AAI用1860-1880年的数据代替
D=size(field1);
field1_add=zeros(D(1,1),D(1,2),D(1,3)+9);
field1_add(:,:,1:9)=0;
field1_add(:,:,10:end)=field1;
%
% 接下来读取前半段
tas_Amon=subf_readNC('./PMIP3_en/tas_Amon_2deg_HadCM3_085001_185012_NH.nc','tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
B=size(tas_yrMean);
%
% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(B(1,2),B(1,1),B(1,3));
for i=1:B(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
field2=geoFormat(:,:,(1000-850+1):end);
%
% 前后两段数据合并
fieldHad=zeros(B(1,2),B(1,1),(2000-1000+1));
fieldHad(:,:,1:(1850-1000+1))=field2;
fieldHad(:,:,(1851-1000+1):(2000-1000+1))=field1_add;
tasField=subf_yrAnomaly(fieldHad,1000,2000,1961,1990);
end

%% 读取CSIROMk3L-1-2的数据
function tasField=readCSIRO()
% CSIROMk3L-1-2的数据时段：100001-200012
% 读取原始数据，并进行空间规则化转化
tas_Amon=subf_readNC('./PMIP3_en/tsc_Amon_2deg_CSIRO_100001-200012_NH.nc','tsc');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));

% 经检测，数据结构和PMIP3的数据一致

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,1999,1961,1990);
end

%% 读取Goosse2012同化的数据，已经是年分辨率
function tasField=readGoosse2012()
% Goosse2012的数据时段：1000-2000，年分辨率，已经做了距平
% 读取原始数据，并进行空间规则化转化
tas_yrMean=subf_readNC('./PMIP3_en/ts_Ayr_2deg_assim_loveclim_mann2009_1000_2000_NH.nc','ts');
A=size(tas_yrMean);

% figure
% imagesc(tas_yrMean(:,:,1));
% 经检测，数据结构和PMIP3的数据一致

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
field=geoFormat;
tasField=subf_yrAnomaly(field,1000,1999,1961,1990);
end

%%
% ******************计算极地放大指数**************************
% 计算极地放大指数，用30年的数据来计算，30年是一个气候基准期
% 1000-1040; 1041-1070;1071-1100,...,
% 注意：因为有1001个数据，所以分配方式：41*1+30*32，总计33个数据
% ***********************************************************
%
% 计算模式的AAI，适用于PDA,BCC,CCSM4,FGOALg2,GISS,IPSL,MPI,CSIRO等模式
function [AAyearValue,trendValue,AA30yrValue]=computeAAI(tasField)
% 输入数据起止时间：1000-2000年
A=size(tasField);
regSeries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    regSeries(y)=mean2(tasField(:,:,y));
end
% 60-90°N，30-60°N，0-30°N，0-90°N四个纬度带
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % 每30°一个纬度带，总计3个纬度带
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end

AAyearValue=zeros(1,1000);  % 1001-2000年
% 1000-1040的系数单独计算
z=zonalSeries(1,:)'; g=regSeries';
for i=2:1001
   AAyearValue(i-1)= abs(z(i)-z(i-1))/abs(g(i)-g(i-1));
end
%
% 对曲线进行线性拟合
tt=1:1000;
[p,S]=polyfit(tt',AAyearValue',1);  % 必须是列向量
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),1000);

AAyearValue=smooth(AAyearValue,31);

plot(tt,AAyearValue,'k',tt,trendValue,'b','linewidth',1.5);

end

%% 计算模式的AAI，适用于FGOALS-g1模式
function [AAyearValue,trendValue,AA30yrValue]=computeAAIFGOALSg1(tasField)
% FGOALS-g1模式模拟结果时间长度：1000-1999年
A=size(tasField);
regSeries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    regSeries(y)=mean2(tasField(:,:,y));
end
% 60-90°N，30-60°N，0-30°N，0-90°N四个纬度带
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % 每30°一个纬度带，总计3个纬度带
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end

AA30yrValue=zeros(33,1);
AAyearValue=zeros(1,A(1,3)+1); % 因为缺少2000年的数据，为了保证数据长度和上面一致，加1
% 1000-1040的系数单独计算
z=zonalSeries(1,1:41)'; x=regSeries(1:41)';
[b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
AA30yrValue(1)=b(2);
AAyearValue(1:41)=AA30yrValue(1);
% 从1041-1970的系数用下面的循环计算
for i=1:31
    head=(i-1)*30+41+1; tail=i*30+41; % 每30年的数据做一次回归
    z=zonalSeries(1,head:tail)'; x=regSeries(head:tail)';
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    AA30yrValue(i+1)=b(2);
    AAyearValue(head:tail)=AA30yrValue(i+1);
end
% 1971-1999的系数单独计算
z=zonalSeries(1,(1971-1000+1):end)'; x=regSeries((1971-1000+1):end)';
[b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
AA30yrValue(33)=b(2);
AAyearValue((1971-1000+1):end)=AA30yrValue(33);
%
% 对曲线进行线性拟合
tt=1:length(AA30yrValue);
[p,S]=polyfit(tt',AA30yrValue,1);  % 必须是列向量
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),A(1,3)+1);
end

%% 计算模式的AAI，适用于HadCM3模式
function [AAyearValue,trendValue,AA30yrValue]=computeAAIHadCM3(tasField)
% HadCM3数据分为两段：185912-200512；085001-185012，缺少1851-1859的数据
% ！！！缺少的1851-1859年这9年的数据暂时不管，用0代替
% ！！！计算AAI的时候，1851-1880年的AAI用1860-1880年的数据代替
A=size(tasField);
regSeries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    regSeries(y)=mean2(tasField(:,:,y));
end
% 60-90°N，30-60°N，0-30°N，0-90°N四个纬度带
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % 每30°一个纬度带，总计3个纬度带
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end

AA30yrValue=zeros(33,1);
AAyearValue=zeros(1,A(1,3));
% --(1)--
% 1000-1040的系数单独计算，数据索引起止：1-41
z=zonalSeries(1,1:41)'; x=regSeries(1:41)';
[b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
AA30yrValue(1)=b(2);
AAyearValue(1:41)=AA30yrValue(1);
% --(2)--
% 从1041-1850的系数用下面的循环计算，数据索引起止：42-851
for i=1:27
    head=(i-1)*30+41+1; tail=i*30+41; % 每30年的数据做一次回归
    z=zonalSeries(1,head:tail)'; x=regSeries(head:tail)';
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    AA30yrValue(i+1)=b(2);
    AAyearValue(head:tail)=AA30yrValue(i+1);
end
% --(3)--
% 1851-1880年（数据索引起止：852-881）的系数单独计算，由于1851-1859年（数据索引起止：852-860）
% 因此用1860-1880年间的数据代替，数据索引起止：861-881
z=zonalSeries(1,(1860-1000+1):(1880-1000+1))'; x=regSeries((1860-1000+1):(1880-1000+1))';
[b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
AA30yrValue(29)=b(2);
AAyearValue((1851-1000+1):(1880-1000+1))=AA30yrValue(29);
% --(4)--
% 1881-2000年的系数单独计算，数据索引起止：882-1001
for i=29:32
    head=(i-1)*30+41+1; tail=i*30+41; % 每30年的数据做一次回归
    z=zonalSeries(1,head:tail)'; x=regSeries(head:tail)';
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    AA30yrValue(i+1)=b(2);
    AAyearValue(head:tail)=AA30yrValue(i+1);
end
%
% 对曲线进行线性拟合
tt=1:length(AA30yrValue);
[p,S]=polyfit(tt',AA30yrValue,1);  % 必须是列向量
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),A(1,3));
end

%% 计算20CR-V2c数据的AAI，20CR-V2c数据：1851-2000
function AAyearValue=computeAAI20CR(field20CRV2c)

D=size(field20CRV2c);
regSeries=zeros(1,D(1,3));
for y=1:D(1,3)    % y: year
    regSeries(y)=mean2(field20CRV2c(:,:,y));
end
% 60-90°N，30-60°N，0-30°N，0-90°N四个纬度带
zonalSeries=zeros(D(1,1)/15,D(1,3));  % 每30°一个纬度带，总计3个纬度带
for i=1:D(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(field20CRV2c(head:tail,:,i));
    end
end
AA30yrValue=zeros(5,1);
AAyearValue=zeros(1,(2000-1000+1));
% 1000-1850的全部是空值
AAyearValue(1:851)=NaN;
% 从1851-2000的系数用下面的循环计算
for i=1:5
    head=(i-1)*30+1; tail=i*30; % 每30年的数据做一次回归
    z=zonalSeries(1,head:tail)'; x=regSeries(head:tail)';
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    AA30yrValue(i)=b(2);
    AAyearValue((head+851):(tail+851))=AA30yrValue(i);
end
disp(mean(AA30yrValue));
end

function plotAAI(t,AAI,trendValue,color,slope,pvalue)
plot(t,AAI,color,t,trendValue,'k--','linewidth',1.5);
name=strcat('trend =',32,num2str(roundn(slope,-3)),'/30 yrs;',32,'p =',32, num2str(roundn(pvalue,-3)));
xlabel('Year'); ylabel({'AA index'}); title(name);
set(gca,'ylim',[1 3],'fontsize',12,'fontname','cambria');
hold on
end
