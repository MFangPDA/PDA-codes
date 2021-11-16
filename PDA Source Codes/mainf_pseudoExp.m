% 输入同化结果到正向模型中，模拟树轮宽度，再和观测的树轮宽度值进行比较
% 检验同化结果的可信性，分为：1850-2000；1000-1849两个时间段，同时把
% MPI-ESM-P也加进来，互相比较，然后做成箱线图
function postf6_pseudoExp()

% 调用子文件夹中的m文件
addpath(genpath(pwd));

%% ------------------声明全局变量-----------------------
% GISTEMP数据起止年份及保存位置
global startGISTEMP   endGISS   pathGISS;
% MPI数据的起止年份及保存位置
global startMPI endMPI pathMPI;
% 代用资料数据起止年份及保存位置
global startProxy endProxy pathProxy;
% 模型训练起止年份及训练的时间长度
global startTrain endTrain lenTrain;
% 正向模型模拟的起止年份及模拟的时间长度
global startSim endSim lenSim;
% 训练数据网格空间分辨率
global interval;
% 定义全局变量
global startDA  endDA  lenDA;
%
global smoothSpan;
%
% 全局变量赋值(固定值)
startProxy=1000; endProxy=2000;
startGISTEMP=1880;  endGISS=2017;
startTrain=1880; endTrain=2000;
lenTrain=endTrain-startTrain+1;
startSim=startMPI;
endSim=endMPI;
lenSim=endSim-startSim+1;
interval=2;
startDA=1000;   endDA=2000;   lenDA=endDA-startDA+1;
smoothSpan=1;
pathGISS='./Data/air.2x2.250.mon.anom.land_NH.nc';
pathProxy='./Data/NHmultiProxies.mat';

% Step 1: 读取GISTEMP数据并规则化
tas_monMean=subf_readNC(pathGISS,'air');
tas_yrMean=subf_yrMean(tas_monMean);    % 距平处理放在下面程序中
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=tas_yrMean(:,:,i)';  % 注意：GISTEMP数据规则化方法和其他数据不一样
end
normGISS=geoFormat;
% 检查地理空间规则化效果
% plotT(normGISS(:,:,end));

% Step 2: 加载同化结果并计算距平值
PDA=load('./Result/PDAResult/PDATas.mat');
PDATas=PDA.ANNTas;
% All time series used here are anomalous relative to the means of 1000-2000,
Anomaly=subf_yrAnomaly(PDATas,startDA,endDA,1880,2000);

% 然后再对距平值进行平滑处理
Smooth=zeros(size(Anomaly));  A=size(Anomaly);
for i=1:A(1,1)
    for j=1:A(1,2)
        Smooth(i,j,:)=smooth(Anomaly(i,j,:),smoothSpan,'lowess');
    end
end
PDA=Smooth;

% Step 3: 读取MPI-ESM-P的过去千年（1000-2000）模拟数据，并计算距平值（相对1000-2000）
tasFieldMPI=readPMIP3('MPI');

% Step 4: 读取研究区内的所有代用资料
multiProxies=load(pathProxy);
% 所有将要被同化的代用资料的观测值序列，注意第一列是年份;第1行是资料类型、第2-3行是经纬度
obsProxies=multiProxies.NHProxies;
A=size(obsProxies);
allProxies=A(1,2);        % 记住：第1列是年份
screenProxy=[];
corMPIvct=[]; ceMPIvct=[];
corPDAvct=[]; cePDAvct=[];

% % Step 5: 定义存放正向模型所模拟的851-1849年的树轮、R以及经纬度的矩阵
% % 第1行：经度；第2行：纬度；第3行：R；第4行：这条模拟的序列其在NHmultiProxies.mat文件中的列号
% % 第5：(4+dataLen)行存放模拟值
% dataLen=endSim-startSim+1;
% LUMData=zeros(4+dataLen,allProxies);    LUMData(1:4,1)=9999;
% ANNData=zeros(4+dataLen,allProxies);    ANNData(1:4,1)=9999;
% for i=5:(4+dataLen)
%     LUMData(i,1)=850+(i-4);
% end
% RMS=zeros(allProxies-1,2);
% Coreff=zeros(allProxies-1,2);

% Step 5: 进入训练的循环过程
for i=2:allProxies     % 记住：第1列是年份
    Index=i;           % index:代用资料在NHmultiProxies.mat文件中的列号
    Proxy=obsProxies(:,Index);
    pxyType=Proxy(1);
    lon=Proxy(2);
    lat=Proxy(3);
    [row,col] = subf_grdNO(lat,lon);
    
    % 挑选边界之内的树轮点，边界之外的一律不要
    if (-180<=lon && lon<=180)
        if (0<=lat && lat<=90)
            
            % 注意：代用资料、GISSTEMP在1880-2000年间的数据，各自都存在空值的弊端
            % 因此，模型训练时，代用资料、GISSTEMP都截取1880-2000年的数据，并判断
            % 然后用两组数据的非空值的交集进行训练，得到训练的模型
            
            % obsProxy：第1行是资料类型、第2-3行是经纬度
            obsProxy=Proxy((startTrain-startProxy+1+3):(endTrain-startProxy+1+3));
            Tas=normGISS(row,col,(startTrain-startGISTEMP+1):(endTrain-startGISTEMP+1));
            obsTas=reshape(Tas,lenTrain,1);
            
            flagGISS=findGISS(obsTas);
            flagProxy=findProxy(obsProxy);
            
            % 提取两条序列非空值的交集
            intersection=findIntersection(flagGISS,flagProxy);
            
            z=obsProxy(intersection);
            x=obsTas(intersection)-mean(obsTas(intersection));  % 记住气温求全时段距平值
            
            if (length(z)>=30 && length(x)>=30)   % 只有交集数据的长度大于30，回归才有意义
                C=corrcoef(z,x);
                screenProxy=[screenProxy;Index,pxyType,lon,lat,C(1,2)];
                
                %  Step 4: 在此处直接添加模型训练的程序效果更好
                fprintf(strcat('The',32,num2str(i),'th column proxy is calibrated now!','\n'))
                % 注：32是空格的ASCII码；strcat会自动删除被链接字符串的头尾的空格
                
                % 训练LUM模型(一元线性回归): Y=b0+b1*X
                [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
                % 得到模型训练的结果与参数  R1=stats(1,4);
                b0=b(1);           % b0
                b1=b(2);           % b1
                
                % 下面读入同化结果和模拟结果，分别进行模拟
                % (1)读取PDA对应网格上的气温值
                tasPDAgrid=PDA(row,col,:);
                % 需要把数据转换成一维的形式
                tasPDA=reshape(tasPDAgrid,lenSim,1);
                % (2)模拟的代用资料序列（1000-2000年之间的）
                proxyPDA=b0+b1*tasPDA;
                
                % (3)读取MPI对应网格上的气温值
                tasMPIgrid=tasFieldMPI(row,col,:);
                % 需要把数据转换成一维的形式
                tasMPI=reshape(tasMPIgrid,lenSim,1);
                % (4)模拟的代用资料序列（1000-2000年之间的）
                proxyMPI=b0+b1*tasMPI;
                
                % 在这里预留位置，可以添加更多的模式结果看看
                
                % (5)读取代用资料序列。注意：代用资料序列的值的时间范围
                % 有可能是1000-2000之间的某段时间，所以需要判断
                obsProxy=Proxy(4:end);
                flagProxy=findProxy(obsProxy);
                flagMPI=ones(1,lenDA);
                flagPDA=flagMPI;
                
                % (6)查找交集，并找到验证数据序列
                overlap=findIntersection(flagProxy,flagPDA);
                
                verifyProxy=obsProxy(overlap);
                verifyMPI=proxyMPI(overlap);
                verifyPDA=proxyPDA(overlap);
                
                % (7)计算相关系数和CE
                [corMPI,ceMPI]=verifyIDX(verifyMPI,verifyProxy);
                [corPDA,cePDA]=verifyIDX(verifyPDA,verifyProxy);
                
                % (8)暂存验证指标
                corMPIvct=[corMPIvct;corMPI];  ceMPIvct=[ceMPIvct;ceMPI];
                corPDAvct=[corPDAvct;corPDA];  cePDAvct=[cePDAvct;cePDA];
            end
        end
    end
end
save('./Result/VerResult/corMPIvct.mat','corMPIvct');
save('./Result/VerResult/corPDAvct.mat','corPDAvct');
save('./Result/VerResult/ceMPIvct.mat','ceMPIvct');
save('./Result/VerResult/cePDAvct.mat','cePDAvct');

% 作图，直方图
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
subplot(2,2,1);
plotHistCOR(corMPIvct,'Correlation coefficient','b'); hold on
plotHistCOR(corPDAvct,'Correlation coefficient','r'); hold on
subplot(2,2,2);
plotHistCE(ceMPIvct,'Coefficient of efficiency','b'); hold on
plotHistCE(cePDAvct,'Coefficient of efficiency','r'); hold on
subplot(2,2,3);
ceChange=corPDAvct-corMPIvct;
plotHistCOR(ceChange,'Change in correlation coefficient','g'); hold on
subplot(2,2,4);
ceChange=cePDAvct-ceMPIvct;
plotHistCE(ceChange,'Change in coefficient of efficiency','g'); hold on
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
text(-0.9,80,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',14,'fontweight','bold'); 
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
text(-4.5,160,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',12,'fontweight','bold'); 
% set(gcf,'color','none');   % 图形背景设为无色
set(gca,'color','none');     % 坐标轴背景为无色，这条很重要，通常图形背景的白色（实际为坐标轴默认背景色），设置为透明，就可以进行多幅图的叠加
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end

% 计算相关系数Cor和CE
function [Cor,CE]=verifyIDX(rec,ref)
Correlation= corrcoef(rec,ref); % corrcoef计算的是输入数据每列之间的相关系数
Cor=Correlation(1,2);
CE=1-(sum((ref-rec).*(ref-rec)))/(sum((ref-mean(ref)).*(ref-mean(ref))));
end

function [row,col] = subf_grdNO(lat,lon)
% 重新编写的根据坐标点(lat,lon)查找其在真实地理空间中的网格号的方法
% 所有NC数据在被提取时间序列之前，都需要先转换成与真实地理空间一致的结构

% 例如：坐标点(89,-179)的点，应该是真实地理空间
% 中左上角的那一个网格点，对应于geoFormation
% 矩阵中，对应的就是geoFormation(1,1)的数据
% (89,-179) --> (1,1)       (89,179) --> (1,180)
% (31,-179)  --> (30,1)    (31,179) -->  (30,180)

interval=2;

% ceil:向上取整
row=ceil((90-lat)/interval);
col=ceil((180+lon)/interval);
end

function flagGISS=findGISS(Series)
len=length(Series);
flagGISS=zeros(1,len);
for i=1:len
    if(Series(i)<=15)
        flagGISS(i)=1;  % 非空值的话，标记为1
    else
        flagGISS(i)=0;
    end
end
end

function flagProxy=findProxy(Series)
len=length(Series);
flagProxy=zeros(1,len);
for i=1:len
    if(~isnan(Series(i)))
        flagProxy(i)=1;
    else
        flagProxy(i)=0;
    end
end
end

function intersection=findIntersection(flag1,flag2)
len=length(flag1);
intersection=[];
for i=1:len
    if (flag1(i)==1 && flag2(i)==1)  % 两个数据都非空
        intersection=[intersection,i];
    end
end
end

%% 读取BCC,CCSM4,FGOALg2,GISS,IPSL,MPI的通用程序
function tasField=readPMIP3(modelname)
% 以上数据时段：085001-200512
% 读取原始数据，并进行空间规则化转化
tas_Amon=subf_readNC(strcat('./PMIP3_en/tas_Amon_2deg_',modelname,'_085001_200512_NH.nc'),'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% 只提取1000-2000的数据
field=geoFormat(:,:,(1000-850+1):(2000-850+1));
tasField=subf_yrAnomaly(field,1000,2000,1880,2000);
end
