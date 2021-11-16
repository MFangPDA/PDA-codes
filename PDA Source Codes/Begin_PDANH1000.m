% ***********************************************************************************
% 说明：本代码是基于Offline-EnSRF同化策略的同化代码
%    (1):  背景场抽取原始资料：MPI的月平均气温（时间段：085101-185012）
%    (2):  背景场抽取方式：随机从0851-1850年的年平均气温中抽取100个作为背景场集合，
%          即集合成员个数为100个。然后，这100个集合就固定下来，所有的实验、所有的
%          年份都用这100个集合作为背景场集合
%    (3):  同化结果验证资料：观测资料；20CR-V2c资料；以及代用资料重建的序列
%    (4):  同化过程中的HX是提前已经模拟好的，时间范围为0851-1850年对应代用资料点的模拟
%    (5)   研究区范围：
%               90
%          -180    180
%               00
%    (6): 所有NC数据在被使用之前，都需要先转换成与真实地理空间一致的结构
% 注意：Hakim et al.,所用的背景场生成方案
% we use the MPI last millennium (MPI-LM) simulation, which is a coupled
% atmosphere-ocean-sea ice simulation covering the 850–1850 C.E. time period.
% Monthly data are averaged to a calendar year,and the temporal mean over the entire data set is removed.
% 是相对于整个时间段的距平值
% **********************************************************************************

clc
clear

% 调用子文件夹Tool中的m文件
addpath(genpath(pwd));

%% -------------------声明全局变量----------------------
% 定义同化的起止时间和同化时间长度
global startDA endDA lenDA;
% MPI数据的相关全局变量
global startMPI endMPI lenMPI pathMPI;
% 代用资料数据起止年份及保存位置
global startProxy endProxy pathProxy;
% 集合个数和MCMC次数
global N MCMC;
% 区域内网格数量
global nGrid nCols nRows;
% 随机数生成方式
global randFlag;  % =0：表示使用已有的随机数; =1：表示重新生成随机数
% 同化的代用资料
global proxyFlag; % =0：表示同化所有代用资料；=1：随机挑选75%的资料; 

% 全局变量赋值(固定值)
startMPI=851;  endMPI=1849;
lenMPI=endMPI-startMPI+1;
startProxy=1000; endProxy=2000;
pathMPI='./Data/tas_Amon_2deg_MPI_085101-184912_NH.nc';
pathProxy='./Data/NHmultiProxies.mat';
% 
% 全局变量赋值(动态值)
startDA=1000; endDA=2000;
lenDA=endDA-startDA+1;
% 
% 背景场集合=100；这种大集合数目的情况下不用做协方差局地化（Setiger et al., 2017）
N=100;    
% 按照Hakim et al.,2016和Dee et al.,2016的说法，MCMC是指随机采样观测资料的次数
% 即：总共进行50次同化，每次同化过程中的观测资料不一样，而背景场始终保持不变
MCMC=50;
% 
randFlag=1;
proxyFlag=1;

%% 第一部分：同化开始之前的准备工作------------------

% 注意：所有读取的NC数据在被使用之前，都需要先转换
% 成与真实地理空间一致的结构，然后再进行相关操作
% -----1.1 读取MPI数据，提取空间网格数量参数-------
tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);clear('tas_yrMean');
% 计算MPI三维数组的结构，等同于区域三维结构
% 对部分全局变量赋值
nCols=A(1,1);  nRows=A(1,2);
nGrid=A(1,1)*A(1,2);

% -----生成从MPI采集背景场集合时需要的随机数-------
if (randFlag==0)  % 使用已有的随机数
    Num=load('./optIndex/optNum.mat');
    randEnsemble=Num.randEnsemble;
    N=length(randEnsemble);
else              % 重新生成从MPI采集背景场集合时需要的随机数
    % 生成不重复的随机数
    randEnsemble=randperm(lenMPI,N);
end
% -----------------------------------------------------

%% 第二部分： 同化试验----------------------------------
% 定义存放MCMC结果的矩阵
% 只保存集合平均值，用于计算时间序列与空间场
Xa_2D=zeros(nGrid,lenDA);    
% 声明如此大的数组，会超过内存的大小，因此必须简化
% 用于计算不确定性区间，即所有成员的结果都保存
% Xa_3D=zeros(nGrid,N*MCMC,lenDA); 
Xa_reg=zeros(lenDA,N*MCMC);  % 只保留每个集合成员的区域平均值
for MC=1:MCMC
    Str=strcat('>> LUM is the ',32, num2str(MC),'th MCMC assimilating process<<');
    fprintf('%12s \n',Str);
    
    % randEnsemble=randperm(lenMPI,N);
    [Xa,Xa_ens]=mainf_PDAProcess(randEnsemble,'LUM');

    Xa_2D=Xa_2D+Xa;
    
    head=(MC-1)*N+1; tail=MC*N;
    Xa_reg(:,head:tail)=Xa_ens;
end
Xa_2D=Xa_2D/MCMC;

% (a)把同化结果从列向量的形式转化成2维空间格式
geoFormat=zeros(nRows,nCols,lenDA);
for i=1:lenDA
    Temp=reshape(Xa_2D(:,i),nRows,nCols);
    geoFormat(:,:,i)=Temp;
end
PDATas=geoFormat;
save('./Result/PDATas.mat','PDATas');

% (b)把通过Xa_3D求同化结果的标准差
% Xa_reg的结构：Xa_reg=zeros(lenDA,N*MC_iter);
stdXa=zeros(1,lenDA);
for i=1:lenDA
    ensXa=Xa_reg(i,:);
    stdXa(i)=std(ensXa);
end
save('./Result/stdXa.mat','stdXa');
