% ****************************************************************************
% 本代码是同化的主要过程
%
%
% ****************************************************************************
function [Xa,Xa_ens]=mainf_PDAProcess(randEnsemble,method)

% 调用子文件夹中的m文件
addpath(genpath(pwd));

%% ------------------声明全局变量-----------------------
% MPI的相关全局变量
global startMPI endMPI pathMPI;
% 代用资料数据起止年份及保存位置
global pathProxy;
% 同化的代用资料选择方式
global proxyFlag; % =1：表示同化所有代用资料；=2：随机挑选75%的资料; =3:选择相关系数大于某值的资料

%% -------------------------------------------------------
% 按照一定的规则选取本轮同化过程中将要被同化的代用资料，也可以全部都同化
% 不同的实验规则，有不同的paraProxies，根据paraProxies从obsProxies和simProxies中选择数据开展同化
mainf_selProxies(proxyFlag);

if (proxyFlag==1)
    % 读取研究区内的所有代用资料
    Para=load('./Data/Parameters.mat');
    paraProxies=Para.Parameters;
else
    % 读取研究区内按照一定的规则选择出来的代用资料
    Para=load('./Data/selParameters.mat');
    paraProxies=Para.selParameters;
end

% 读取研究区内的代用资料观测值保存文件.mat
multiProxies=load(pathProxy);
% 所有将要被同化的代用资料的观测值序列，注意第一列是年份
obsProxies=multiProxies.NHProxies;   

% 导入模拟的代用资料以及相关的参数文件
simData=load(strcat('./HX/','simProxies',method,'.mat'));
% 注意：通过load导入的simProxies是一个结构体变量
% 需要通过如下方式提取其中的数据
simProxies=simData.(strcat(method,'Data'));

% 读取背景场资料
% 读取MPI数据，并进行空间规则化转化
tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
% 计算相对于整个时间段的距平值（参考Hakim et al., 2016）
tas_yrAnomaly=subf_yrAnomaly(geoFormat,startMPI,endMPI,startMPI,endMPI);

% 启动同化程序
[Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies);
end

%% 古气候数据同化的主函数
function [Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies)

% ------------------声明全局变量-----------------------
% 定义同化的起止时间和同化时间长度
global startDA endDA lenDA;
% 集合成员个数
global N;
% 区域内网格数量
global nGrid;
%
D=size(paraProxies);
nProxies=D(1,1);
%
% 接下来构造背景场集合，背景场集合的结构如下：
% 注意：按照Hakim et al., 2016的方法，在所有同化实验、
% 所有同化时间段内，每一年同化过程中都用这一个背景场集合，即：
% 从1000-2000年，每一年的背景场都一样，这样同化结果中所体现
% 出的波动信号就都来源于被同化的观测资料的信号。（很好！）
%
% 定义同化过程中的相关矩阵
Xb_ens=zeros(nGrid,N);   % 背景场集合
Xa=zeros(nGrid,lenDA);   % 只保存集合平均值，用于计算空间场与时间序列
Xa_ens=zeros(lenDA,N);   % 保存每个成员的区域均值，用于计算时间序列的不确定范围
Ye_ens=zeros(1,N);       % 第year年第i个代用的集合模拟值

for i=1:N
    vctXb=reshape(tas_yrAnomaly(:,:,randEnsemble(i)),nGrid,1);
    Xb_ens(:,i)=vctXb;
end

for year=startDA:endDA
    fprintf(strcat('****it is the',32,num2str(year),32,'year assimilating process****','\n'))
    % 采用对代用资料逐个同化的策略
    flag=0;
    for k=1:nProxies
        % obsProxies的第1行是资料类型，第2-3行是经纬度；第1列是年份值
        yearInx=year-1000+1+3;       % 代用资料的起始时间是1000年
        proxyIndex=paraProxies(k,1); % proxyIndex是这条代用资料在NHmultiProxies.mat文件中对应的列号
        Yo=obsProxies(yearInx,proxyIndex);
        
        % 注意：因为是过去1000年重建，所以有可能某些资料长度达不到，导致
        % 在前面某些年份中读出的Yo是NaN，因此需要判断一下
        if (isnan(Yo))
            continue                 % 直接进行下一条资料
        end
        % 只有第k条资料不为NaN，程序才能运行到这里，记录下非空的资料观测值的个数
        flag=flag+1;
        
        % 第1行：经度；第2行：纬度；第3行：R；第4行：这条模拟的序列其在NHmultiProxies.mat文件中的列号
        % 第5：(4+dataLen)行存放模拟值
        for i=1:N
            Ye_ens(1,i)=simProxies(randEnsemble(i)+4,proxyIndex);
        end
        
        if (flag==1)
            PH=(Xb_ens(:,:)-repmat(mean(Xb_ens,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        else
            PH=(Xa_ens_yr(:,:)-repmat(mean(Xa_ens_yr,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        end
        PH=(PH/(N-1));    % 此处的PH--> nGrid x 1
        
        % 因为所用集合数目达到了100，很大程度上不会有欠采样
        % 的问题，因此从理论上来讲可以不用再考虑局地化。（Steiger et al.,2017）
        
        % 因为是逐年逐个同化观测数据，所以此处的HPH其实是标量
        HPH=(Ye_ens(1,:)-mean(Ye_ens))*(Ye_ens(1,:)-mean(Ye_ens))';
        HPH=HPH/(N-1);
        
        % 计算卡尔曼增益；根据proxyIndex从simProxies里面提取R
        R=simProxies(3,proxyIndex);
        K=PH/(HPH+R);
        
        % Xa_ens_yr存放的是第year年的分析值的集合,结构：[nGrid x N]
        if (flag==1)
            Xa_ens_yr(:,:)=Xb_ens(:,:)+K*(Yo-Ye_ens(1,:));
        else
            Xa_ens_yr(:,:)=Xa_ens_yr(:,:)+K*(Yo-Ye_ens(1,:));
        end
    end
    
    layerIndex=year-startDA+1;
    Xa(:,layerIndex)=mean(Xa_ens_yr,2);       % mean(A,2)是矩阵求各行的均值
    Xa_ens(layerIndex,:)=mean(Xa_ens_yr,1);   % mean(A,1)是矩阵求各列的均值  % Xa_ens_yr：nGrid*N
end
end

%% -----------------------------------------------

%% 数据同化中各种中间矩阵的大小说明：
%   Xb: [M_grid, N]
%    R: [M_obs,M_obs]
%   Pb: [M_grid, M_grid]
%   HP: [M_grid,M_obs]
%  HPH: [M_obs,M_obs]
%    K =PH/(HPH+R)
%    K: [M_grid,M_obs)/([M_obs,M_obs]+[M_obs,M_obs])=[M_grid,M_obs]
%% ------------------------------------------------
