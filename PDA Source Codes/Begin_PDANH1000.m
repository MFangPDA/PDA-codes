% ***********************************************************************************
% Description: 
% This is the PDA source code with Offline category from Steiger et al. 2014
%    (1) the 500 Background ensemble members are sampled from MPI-ESM-P of past1000 experiment in PMIP3
%    (2) PDA reconstruction are verified to 20CR-V2c, HadCRUT4 and proxy-based reconstructions
%    (3) the HX were modelled in advance in PSM trainning process
%    (4) study area
%               90
%          -180    180
%               00
% **********************************************************************************

clc
clear

addpath(genpath(pwd));

%% Global variables
global startDA endDA lenDA;
global startMPI endMPI lenMPI pathMPI;
global startProxy endProxy pathProxy;
global N MCMC;
global nGrid nCols nRows;

% ways of producing randoms
global randFlag;  % =0£ºexisting randoms; =1£ºnew randoms
% ways of assimilating proxy
global proxyFlag; % =0£ºall proxy£»=1£º75% proxy which are randomly selected from the all proxy 

% global values
startMPI=851;  endMPI=1849;
lenMPI=endMPI-startMPI+1;
startProxy=1000; endProxy=2000;
pathMPI='./Data/tas_Amon_2deg_MPI_085101-184912_NH.nc';
pathProxy='./Data/NHmultiProxies.mat';
startDA=1000; endDA=2000;
lenDA=endDA-startDA+1;

% we used 500-members background ensemble in the PDA. In this case, it don't need to consider the CL (Setiger et al., 2017£©
N=500;    
% 50 iterations in PDA,
MCMC=50;
randFlag=1;
proxyFlag=1;

%% Part one
% load MPI-ESM-P simulation
tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);clear('tas_yrMean');

nCols=A(1,1);  nRows=A(1,2);
nGrid=A(1,1)*A(1,2);

% randoms 
if (randFlag==0)  % existing randoms
    Num=load('./optIndex/optNum.mat');
    randEnsemble=Num.randEnsemble;
    N=length(randEnsemble);
else              % new randoms
    randEnsemble=randperm(lenMPI,N);
end


%% Part two
Xa_2D=zeros(nGrid,lenDA);    
% Xa_3D=zeros(nGrid,N*MCMC,lenDA);
% Defining such a large matrix, if we want to save all members in all MCMC, would exceed the size of memory.

% a feasible way is to save the NH series
Xa_reg=zeros(lenDA,N*MCMC);  
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

% To transform the vector to 2-D matrix
geoFormat=zeros(nRows,nCols,lenDA);
for i=1:lenDA
    Temp=reshape(Xa_2D(:,i),nRows,nCols);
    geoFormat(:,:,i)=Temp;
end
PDATas=geoFormat;
save('./Result/PDATas.mat','PDATas');

% To calculate the SD of PDA NH series
% Xa_reg=zeros(lenDA,N*MC_iter);
stdXa=zeros(1,lenDA);
for i=1:lenDA
    ensXa=Xa_reg(i,:);
    stdXa(i)=std(ensXa);
end
save('./Result/stdXa.mat','stdXa');
