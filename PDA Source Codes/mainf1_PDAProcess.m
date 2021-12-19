% ****************************************************************************
% Description:
%    the main process of PDA
%
% ****************************************************************************

function [Xa,Xa_ens]=mainf_PDAProcess(randEnsemble,method)

addpath(genpath(pwd));

%% Global variables
global startMPI endMPI pathMPI;
global pathProxy;
global proxyFlag; % =1£ºall proxy£»=2£º75% proxy; =3:proxy with correlation greater than a constant

%% Select proxy
mainf_selProxies(proxyFlag);
if (proxyFlag==1)
    Para=load('./Data/Parameters.mat');
    paraProxies=Para.Parameters;
else
    Para=load('./Data/selParameters.mat');
    paraProxies=Para.selParameters;
end

multiProxies=load(pathProxy);
obsProxies=multiProxies.NHProxies;   

simData=load(strcat('./HX/','simProxies',method,'.mat'));
simProxies=simData.(strcat(method,'Data'));

tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% reformation
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
% anomalies relative to the last millennium£¨Hakim et al., 2016£©
tas_yrAnomaly=subf_yrAnomaly(geoFormat,startMPI,endMPI,startMPI,endMPI);

[Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies);
end

%% main function
function [Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies)

% global variables
global startDA endDA lenDA;
global N;
global nGrid;

D=size(paraProxies);
nProxies=D(1,1);

Xb_ens=zeros(nGrid,N);   % background ensemble
Xa=zeros(nGrid,lenDA);   % analysis all
Xa_ens=zeros(lenDA,N);   % analysis in each members
Ye_ens=zeros(1,N);       % the modelled Y of the i th proxy in the year th year

for i=1:N
    vctXb=reshape(tas_yrAnomaly(:,:,randEnsemble(i)),nGrid,1);
    Xb_ens(:,i)=vctXb;
end

for year=startDA:endDA
    fprintf(strcat('****it is the',32,num2str(year),32,'year assimilating process****','\n'))
    % assimilate one-by-one
    flag=0;
    for k=1:nProxies
        yearInx=year-1000+1+3;       
        proxyIndex=paraProxies(k,1); 
        Yo=obsProxies(yearInx,proxyIndex);
        
        if (isnan(Yo))
            continue                 
        end
        flag=flag+1;
        
        for i=1:N
            Ye_ens(1,i)=simProxies(randEnsemble(i)+4,proxyIndex);
        end
        
        if (flag==1)
            PH=(Xb_ens(:,:)-repmat(mean(Xb_ens,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        else
            PH=(Xa_ens_yr(:,:)-repmat(mean(Xa_ens_yr,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        end
        PH=(PH/(N-1));    % PH----nGrid x 1
        
        % HPH is a scalar
        HPH=(Ye_ens(1,:)-mean(Ye_ens))*(Ye_ens(1,:)-mean(Ye_ens))';
        HPH=HPH/(N-1);
        
        R=simProxies(3,proxyIndex);
        K=PH/(HPH+R);
        
        % Xa_ens_yr----[nGrid x N]
        if (flag==1)
            Xa_ens_yr(:,:)=Xb_ens(:,:)+K*(Yo-Ye_ens(1,:));
        else
            Xa_ens_yr(:,:)=Xa_ens_yr(:,:)+K*(Yo-Ye_ens(1,:));
        end
    end
    
    layerIndex=year-startDA+1;
    Xa(:,layerIndex)=mean(Xa_ens_yr,2);       
    Xa_ens(layerIndex,:)=mean(Xa_ens_yr,1);   % Xa_ens_yr----nGrid*N
end
end
