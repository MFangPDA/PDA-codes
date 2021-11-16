% ***********************************************************************************
% ˵�����������ǻ���Offline-EnSRFͬ�����Ե�ͬ������
%    (1):  ��������ȡԭʼ���ϣ�MPI����ƽ�����£�ʱ��Σ�085101-185012��
%    (2):  ��������ȡ��ʽ�������0851-1850�����ƽ�������г�ȡ100����Ϊ���������ϣ�
%          �����ϳ�Ա����Ϊ100����Ȼ����100�����Ͼ͹̶����������е�ʵ�顢���е�
%          ��ݶ�����100��������Ϊ����������
%    (3):  ͬ�������֤���ϣ��۲����ϣ�20CR-V2c���ϣ��Լ����������ؽ�������
%    (4):  ͬ�������е�HX����ǰ�Ѿ�ģ��õģ�ʱ�䷶ΧΪ0851-1850���Ӧ�������ϵ��ģ��
%    (5)   �о�����Χ��
%               90
%          -180    180
%               00
%    (6): ����NC�����ڱ�ʹ��֮ǰ������Ҫ��ת��������ʵ����ռ�һ�µĽṹ
% ע�⣺Hakim et al.,���õı��������ɷ���
% we use the MPI last millennium (MPI-LM) simulation, which is a coupled
% atmosphere-ocean-sea ice simulation covering the 850�C1850 C.E. time period.
% Monthly data are averaged to a calendar year,and the temporal mean over the entire data set is removed.
% �����������ʱ��εľ�ƽֵ
% **********************************************************************************

clc
clear

% �������ļ���Tool�е�m�ļ�
addpath(genpath(pwd));

%% -------------------����ȫ�ֱ���----------------------
% ����ͬ������ֹʱ���ͬ��ʱ�䳤��
global startDA endDA lenDA;
% MPI���ݵ����ȫ�ֱ���
global startMPI endMPI lenMPI pathMPI;
% ��������������ֹ��ݼ�����λ��
global startProxy endProxy pathProxy;
% ���ϸ�����MCMC����
global N MCMC;
% ��������������
global nGrid nCols nRows;
% ��������ɷ�ʽ
global randFlag;  % =0����ʾʹ�����е������; =1����ʾ�������������
% ͬ���Ĵ�������
global proxyFlag; % =0����ʾͬ�����д������ϣ�=1�������ѡ75%������; 

% ȫ�ֱ�����ֵ(�̶�ֵ)
startMPI=851;  endMPI=1849;
lenMPI=endMPI-startMPI+1;
startProxy=1000; endProxy=2000;
pathMPI='./Data/tas_Amon_2deg_MPI_085101-184912_NH.nc';
pathProxy='./Data/NHmultiProxies.mat';
% 
% ȫ�ֱ�����ֵ(��ֵ̬)
startDA=1000; endDA=2000;
lenDA=endDA-startDA+1;
% 
% ����������=100�����ִ󼯺���Ŀ������²�����Э����ֵػ���Setiger et al., 2017��
N=100;    
% ����Hakim et al.,2016��Dee et al.,2016��˵����MCMC��ָ��������۲����ϵĴ���
% �����ܹ�����50��ͬ����ÿ��ͬ�������еĹ۲����ϲ�һ������������ʼ�ձ��ֲ���
MCMC=50;
% 
randFlag=1;
proxyFlag=1;

%% ��һ���֣�ͬ����ʼ֮ǰ��׼������------------------

% ע�⣺���ж�ȡ��NC�����ڱ�ʹ��֮ǰ������Ҫ��ת��
% ������ʵ����ռ�һ�µĽṹ��Ȼ���ٽ�����ز���
% -----1.1 ��ȡMPI���ݣ���ȡ�ռ�������������-------
tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);clear('tas_yrMean');
% ����MPI��ά����Ľṹ����ͬ��������ά�ṹ
% �Բ���ȫ�ֱ�����ֵ
nCols=A(1,1);  nRows=A(1,2);
nGrid=A(1,1)*A(1,2);

% -----���ɴ�MPI�ɼ�����������ʱ��Ҫ�������-------
if (randFlag==0)  % ʹ�����е������
    Num=load('./optIndex/optNum.mat');
    randEnsemble=Num.randEnsemble;
    N=length(randEnsemble);
else              % �������ɴ�MPI�ɼ�����������ʱ��Ҫ�������
    % ���ɲ��ظ��������
    randEnsemble=randperm(lenMPI,N);
end
% -----------------------------------------------------

%% �ڶ����֣� ͬ������----------------------------------
% ������MCMC����ľ���
% ֻ���漯��ƽ��ֵ�����ڼ���ʱ��������ռ䳡
Xa_2D=zeros(nGrid,lenDA);    
% ������˴�����飬�ᳬ���ڴ�Ĵ�С����˱����
% ���ڼ��㲻ȷ�������䣬�����г�Ա�Ľ��������
% Xa_3D=zeros(nGrid,N*MCMC,lenDA); 
Xa_reg=zeros(lenDA,N*MCMC);  % ֻ����ÿ�����ϳ�Ա������ƽ��ֵ
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

% (a)��ͬ�����������������ʽת����2ά�ռ��ʽ
geoFormat=zeros(nRows,nCols,lenDA);
for i=1:lenDA
    Temp=reshape(Xa_2D(:,i),nRows,nCols);
    geoFormat(:,:,i)=Temp;
end
PDATas=geoFormat;
save('./Result/PDATas.mat','PDATas');

% (b)��ͨ��Xa_3D��ͬ������ı�׼��
% Xa_reg�Ľṹ��Xa_reg=zeros(lenDA,N*MC_iter);
stdXa=zeros(1,lenDA);
for i=1:lenDA
    ensXa=Xa_reg(i,:);
    stdXa(i)=std(ensXa);
end
save('./Result/stdXa.mat','stdXa');
