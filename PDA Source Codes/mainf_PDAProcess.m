% ****************************************************************************
% ��������ͬ������Ҫ����
%
%
% ****************************************************************************
function [Xa,Xa_ens]=mainf_PDAProcess(randEnsemble,method)

% �������ļ����е�m�ļ�
addpath(genpath(pwd));

%% ------------------����ȫ�ֱ���-----------------------
% MPI�����ȫ�ֱ���
global startMPI endMPI pathMPI;
% ��������������ֹ��ݼ�����λ��
global pathProxy;
% ͬ���Ĵ�������ѡ��ʽ
global proxyFlag; % =1����ʾͬ�����д������ϣ�=2�������ѡ75%������; =3:ѡ�����ϵ������ĳֵ������

%% -------------------------------------------------------
% ����һ���Ĺ���ѡȡ����ͬ�������н�Ҫ��ͬ���Ĵ������ϣ�Ҳ����ȫ����ͬ��
% ��ͬ��ʵ������в�ͬ��paraProxies������paraProxies��obsProxies��simProxies��ѡ�����ݿ�չͬ��
mainf_selProxies(proxyFlag);

if (proxyFlag==1)
    % ��ȡ�о����ڵ����д�������
    Para=load('./Data/Parameters.mat');
    paraProxies=Para.Parameters;
else
    % ��ȡ�о����ڰ���һ���Ĺ���ѡ������Ĵ�������
    Para=load('./Data/selParameters.mat');
    paraProxies=Para.selParameters;
end

% ��ȡ�о����ڵĴ������Ϲ۲�ֵ�����ļ�.mat
multiProxies=load(pathProxy);
% ���н�Ҫ��ͬ���Ĵ������ϵĹ۲�ֵ���У�ע���һ�������
obsProxies=multiProxies.NHProxies;   

% ����ģ��Ĵ��������Լ���صĲ����ļ�
simData=load(strcat('./HX/','simProxies',method,'.mat'));
% ע�⣺ͨ��load�����simProxies��һ���ṹ�����
% ��Ҫͨ�����·�ʽ��ȡ���е�����
simProxies=simData.(strcat(method,'Data'));

% ��ȡ����������
% ��ȡMPI���ݣ������пռ����ת��
tas_Amon=subf_readNC(pathMPI,'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
% �������������ʱ��εľ�ƽֵ���ο�Hakim et al., 2016��
tas_yrAnomaly=subf_yrAnomaly(geoFormat,startMPI,endMPI,startMPI,endMPI);

% ����ͬ������
[Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies);
end

%% ����������ͬ����������
function [Xa,Xa_ens]=PDA(randEnsemble,tas_yrAnomaly,obsProxies,paraProxies,simProxies)

% ------------------����ȫ�ֱ���-----------------------
% ����ͬ������ֹʱ���ͬ��ʱ�䳤��
global startDA endDA lenDA;
% ���ϳ�Ա����
global N;
% ��������������
global nGrid;
%
D=size(paraProxies);
nProxies=D(1,1);
%
% ���������챳�������ϣ����������ϵĽṹ���£�
% ע�⣺����Hakim et al., 2016�ķ�����������ͬ��ʵ�顢
% ����ͬ��ʱ����ڣ�ÿһ��ͬ�������ж�����һ�����������ϣ�����
% ��1000-2000�꣬ÿһ��ı�������һ��������ͬ�������������
% ���Ĳ����źžͶ���Դ�ڱ�ͬ���Ĺ۲����ϵ��źš����ܺã���
%
% ����ͬ�������е���ؾ���
Xb_ens=zeros(nGrid,N);   % ����������
Xa=zeros(nGrid,lenDA);   % ֻ���漯��ƽ��ֵ�����ڼ���ռ䳡��ʱ������
Xa_ens=zeros(lenDA,N);   % ����ÿ����Ա�������ֵ�����ڼ���ʱ�����еĲ�ȷ����Χ
Ye_ens=zeros(1,N);       % ��year���i�����õļ���ģ��ֵ

for i=1:N
    vctXb=reshape(tas_yrAnomaly(:,:,randEnsemble(i)),nGrid,1);
    Xb_ens(:,i)=vctXb;
end

for year=startDA:endDA
    fprintf(strcat('****it is the',32,num2str(year),32,'year assimilating process****','\n'))
    % ���öԴ����������ͬ���Ĳ���
    flag=0;
    for k=1:nProxies
        % obsProxies�ĵ�1�����������ͣ���2-3���Ǿ�γ�ȣ���1�������ֵ
        yearInx=year-1000+1+3;       % �������ϵ���ʼʱ����1000��
        proxyIndex=paraProxies(k,1); % proxyIndex����������������NHmultiProxies.mat�ļ��ж�Ӧ���к�
        Yo=obsProxies(yearInx,proxyIndex);
        
        % ע�⣺��Ϊ�ǹ�ȥ1000���ؽ��������п���ĳЩ���ϳ��ȴﲻ��������
        % ��ǰ��ĳЩ����ж�����Yo��NaN�������Ҫ�ж�һ��
        if (isnan(Yo))
            continue                 % ֱ�ӽ�����һ������
        end
        % ֻ�е�k�����ϲ�ΪNaN������������е������¼�·ǿյ����Ϲ۲�ֵ�ĸ���
        flag=flag+1;
        
        % ��1�У����ȣ���2�У�γ�ȣ���3�У�R����4�У�����ģ�����������NHmultiProxies.mat�ļ��е��к�
        % ��5��(4+dataLen)�д��ģ��ֵ
        for i=1:N
            Ye_ens(1,i)=simProxies(randEnsemble(i)+4,proxyIndex);
        end
        
        if (flag==1)
            PH=(Xb_ens(:,:)-repmat(mean(Xb_ens,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        else
            PH=(Xa_ens_yr(:,:)-repmat(mean(Xa_ens_yr,2),1,N))*(Ye_ens(1,:)-mean(Ye_ens))';
        end
        PH=(PH/(N-1));    % �˴���PH--> nGrid x 1
        
        % ��Ϊ���ü�����Ŀ�ﵽ��100���ܴ�̶��ϲ�����Ƿ����
        % �����⣬��˴��������������Բ����ٿ��Ǿֵػ�����Steiger et al.,2017��
        
        % ��Ϊ���������ͬ���۲����ݣ����Դ˴���HPH��ʵ�Ǳ���
        HPH=(Ye_ens(1,:)-mean(Ye_ens))*(Ye_ens(1,:)-mean(Ye_ens))';
        HPH=HPH/(N-1);
        
        % ���㿨�������棻����proxyIndex��simProxies������ȡR
        R=simProxies(3,proxyIndex);
        K=PH/(HPH+R);
        
        % Xa_ens_yr��ŵ��ǵ�year��ķ���ֵ�ļ���,�ṹ��[nGrid x N]
        if (flag==1)
            Xa_ens_yr(:,:)=Xb_ens(:,:)+K*(Yo-Ye_ens(1,:));
        else
            Xa_ens_yr(:,:)=Xa_ens_yr(:,:)+K*(Yo-Ye_ens(1,:));
        end
    end
    
    layerIndex=year-startDA+1;
    Xa(:,layerIndex)=mean(Xa_ens_yr,2);       % mean(A,2)�Ǿ�������еľ�ֵ
    Xa_ens(layerIndex,:)=mean(Xa_ens_yr,1);   % mean(A,1)�Ǿ�������еľ�ֵ  % Xa_ens_yr��nGrid*N
end
end

%% -----------------------------------------------

%% ����ͬ���и����м����Ĵ�С˵����
%   Xb: [M_grid, N]
%    R: [M_obs,M_obs]
%   Pb: [M_grid, M_grid]
%   HP: [M_grid,M_obs]
%  HPH: [M_obs,M_obs]
%    K =PH/(HPH+R)
%    K: [M_grid,M_obs)/([M_obs,M_obs]+[M_obs,M_obs])=[M_grid,M_obs]
%% ------------------------------------------------
