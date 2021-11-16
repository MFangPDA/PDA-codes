% ****************************************************
% ����˵����ѡ����ͬ�������н���ͬ���Ĵ�������
%
%
% ****************************************************

function mainf_selProxies(proxyFlag)

% ֱ�Ӵ�screenProxy.mat����ѡ����ͬ��������
Data=load('./Data/screenProxy.mat');
Parameters=Data.screenProxy;

% ����ѡ�񷽰���
%��1�������д������϶�ͬ��
if (proxyFlag==0)
    save('./Data/Parameters.mat','Parameters');
%��2������Hakim et al.,2016��Dee et al.,2016,���ѡ��75%������
else
    selParameters=rand75Select(Parameters, proxyFlag);
    save('./Data/selParameters.mat','selParameters');
end
end

function selParameters=rand75Select(orgParameters)
A=size(orgParameters);
% round ��������ȡ��
N=round(A(1,1)*(75/100));

% ���ɲ��ظ��������
randProxy=randperm(A(1,1),N);
selParameters=zeros(N,A(1,2));
for i=1:N
    selParameters(i,:)=orgParameters(randProxy(i),:);
end
end