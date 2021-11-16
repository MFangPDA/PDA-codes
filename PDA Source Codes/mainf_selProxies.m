% ****************************************************
% 函数说明：选择本轮同化过程中将被同化的代用资料
%
%
% ****************************************************

function mainf_selProxies(proxyFlag)

% 直接从screenProxy.mat里面选择本轮同化的资料
Data=load('./Data/screenProxy.mat');
Parameters=Data.screenProxy;

% 资料选择方案：
%（1）把所有代用资料都同化
if (proxyFlag==0)
    save('./Data/Parameters.mat','Parameters');
%（2）参照Hakim et al.,2016和Dee et al.,2016,随机选择75%的资料
else
    selParameters=rand75Select(Parameters, proxyFlag);
    save('./Data/selParameters.mat','selParameters');
end
end

function selParameters=rand75Select(orgParameters)
A=size(orgParameters);
% round 四舍五入取整
N=round(A(1,1)*(75/100));

% 生成不重复的随机数
randProxy=randperm(A(1,1),N);
selParameters=zeros(N,A(1,2));
for i=1:N
    selParameters(i,:)=orgParameters(randProxy(i),:);
end
end