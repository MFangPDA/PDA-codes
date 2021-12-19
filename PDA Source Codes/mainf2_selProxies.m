% ****************************************************
% Description: 
%    function used to select proxy to be assimilated in this iteration
%
%
% ****************************************************

function mainf_selProxies(proxyFlag)

Data=load('./Data/screenProxy.mat');
Parameters=Data.screenProxy;

% all proxy
if (proxyFlag==0)
    save('./Data/Parameters.mat','Parameters');
% 75% proxy
else
    selParameters=rand75Select(Parameters, proxyFlag);
    save('./Data/selParameters.mat','selParameters');
end
end

function selParameters=rand75Select(orgParameters)
A=size(orgParameters);
N=round(A(1,1)*(75/100));

randProxy=randperm(A(1,1),N);
selParameters=zeros(N,A(1,2));
for i=1:N
    selParameters(i,:)=orgParameters(randProxy(i),:);
end
end