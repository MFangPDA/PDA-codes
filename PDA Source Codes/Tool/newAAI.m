function newAAI()

% �������ļ���Tool�е�m�ļ�
addpath(genpath(pwd));

% 1. Loading PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;

anomalyPDA=subf_yrAnomaly(PDATas,1000,2000,1961,1990);

PDA=anomalyPDA;

% ���㱱���Ŵ�ָ��
computeAAI(PDA);
end

% ����ģʽ��AAI��������PDA,BCC,CCSM4,FGOALg2,GISS,IPSL,MPI,CSIRO��ģʽ
function computeAAI(tasField)
% ����������ֹʱ�䣺1000-2000��
A=size(tasField);
regSeries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    regSeries(y)=mean2(tasField(:,:,y));
end
% 60-90��N��30-60��N��0-30��N��0-90��N�ĸ�γ�ȴ�
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % ÿ30��һ��γ�ȴ����ܼ�3��γ�ȴ�
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
arcSeries=zonalSeries(1,:);

% �û����ع�ķ�������AAI������������Ϊ31��
lenAAI=(1985-1015+1);
yrAAI=zeros(1,lenAAI);

for i=1:lenAAI
    head=i; tail=i+30;
    z=zonalSeries(1,head:tail)'; x=regSeries(head:tail)';
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    yrAAI(i)=b(2);
end
t=1015:1985;
plot(t,yrAAI);

xlswrite('yrAAI.xlsx',yrAAI');

[slope,pValue]=MKTrendDetection(yrAAI,'PDA');
disp(slope);
disp(pValue);


% %����ȥ�������ݣ����Ҵ�ԭʼ�������Ƴ�
% detrend_sdata = detrend(sdata);
% Annual time-scale 
% Preindustrial Era
[Num,txt] = xlsread('./Data/Forcings/forcingsData.xlsx');

[corMatrix,pValue]=corrcoef(Num(50:851,2:end));
disp(corMatrix);  disp(pValue);

% Industrial Era
[Num,txt] = xlsread('./Data/Forcings/forcingsData.xlsx');
[corMatrix,pValue]=corrcoef(Num(852:986,2:end));
disp(corMatrix);  disp(pValue);

end
