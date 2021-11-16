function newAAI()

% 调用子文件夹Tool中的m文件
addpath(genpath(pwd));

% 1. Loading PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;

anomalyPDA=subf_yrAnomaly(PDATas,1000,2000,1961,1990);

PDA=anomalyPDA;

% 计算北极放大指数
computeAAI(PDA);
end

% 计算模式的AAI，适用于PDA,BCC,CCSM4,FGOALg2,GISS,IPSL,MPI,CSIRO等模式
function computeAAI(tasField)
% 输入数据起止时间：1000-2000年
A=size(tasField);
regSeries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    regSeries(y)=mean2(tasField(:,:,y));
end
% 60-90°N，30-60°N，0-30°N，0-90°N四个纬度带
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  % 每30°一个纬度带，总计3个纬度带
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
arcSeries=zonalSeries(1,:);

% 用滑动回归的方法计算AAI，滑动窗口设为31年
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


% %计算去趋势数据，并且从原始数据中移除
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
