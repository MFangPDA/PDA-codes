clc
clear all
close all

%����һ��ģ�����ݼ���������ƽ��ֵ�� sdata��ʾ��Ʊ��ÿ�ռ۸�仯��
t = 0:300;
dailyFluct = gallery('normaldata',size(t),2);
sdata = cumsum(dailyFluct) + 20 + t/100;
%�����ֵ
mean(sdata)

figure
plot(t,sdata);
legend('Original Data','Location','northwest');
xlabel('Time (days)');
ylabel('Stock Price (dollars)');

%����ȥ�������ݣ����Ҵ�ԭʼ�������Ƴ�
detrend_sdata = detrend(sdata);
trend = sdata - detrend_sdata;
mean(detrend_sdata)

hold on
plot(t,trend,':r')
plot(t,detrend_sdata,'m')
plot(t,zeros(size(t)),':k')
legend('Original Data','Trend','Detrended Data',...
       'Mean of Detrended Data','Location','northwest')
xlabel('Time (days)');
ylabel('Stock Price (dollars)');
% ��������������������������������
% ��Ȩ����������ΪCSDN��������С��������ԭ�����£���ѭ CC 4.0 BY-SA ��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
% ԭ�����ӣ�https://blog.csdn.net/wokaowokaowokao12345/article/details/60138308