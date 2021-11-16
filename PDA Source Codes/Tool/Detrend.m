clc
clear all
close all

%创建一个模拟数据集并计算其平均值。 sdata表示股票的每日价格变化。
t = 0:300;
dailyFluct = gallery('normaldata',size(t),2);
sdata = cumsum(dailyFluct) + 20 + t/100;
%计算均值
mean(sdata)

figure
plot(t,sdata);
legend('Original Data','Location','northwest');
xlabel('Time (days)');
ylabel('Stock Price (dollars)');

%计算去趋势数据，并且从原始数据中移除
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
% ————————————————
% 版权声明：本文为CSDN博主「吕小猪不坏」的原创文章，遵循 CC 4.0 BY-SA 版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/wokaowokaowokao12345/article/details/60138308