% ****************************************************************************
% Description:
%     (1) Computing the AA index values based on different datasets,
%     including PDA, multi-models millennium simulations in PMIP3
%     (2) All of the time series are anomalous relative to their individual
%     means of 1961-1990
%     Note that: the choice of the reference time span do not affect the AAI computation
% ****************************************************************************

function AAISensitive2Window()

addpath(genpath(pwd));

data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;
anomalyPDA=subf_yrAnomaly(PDATas,1000,2000,1961,1990);
PDA=anomalyPDA;

figure
set(gcf,'unit','inches','position',[0.2,0.2,7.3,7.3]);
ha = tight_subplot(3,2,[.1 .05],[.06 .04],[.08 .04]);
%%
%
% 2.To compute AA index based on PDA reconstruction and PMIP3 simulations
% Defining a matrix to store all AA index series

intervals=[10,20,30,40,50,60];
startAAs=[1005,1010,1015,1020,1025,1030]; endAAs=[1995,1990,1985,1980,1975,1970];
for i=1:6
    interval=intervals(i);
    startAA=startAAs(i);
    endAA=endAAs(i);
    
    % 2.1 AA index derived from PDA reconstruction
    [NHseries,ARCseries,yearAAI,trendValue]=computeAAI(PDA,interval,startAA,endAA);
   
    % 2.2 To call the Mann-Kendall trend test function to check the trend of AA index
    % during the past millennium
    [slope,pValue]=MKTrendDetection(yearAAI,'');
    
    % 2.3 Poltting the millennial AA index series and its trend
    axes(ha(i));
    t=startAA:endAA;
    plotAAI(t,yearAAI,trendValue,slope,pValue,'A',strcat('Moving Window = ',num2str(intervals(i)+1),'yrs')); hold on
    xlabel('Year (AD)'); ylabel({'AA index'});
    set(gca,'ylim',[0 3.5],'fontsize',9,'fontname','Arial','linewidth',1,'Tickdir','out');
end
end


function [NHseries,ARCseries,yrAAI,trendValue]=computeAAI(tasField,intervalAA,startAA,endAA)

A=size(tasField);
NHseries=zeros(1,A(1,3));
for y=1:A(1,3)    % y: year
    NHseries(y)=mean2(tasField(:,:,y));
end
% 60-90°„N£¨30-60°„N£¨0-30°„N£¨0-90°„N
D=size(tasField);
zonalSeries=zeros(D(1,1)/15,A(1,3));  
for i=1:A(1,3)
    for j=1:D(1,1)/15
        head=(j-1)*15+1; tail=j*15;
        zonalSeries(j,i)=mean2(tasField(head:tail,:,i));
    end
end
ARCseries=zonalSeries(1,:);

lenAAI=endAA-startAA+1;
yrAAI=zeros(lenAAI,1);

for i=1:lenAAI
    head=i; tail=i+intervalAA;       % 1:31, 2:32, 3:33, ..., 971:1001
    z=ARCseries(head:tail)'; % Arctic TAS series
    x=NHseries(head:tail)';  % NH TAS series
    [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
    yrAAI(i)=b(2);
end

tt=1:length(yrAAI);
[p,S]=polyfit(tt',yrAAI,1);  
fitValue=p(1)*tt+p(2);
trendValue=linspace(fitValue(1),fitValue(end),length(yrAAI));
end

function plotAAI(t,AAI,trendValue,slope,pvalue,index,model)
plot(t,AAI,'color',[0.5451 0 0],'linewidth',1.5); hold on
plot(t,trendValue,'color',[0 0 0.5451],'linewidth',1.5);
% name=strcat('Trend =',32,num2str(roundn(slope,-3)),'/30 yrs;','%12s \n','p =',32, num2str(roundn(pvalue,-3)));

if pvalue<0.001
    name=strcat('Trend =',32,num2str(roundn(slope,-2)),'/100 yrs,','p < 0.001');
else
    name=strcat('Trend =',32,num2str(roundn(slope,-2)),'/100 yrs,','p =',32, num2str(roundn(pvalue,-3)));
end
xlabel('Year (AD)');
ylabel({'AA index'});
set(gca,'ylim',[1 3],'fontsize',9,'fontname','Arial','linewidth',1,'Tickdir','out');
text(1015,3.1,index,'fontsize',9,'fontweight','bold');
text(1135,3.1,model,'fontsize',9);
text(1020,0.7,name,'fontsize',8,'fontname','Arial','color',[0 0 0.5451]);
hold on
end
