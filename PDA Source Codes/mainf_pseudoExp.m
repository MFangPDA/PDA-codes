% ����ͬ�����������ģ���У�ģ�����ֿ�ȣ��ٺ͹۲�����ֿ��ֵ���бȽ�
% ����ͬ������Ŀ����ԣ���Ϊ��1850-2000��1000-1849����ʱ��Σ�ͬʱ��
% MPI-ESM-PҲ�ӽ���������Ƚϣ�Ȼ����������ͼ
function postf6_pseudoExp()

% �������ļ����е�m�ļ�
addpath(genpath(pwd));

%% ------------------����ȫ�ֱ���-----------------------
% GISTEMP������ֹ��ݼ�����λ��
global startGISTEMP   endGISS   pathGISS;
% MPI���ݵ���ֹ��ݼ�����λ��
global startMPI endMPI pathMPI;
% ��������������ֹ��ݼ�����λ��
global startProxy endProxy pathProxy;
% ģ��ѵ����ֹ��ݼ�ѵ����ʱ�䳤��
global startTrain endTrain lenTrain;
% ����ģ��ģ�����ֹ��ݼ�ģ���ʱ�䳤��
global startSim endSim lenSim;
% ѵ����������ռ�ֱ���
global interval;
% ����ȫ�ֱ���
global startDA  endDA  lenDA;
%
global smoothSpan;
%
% ȫ�ֱ�����ֵ(�̶�ֵ)
startProxy=1000; endProxy=2000;
startGISTEMP=1880;  endGISS=2017;
startTrain=1880; endTrain=2000;
lenTrain=endTrain-startTrain+1;
startSim=startMPI;
endSim=endMPI;
lenSim=endSim-startSim+1;
interval=2;
startDA=1000;   endDA=2000;   lenDA=endDA-startDA+1;
smoothSpan=1;
pathGISS='./Data/air.2x2.250.mon.anom.land_NH.nc';
pathProxy='./Data/NHmultiProxies.mat';

% Step 1: ��ȡGISTEMP���ݲ�����
tas_monMean=subf_readNC(pathGISS,'air');
tas_yrMean=subf_yrMean(tas_monMean);    % ��ƽ����������������
A=size(tas_yrMean);

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=tas_yrMean(:,:,i)';  % ע�⣺GISTEMP���ݹ��򻯷������������ݲ�һ��
end
normGISS=geoFormat;
% ������ռ����Ч��
% plotT(normGISS(:,:,end));

% Step 2: ����ͬ������������ƽֵ
PDA=load('./Result/PDAResult/PDATas.mat');
PDATas=PDA.ANNTas;
% All time series used here are anomalous relative to the means of 1000-2000,
Anomaly=subf_yrAnomaly(PDATas,startDA,endDA,1880,2000);

% Ȼ���ٶԾ�ƽֵ����ƽ������
Smooth=zeros(size(Anomaly));  A=size(Anomaly);
for i=1:A(1,1)
    for j=1:A(1,2)
        Smooth(i,j,:)=smooth(Anomaly(i,j,:),smoothSpan,'lowess');
    end
end
PDA=Smooth;

% Step 3: ��ȡMPI-ESM-P�Ĺ�ȥǧ�꣨1000-2000��ģ�����ݣ��������ƽֵ�����1000-2000��
tasFieldMPI=readPMIP3('MPI');

% Step 4: ��ȡ�о����ڵ����д�������
multiProxies=load(pathProxy);
% ���н�Ҫ��ͬ���Ĵ������ϵĹ۲�ֵ���У�ע���һ�������;��1�����������͡���2-3���Ǿ�γ��
obsProxies=multiProxies.NHProxies;
A=size(obsProxies);
allProxies=A(1,2);        % ��ס����1�������
screenProxy=[];
corMPIvct=[]; ceMPIvct=[];
corPDAvct=[]; cePDAvct=[];

% % Step 5: ����������ģ����ģ���851-1849������֡�R�Լ���γ�ȵľ���
% % ��1�У����ȣ���2�У�γ�ȣ���3�У�R����4�У�����ģ�����������NHmultiProxies.mat�ļ��е��к�
% % ��5��(4+dataLen)�д��ģ��ֵ
% dataLen=endSim-startSim+1;
% LUMData=zeros(4+dataLen,allProxies);    LUMData(1:4,1)=9999;
% ANNData=zeros(4+dataLen,allProxies);    ANNData(1:4,1)=9999;
% for i=5:(4+dataLen)
%     LUMData(i,1)=850+(i-4);
% end
% RMS=zeros(allProxies-1,2);
% Coreff=zeros(allProxies-1,2);

% Step 5: ����ѵ����ѭ������
for i=2:allProxies     % ��ס����1�������
    Index=i;           % index:����������NHmultiProxies.mat�ļ��е��к�
    Proxy=obsProxies(:,Index);
    pxyType=Proxy(1);
    lon=Proxy(2);
    lat=Proxy(3);
    [row,col] = subf_grdNO(lat,lon);
    
    % ��ѡ�߽�֮�ڵ����ֵ㣬�߽�֮���һ�ɲ�Ҫ
    if (-180<=lon && lon<=180)
        if (0<=lat && lat<=90)
            
            % ע�⣺�������ϡ�GISSTEMP��1880-2000�������ݣ����Զ����ڿ�ֵ�ı׶�
            % ��ˣ�ģ��ѵ��ʱ���������ϡ�GISSTEMP����ȡ1880-2000������ݣ����ж�
            % Ȼ�����������ݵķǿ�ֵ�Ľ�������ѵ�����õ�ѵ����ģ��
            
            % obsProxy����1�����������͡���2-3���Ǿ�γ��
            obsProxy=Proxy((startTrain-startProxy+1+3):(endTrain-startProxy+1+3));
            Tas=normGISS(row,col,(startTrain-startGISTEMP+1):(endTrain-startGISTEMP+1));
            obsTas=reshape(Tas,lenTrain,1);
            
            flagGISS=findGISS(obsTas);
            flagProxy=findProxy(obsProxy);
            
            % ��ȡ�������зǿ�ֵ�Ľ���
            intersection=findIntersection(flagGISS,flagProxy);
            
            z=obsProxy(intersection);
            x=obsTas(intersection)-mean(obsTas(intersection));  % ��ס������ȫʱ�ξ�ƽֵ
            
            if (length(z)>=30 && length(x)>=30)   % ֻ�н������ݵĳ��ȴ���30���ع��������
                C=corrcoef(z,x);
                screenProxy=[screenProxy;Index,pxyType,lon,lat,C(1,2)];
                
                %  Step 4: �ڴ˴�ֱ�����ģ��ѵ���ĳ���Ч������
                fprintf(strcat('The',32,num2str(i),'th column proxy is calibrated now!','\n'))
                % ע��32�ǿո��ASCII�룻strcat���Զ�ɾ���������ַ�����ͷβ�Ŀո�
                
                % ѵ��LUMģ��(һԪ���Իع�): Y=b0+b1*X
                [b,bint,r,rint,stats] = regress(z,[ones(size(x)),x]);
                % �õ�ģ��ѵ���Ľ�������  R1=stats(1,4);
                b0=b(1);           % b0
                b1=b(2);           % b1
                
                % �������ͬ�������ģ�������ֱ����ģ��
                % (1)��ȡPDA��Ӧ�����ϵ�����ֵ
                tasPDAgrid=PDA(row,col,:);
                % ��Ҫ������ת����һά����ʽ
                tasPDA=reshape(tasPDAgrid,lenSim,1);
                % (2)ģ��Ĵ����������У�1000-2000��֮��ģ�
                proxyPDA=b0+b1*tasPDA;
                
                % (3)��ȡMPI��Ӧ�����ϵ�����ֵ
                tasMPIgrid=tasFieldMPI(row,col,:);
                % ��Ҫ������ת����һά����ʽ
                tasMPI=reshape(tasMPIgrid,lenSim,1);
                % (4)ģ��Ĵ����������У�1000-2000��֮��ģ�
                proxyMPI=b0+b1*tasMPI;
                
                % ������Ԥ��λ�ã�������Ӹ����ģʽ�������
                
                % (5)��ȡ�����������С�ע�⣺�����������е�ֵ��ʱ�䷶Χ
                % �п�����1000-2000֮���ĳ��ʱ�䣬������Ҫ�ж�
                obsProxy=Proxy(4:end);
                flagProxy=findProxy(obsProxy);
                flagMPI=ones(1,lenDA);
                flagPDA=flagMPI;
                
                % (6)���ҽ��������ҵ���֤��������
                overlap=findIntersection(flagProxy,flagPDA);
                
                verifyProxy=obsProxy(overlap);
                verifyMPI=proxyMPI(overlap);
                verifyPDA=proxyPDA(overlap);
                
                % (7)�������ϵ����CE
                [corMPI,ceMPI]=verifyIDX(verifyMPI,verifyProxy);
                [corPDA,cePDA]=verifyIDX(verifyPDA,verifyProxy);
                
                % (8)�ݴ���ָ֤��
                corMPIvct=[corMPIvct;corMPI];  ceMPIvct=[ceMPIvct;ceMPI];
                corPDAvct=[corPDAvct;corPDA];  cePDAvct=[cePDAvct;cePDA];
            end
        end
    end
end
save('./Result/VerResult/corMPIvct.mat','corMPIvct');
save('./Result/VerResult/corPDAvct.mat','corPDAvct');
save('./Result/VerResult/ceMPIvct.mat','ceMPIvct');
save('./Result/VerResult/cePDAvct.mat','cePDAvct');

% ��ͼ��ֱ��ͼ
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.48,0.48]);
subplot(2,2,1);
plotHistCOR(corMPIvct,'Correlation coefficient','b'); hold on
plotHistCOR(corPDAvct,'Correlation coefficient','r'); hold on
subplot(2,2,2);
plotHistCE(ceMPIvct,'Coefficient of efficiency','b'); hold on
plotHistCE(cePDAvct,'Coefficient of efficiency','r'); hold on
subplot(2,2,3);
ceChange=corPDAvct-corMPIvct;
plotHistCOR(ceChange,'Change in correlation coefficient','g'); hold on
subplot(2,2,4);
ceChange=cePDAvct-ceMPIvct;
plotHistCE(ceChange,'Change in coefficient of efficiency','g'); hold on
end

function plotHistCOR(vct,paraName,colorName)
X=-1:0.05:1;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=-1+(i-1)*0.05;
    tail=-1+i*0.05;
    for j=1:length(vct)
        if vct(j)>=head && vct(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,colorName,'linewidth',1.5);
set(gca,'XLim',[-1 1]); set(gca,'XTick',-1:0.5:1);
set(gca,'YLim',[0 100]); set(gca,'YTick',0:20:100);
xlabel(paraName); ylabel('Frequency');
text(-0.9,80,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',14,'fontweight','bold'); 
% set(gcf,'color','none');   % ͼ�α�����Ϊ��ɫ
set(gca,'color','none');     % �����ᱳ��Ϊ��ɫ����������Ҫ��ͨ��ͼ�α����İ�ɫ��ʵ��Ϊ������Ĭ�ϱ���ɫ��������Ϊ͸�����Ϳ��Խ��ж��ͼ�ĵ���
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end

function plotHistCE(vct,paraName,colorName)
X=-5:0.15:1;
N=length(X); disp(N);
Count=zeros(1,N);
for i=1:N
    head=-5+(i-1)*0.15;
    tail=-5+i*0.15;
    for j=1:length(vct)
        if vct(j)>=head && vct(j)<tail
            Count(i)=Count(i)+1;
        end
    end
end
stairs(X,Count,colorName,'linewidth',1.5);
set(gca,'XLim',[-5 1]); set(gca,'XTick',-5:1:1);
set(gca,'YLim',[0 200]); set(gca,'YTick',0:40:200);
xlabel(paraName); ylabel('Frequency');
text(-4.5,160,strcat('Mean: ',num2str(mean(vct),'%.3f')),'color',colorName,'fontsize',12,'fontweight','bold'); 
% set(gcf,'color','none');   % ͼ�α�����Ϊ��ɫ
set(gca,'color','none');     % �����ᱳ��Ϊ��ɫ����������Ҫ��ͨ��ͼ�α����İ�ɫ��ʵ��Ϊ������Ĭ�ϱ���ɫ��������Ϊ͸�����Ϳ��Խ��ж��ͼ�ĵ���
set(gca,'linewidth',1.5,'fontname','cambria','fontsize',12);
box off;
hold on;
end

% �������ϵ��Cor��CE
function [Cor,CE]=verifyIDX(rec,ref)
Correlation= corrcoef(rec,ref); % corrcoef���������������ÿ��֮������ϵ��
Cor=Correlation(1,2);
CE=1-(sum((ref-rec).*(ref-rec)))/(sum((ref-mean(ref)).*(ref-mean(ref))));
end

function [row,col] = subf_grdNO(lat,lon)
% ���±�д�ĸ��������(lat,lon)����������ʵ����ռ��е�����ŵķ���
% ����NC�����ڱ���ȡʱ������֮ǰ������Ҫ��ת��������ʵ����ռ�һ�µĽṹ

% ���磺�����(89,-179)�ĵ㣬Ӧ������ʵ����ռ�
% �����Ͻǵ���һ������㣬��Ӧ��geoFormation
% �����У���Ӧ�ľ���geoFormation(1,1)������
% (89,-179) --> (1,1)       (89,179) --> (1,180)
% (31,-179)  --> (30,1)    (31,179) -->  (30,180)

interval=2;

% ceil:����ȡ��
row=ceil((90-lat)/interval);
col=ceil((180+lon)/interval);
end

function flagGISS=findGISS(Series)
len=length(Series);
flagGISS=zeros(1,len);
for i=1:len
    if(Series(i)<=15)
        flagGISS(i)=1;  % �ǿ�ֵ�Ļ������Ϊ1
    else
        flagGISS(i)=0;
    end
end
end

function flagProxy=findProxy(Series)
len=length(Series);
flagProxy=zeros(1,len);
for i=1:len
    if(~isnan(Series(i)))
        flagProxy(i)=1;
    else
        flagProxy(i)=0;
    end
end
end

function intersection=findIntersection(flag1,flag2)
len=length(flag1);
intersection=[];
for i=1:len
    if (flag1(i)==1 && flag2(i)==1)  % �������ݶ��ǿ�
        intersection=[intersection,i];
    end
end
end

%% ��ȡBCC,CCSM4,FGOALg2,GISS,IPSL,MPI��ͨ�ó���
function tasField=readPMIP3(modelname)
% ��������ʱ�Σ�085001-200512
% ��ȡԭʼ���ݣ������пռ����ת��
tas_Amon=subf_readNC(strcat('./PMIP3_en/tas_Amon_2deg_',modelname,'_085001_200512_NH.nc'),'tas');
tas_yrMean=subf_yrMean(tas_Amon)-273.5;
A=size(tas_yrMean);

% �Ѷ�������ԭʼ������֯������ʵ����ռ�һ�µĽṹ
% �����������ݵ�1�д�����ʵ����ռ����������1��
geoFormat=zeros(A(1,2),A(1,1),A(1,3));
for i=1:A(1,3)
    geoFormat(:,:,i)=flipud(tas_yrMean(:,:,i)');
end
%
% ֻ��ȡ1000-2000������
field=geoFormat(:,:,(1000-850+1):(2000-850+1));
tasField=subf_yrAnomaly(field,1000,2000,1880,2000);
end
