function postf6_toNetCDF

% To load subfunction in tool
addpath(genpath(pwd));

% 1. Load the PDA reconstruction
data=load('./Result/PDAResult/PDATas.mat');
PDATas=data.ANNTas;

% Anomaly with respective to the mean of 1961-1990 
headyr=1000; endyr=2000; times=(endyr-headyr+1);
Anomaly=subf_yrAnomaly(PDATas,headyr,endyr,headyr,endyr);
fieldPDA=Anomaly;

% ����߶ȵ�����д��netCDF�ļ�
% ��ɾ����һ�ֹ����д洢�Ľ��
delete('./Result/tas_anomaly_wrt_1000-2000_Ayr_2deg_PDA_1000_2000_multiProxy.nc');
row=45; col=180;  resolution=2; 
filename='./Result/tas_anomaly_wrt_1000-2000_Ayr_2deg_PDA_1000_2000_multiProxy.nc';
write2NC(filename,fieldPDA,row,col,times,resolution,headyr,endyr);

end

function write2NC(filename,data,row,col,yrLen,resolution,headyr,endyr)

% �ȶ���д��Ĳ���
lon=(-180+(resolution/2)):resolution:(180-(resolution/2));
lat=(90-(resolution/2)):(-resolution):(0+(resolution/2));  % ע�⣺��д��������Ѿ��������Ŀռ��ʽ�������ӱ��������е�
time=headyr:endyr;

myData=data;

% ��ʼ����netCDF�ļ�
nccreate(filename,'tas','Datatype','double',...
    'Dimensions',{'lat',row,'lon',col,'time',yrLen},'Format','classic');

nccreate(filename,'lat','Datatype','double','Dimensions',{'lat',row},'Format','classic');
nccreate(filename,'lon','Datatype','double','Dimensions',{'lon',col},'Format','classic');
nccreate(filename,'time','Datatype','int16','Dimensions',{'time',yrLen},'Format','classic');

% variable attribute
% ncwriteatt(filename,varname,attname,attvalue)
ncwriteatt(filename,'lat','name','longtitude');
ncwriteatt(filename,'lat','units','degree_north');
ncwriteatt(filename,'lat','axis','Y');

ncwriteatt(filename,'lon','name','longtitude');
ncwriteatt(filename,'lon','units','degree_east');
ncwriteatt(filename,'lon','axis','X');

ncwriteatt(filename,'time','name','time');
ncwriteatt(filename,'time','units','year');
ncwriteatt(filename,'time','axis','T');

ncwriteatt(filename,'tas','long_name','Annually Temperature Anomalies');
ncwriteatt(filename,'tas','units','degC');
ncwriteatt(filename,'tas','missing_value',NaN);
ncwriteatt(filename,'tas','statistic','Anomaly');
ncwriteatt(filename,'tas','scale_value',1);
ncwriteatt(filename,'tas','add_offset',0);

% Global attribute   '/'��ȫ�����Ա�־����
ncwriteatt(filename,'/','Description:','Annual mean near_surface_temperature_anomaly');
ncwriteatt(filename,'/','Content:','Millennial PDA outputs with assimilating multi-Proxy');
ncwriteatt(filename,'/','Anomaly period:','w.r.t. the millennial mean');
ncwriteatt(filename,'/','Production date:','May-1-2020');
ncwriteatt(filename,'/','Producer:','Miao Fang');
ncwriteatt(filename,'/','E-Mail:','mfang@lzb.ac.cn');
ncwriteatt(filename,'/','Institute:','Northwest Institute of Eco-Environment and Resources, CAS');
ncwriteatt(filename,'/','Reference:',' ');
ncwriteatt(filename,'/','Funding:','The Strategic Priority Research Program of the Chinese Academy of Sciences(XDA19070103)');

% ���д��
ncwrite(filename,'lon',lon);
ncwrite(filename,'lat',lat);
ncwrite(filename,'time',time);
ncwrite(filename,'tas',myData);

% ��ʾ������Ϣ
S = ncinfo(filename);
datafmt = S.Format;
end

