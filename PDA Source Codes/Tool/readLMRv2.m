% 读取LMR2018的数据并进行格式标准化
% LMR2018的数据1-2000 AD
function [LMRNHseries,LMRArcNHseries]=readLMRv2()

tasField4D=subf_readNC('.\Data\air_MCruns_ensemble_mean_LMRv2.1.nc','air');
A=size(tasField4D);  % 180,91,20,2001

% LMR的数据先要上下翻转，然后再向左平移180°，以格林威治天文台为对称轴
% 把读出来的原始数据组织成与真实地理空间一致的结构
% 即：重组数据第1行代表真实地理空间中最上面的1行

geoFormat=zeros(A(1,2),A(1,1),A(1,4));
for y=1:2001
    
    tasField3D=tasField4D(:,:,:,y);
    
    tasField=mean(tasField3D,3);
    
    geoFormat(:,:,y)=flipud(tasField(:,:)');
    
end

movegeoFormat=zeros(size(geoFormat));
movegeoFormat(:,1:90,:)=geoFormat(:,91:180,:);
movegeoFormat(:,91:180,:)=geoFormat(:,1:90,:);

figure
subplot(1,2,1);
imagesc(geoFormat(:,:,end));  title('平移之前');  hold on
subplot(1,2,2);
imagesc(movegeoFormat(:,:,end));  title('平移之后');  hold on

% % 96*144*2000
% figure
% imagesc(movegeoFormat(1:48,:,end));

% 然后和提取PMIP3的模拟结果一样，计算数据的距平值
tasFieldAnomaly=subf_yrAnomaly(movegeoFormat(:,:,1000:2000),1000,2000,1961,1990);

% % 提取出1000-2000的数据
% LMRTAS=data2deg(1:45,:,1000:2000);

% NH和Arctic气温序列
NHseries=zeros(1,1001);
arcNHseries=zeros(1,1001);

for i=1:1001
    NHseries(i)=mean2(tasFieldAnomaly(1:45,:,i));
    arcNHseries(i)=mean2(tasFieldAnomaly(1:15,:,i));
end
LMRNHseries=roundn(NHseries,-3);
LMRArcNHseries=roundn(arcNHseries,-3);

t=1000:2000;
figure
plot(t,NHseries,'r','linewidth',1.5); hold on

%
% LMRNHfield=LMRTAS;
%
% LMRTropicNHseries=tropicNHseries;
end


