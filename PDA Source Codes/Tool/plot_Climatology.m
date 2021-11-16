function plot_Climatology(data,titleName)

Z=data;
Z=flipud(Z);
% ע�⣺���ñ�����ͼ֮ǰ����Ҫ�������ԭʼ���ݽ������·�ת������Ϊ
% ���к����ڻ�ͼ��ʱ�����о������±߿�ʼ��䣬��˶����ĵ�1������
% �ᱻ�������±ߡ������ǵ���������ݣ��Ѿ�����ʵ����ռ�����һ����

% ****************************************************
R = georasterref('RasterSize',size(Z),...
    'Latlim',[0 90],'Lonlim',[-180 180]);             % ����դ�����ݲο�����(��)

ax = axesm('MapProjection','stereo'...                 % 'lambert':������ͶӰ; 'eqdcylin': �Ⱦ���Բ��ͶӰ
    ,'maplatlimit',[-0.01 90],'maplonlimit',[-180 180]);    
axis off  %�رձ���������ϵͳ
tightmap
setm(ax,'GLineStyle','--','Grid','on','Frame','on','FLineWidth',1,'FFaceColor',[0.8 0.8 0.8],'fontsize',9)   %  ָ����������,����frame���
setm(ax,...
    'MlabelLocation',90,...                            % ÿ��90�Ȼ��ƾ��ȿ̶ȱ�ǩ
    'MlineLocation',90,...                             % ÿ��90�Ȼ��ƾ�����
    'MeridianLabel','off',...                          % ��ʱ����ʾ���ȿ̶ȱ�ǩ
    'PlabelLocation',0:30:90,...                       % ֻ��ָ��ֵ������γ�ȿ̶ȱ�ǩ
    'ParallelLabel','off',...                          % ��ʱ����ʾγ�ȿ̶ȱ�ǩ
    'PlineLocation',0:30:90,...                        % ��ָ��ֵ������γ����
    'MLabelParallel','south' ...                       % �����ȿ̶ȱ�ǩ�����Ϸ�,���²�
    );
geoshow(Z,R,'DisplayType','texturemap','ZData',zeros(size(Z)),'CData',Z);   % ��ʾ��������
demcmap(Z);

% ��Ӻ���������-------
load coast
plotm(lat,long,'k-','linewidth',0.5);
% -----------------------
title(titleName,'Fontname','Arial','Fontsize',9)

% ���ͼ��
myColor=genColor();
colormap(myColor);
hcb=colorbar('Ticks',-1:0.5:1,'location','SouthOutside','fontname','Arial','fontsize',9);
caxis([-1,1]);
title(hcb,'Temperature anomaly (\circC)','Fontname','Arial','Fontsize',9)

hold on
end

function myColor=genColor()
myColor=zeros(64,3);
% matlab colourmap��ÿһ����ɫ����64*3�ľ�����ʽ�洢
flag=0.5451;  % 0-1֮���һ����   DarkBlue:[0 0 0.5451]   DarkRed:[0.5451 0 0]
% ��ɫ����
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=linspace(flag,1,32);
% ��ɫ����
myColor(33:64,1)=flip(linspace(flag,1,32)); myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32));
end

