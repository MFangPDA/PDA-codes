function plot_corField(data,titleName)

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
setm(ax,'GLineStyle','--','Grid','on','Frame','off','fontsize',14,'fontweight','bold')   %  ָ����������,����frame���
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
plotm(lat,long,'k-','linewidth',1);
% -----------------------
% title(titleName)

% ���ͼ��
myColor=genColor();
colormap(myColor);
hcb=colorbar('location','SouthOutside','fontname','cambria','fontsize',12);
caxis([-1,1]);
% title(hcb,'Correlation coefficient','fontname','cambria','fontsize',10);

print('-dtiff','-r600',strcat('./Result/VerResult/',titleName));
hold on
end

function myColor=genColor()
myColor=zeros(64,3);
% matlab colourmap��ÿһ����ɫ����64*3�ľ�����ʽ�洢
% ��ɫ����
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=1;
% ��ɫ����
myColor(33:64,1)=1; myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32)); 
end

