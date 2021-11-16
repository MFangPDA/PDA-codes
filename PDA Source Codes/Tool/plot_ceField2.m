function plot_ceField2(data,titleName)

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
setm(ax,'GLineStyle','--','Grid','on','Frame','off','fontsize',10,'fontweight','bold')   %  ָ����������,����frame���
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
colormap(flipud(hot));
hcb=colorbar('location','SouthOutside','fontname','cambria','fontsize',12);
caxis([-5,1]);
% title(hcb,'Root mean square error','fontname','cambria','fontsize',10);

print('-dtiff','-r600',strcat('./Result/VerResult/',titleName));
hold on
end
