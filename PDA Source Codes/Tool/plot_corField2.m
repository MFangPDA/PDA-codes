function plot_corField(data,titleName)

Z=data;
Z=flipud(Z);
% 注意：利用本程序画图之前，需要对输入的原始数据进行上下反转处理，因为
% 下列函数在画图的时候会从研究区最下边开始填充，因此读出的第1行数据
% 会被放在最下边。而我们的输入的数据，已经和真实地理空间世界一致了

% ****************************************************
R = georasterref('RasterSize',size(Z),...
    'Latlim',[0 90],'Lonlim',[-180 180]);             % 地理栅格数据参考对象(类)

ax = axesm('MapProjection','stereo'...                 % 'lambert':兰伯特投影; 'eqdcylin': 等距离圆柱投影
    ,'maplatlimit',[-0.01 90],'maplonlimit',[-180 180]);    

axis off  %关闭本地坐标轴系统
setm(ax,'GLineStyle','--','Grid','on','Frame','off','fontsize',14,'fontweight','bold')   %  指定网格线形,绘制frame框架
setm(ax,...
    'MlabelLocation',90,...                            % 每隔90度绘制经度刻度标签
    'MlineLocation',90,...                             % 每隔90度绘制经度线
    'MeridianLabel','off',...                          % 暂时不显示经度刻度标签
    'PlabelLocation',0:30:90,...                       % 只在指定值处绘制纬度刻度标签
    'ParallelLabel','off',...                          % 暂时不显示纬度刻度标签
    'PlineLocation',0:30:90,...                        % 在指定值处绘制纬度线
    'MLabelParallel','south' ...                       % 将经度刻度标签放在南方,即下部
    );
geoshow(Z,R,'DisplayType','texturemap','ZData',zeros(size(Z)),'CData',Z);   % 显示地理数据
demcmap(Z);

% 添加海岸线数据-------
load coast
plotm(lat,long,'k-','linewidth',1);
% -----------------------
% title(titleName)

% 添加图标
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
% matlab colourmap中每一组颜色都是64*3的矩阵形式存储
% 蓝色条带
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=1;
% 红色条带
myColor(33:64,1)=1; myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32)); 
end

