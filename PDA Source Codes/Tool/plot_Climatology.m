function plot_Climatology(data,titleName)

Z=data;
Z=flipud(Z);
 
R = georasterref('RasterSize',size(Z),...
    'Latlim',[0 90],'Lonlim',[-180 180]);              

ax = axesm('MapProjection','stereo'...                 % 'lambert'; 'eqdcylin'
    ,'maplatlimit',[-0.01 90],'maplonlimit',[-180 180]);    
axis off  
tightmap
setm(ax,'GLineStyle','--','Grid','on','Frame','on','FLineWidth',1,'FFaceColor',[0.8 0.8 0.8],'fontsize',9)   
setm(ax,...
    'MlabelLocation',90,...                            
    'MlineLocation',90,...                             
    'MeridianLabel','off',...                          
    'PlabelLocation',0:30:90,...                       
    'ParallelLabel','off',...                          
    'PlineLocation',0:30:90,...                        
    'MLabelParallel','south' ...                       
    );
geoshow(Z,R,'DisplayType','texturemap','ZData',zeros(size(Z)),'CData',Z);   
demcmap(Z);

load coast
plotm(lat,long,'k-','linewidth',0.5);
title(titleName,'Fontname','Arial','Fontsize',9)

myColor=genColor();
colormap(myColor);
hcb=colorbar('Ticks',-1:0.5:1,'location','SouthOutside','fontname','Arial','fontsize',9);
caxis([-1,1]);
title(hcb,'Temperature anomaly (\circC)','Fontname','Arial','Fontsize',9)

hold on
end

function myColor=genColor()
myColor=zeros(64,3);
% matlab colourmap    64*3 matrix
flag=0.5451;  % 0-1  DarkBlue:[0 0 0.5451]   DarkRed:[0.5451 0 0]
% blue bar
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=linspace(flag,1,32);
% red bar
myColor(33:64,1)=flip(linspace(flag,1,32)); myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32));
end

