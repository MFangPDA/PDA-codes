function plot_corField(data,titleName)

Z=data;
Z=flipud(Z);

R = georasterref('RasterSize',size(Z),...
    'Latlim',[0 90],'Lonlim',[-180 180]);             

ax = axesm('MapProjection','stereo'...                 % 'lambert'; 'eqdcylin'
    ,'maplatlimit',[-0.01 90],'maplonlimit',[-180 180]);    

axis off  
setm(ax,'GLineStyle','--','Grid','on','Frame','off','fontsize',14,'fontweight','bold')   
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
plotm(lat,long,'k-','linewidth',1);
% title(titleName)

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
% matlab colourmap----64*3 matrix
% blue bar
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=1;
% red bar
myColor(33:64,1)=1; myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32)); 
end

