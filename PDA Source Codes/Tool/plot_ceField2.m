function plot_ceField2(data,titleName)

Z=data;
Z=flipud(Z);
 
R = georasterref('RasterSize',size(Z),...
    'Latlim',[0 90],'Lonlim',[-180 180]);              

ax = axesm('MapProjection','stereo'...                 % 'lambert'; 'eqdcylin'
    ,'maplatlimit',[-0.01 90],'maplonlimit',[-180 180]);    

axis off   
setm(ax,'GLineStyle','--','Grid','on','Frame','off','fontsize',10,'fontweight','bold')    
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

colormap(flipud(hot));
hcb=colorbar('location','SouthOutside','fontname','cambria','fontsize',12);
caxis([-5,1]);
% title(hcb,'Root mean square error','fontname','cambria','fontsize',10);

print('-dtiff','-r600',strcat('./',titleName));
hold on
end
