function plot_corField(data,titleName)

A=size(data);
latm=zeros(A);
lonm=zeros(A);

for i=1:A(1,1)
    latm(i,:)=90-((i-1)*2+1);
end

for j=1:A(1,2)
    lonm(:,j)=-180+((j-1)*2+1);
end

ax = axesm('MapProjection','stereo','maplatlimit',[-0.001 90],'maplonlimit',[-180 180]);    
axis off  
tightmap
setm(ax,'GLineStyle','--','Grid','on','Frame','on','FLineWidth',1,'FFaceColor',[0.8 0.8 0.8],'fontsize',14,'fontweight','bold')   
setm(ax,...
    'MlabelLocation',45,...                           
    'MlineLocation',45,...                             
    'MeridianLabel','off',...                          
    'PlabelLocation',0:30:90,...                      
    'ParallelLabel','off',...                          
    'PlineLocation',0:30:90,...                       
    'MLabelParallel','south' ...                       
    );

geoimg=geoshow(latm,lonm,data,'DisplayType','texturemap');
geoimg.AlphaDataMapping = 'none';   % interpet alpha values as transparency values
geoimg.FaceAlpha = 'texturemap';    % Indicate that the transparency can be different each pixel
alpha(geoimg,double(~isnan(data))); % Change transparency to matrix where if data==NaN --> transparency = 1, else 0.
hold on

load coast
plotm(lat,long,'k-','linewidth',1);
% title(titleName,'position',[-180 90],'fontname','cambria','fontsize',12)

myColor=genColor();
colormap(myColor);
hcb=colorbar('Ticks',-1:0.5:1,'location','SouthOutside','fontname','cambria','fontsize',12);
caxis([-1,1]);
% title(hcb,'Correlation coefficient','fontname','cambria','fontsize',14);
end

function myColor=genColor()
myColor=zeros(64,3);
% matlab colourmap----64*3 matrix
% blue bar
myColor(1:32,1)=linspace(0,1,32);  myColor(1:32,2)=linspace(0,1,32);  myColor(1:32,3)=1;
% red bar
myColor(33:64,1)=1; myColor(33:64,2)=flip(linspace(0,1,32));  myColor(33:64,3)=flip(linspace(0,1,32)); 
end



