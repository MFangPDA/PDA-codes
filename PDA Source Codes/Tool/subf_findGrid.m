function [col, row]=subf_findGrid(lon, lat,interval)
    col=ceil((lon-(-180))/interval);  
    row=ceil((lat-0)/interval);
end

