function [col, row]=subf_findGrid(lon, lat,interval)
       % 即东经
    col=ceil((lon-(-180))/interval);  % ceil:向上取整
    row=ceil((lat-0)/interval);

end

