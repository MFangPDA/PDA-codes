function [col, row]=subf_findGrid(lon, lat,interval)
       % ������
    col=ceil((lon-(-180))/interval);  % ceil:����ȡ��
    row=ceil((lat-0)/interval);

end

