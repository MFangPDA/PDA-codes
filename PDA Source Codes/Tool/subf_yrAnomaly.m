function nc_Anomaly=subf_yrAnomaly(nc_Data,startyear,endyear,headyear,tailyear)
% startyear,endyear£ºnc_Data time range
% headyear,tailyear£ºthe period used to calculate anomaly
A=size(nc_Data);
col_length=A(1);
row_length=A(2);
time_length=A(3);

head=headyear-startyear+1;
tail=tailyear-startyear+1;

nc_Anomaly=zeros(size(nc_Data));
nc_Timemean=zeros(col_length,row_length);
for i=head:tail
    nc_Timemean=nc_Timemean+nc_Data(:,:,i);
end
nc_Timemean=nc_Timemean/(tailyear-headyear+1);

for i=1:time_length
    nc_Anomaly(:,:,i)=nc_Data(:,:,i)-nc_Timemean;
end
end