%% ����ά���飨nc���ݣ���ƽ��ֵ���������ƽ����ͨ�����ڼ�����ƽ������
function nc_YearMean=subf_yrMean(nc_Data_Amon)
I=size(nc_Data_Amon);
col_length=I(1);
row_length=I(2);
year_length=I(3)/12;

nc_YearMean=zeros(col_length,row_length,year_length);
for i=1:year_length
    head=(i-1)*12+1; tail=i*12;
    nc_thisyear=mean(nc_Data_Amon(:,:,head:tail),3);
    nc_YearMean(:,:,i)=nc_thisyear;
end
end