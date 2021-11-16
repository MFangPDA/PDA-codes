%% ����ά���飨nc���ݣ��·ֱ��ʣ���������ۼ�����ͨ�����ڼ������ۼƽ�ˮ
function nc_YearSum=subf_yrSum(nc_Data_Amon)
I=size(nc_Data_Amon);
col_length=I(1);
row_length=I(2);
year_length=I(3)/12;

nc_YearSum=zeros(col_length,row_length,year_length);
for i=1:year_length
    head=(i-1)*12+1; tail=i*12;
    nc_thisyear=nc_Data_Amon(:,:,head:tail);
    Temp=0;
    for j=1:12
        % ��ˮ����������NaN,��Ҫ�жϲ���NaN�ĳ�0��������ӻ����
        % missing_value = -9.96921E36f
        Temp=Temp+nc_thisyear(:,:,j);
    end
    nc_YearSum(:,:,i)=Temp;
end
end