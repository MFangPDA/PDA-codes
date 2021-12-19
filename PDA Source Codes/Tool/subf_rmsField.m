function rmsField=subf_rmsField(prediction,reference)
A=size(prediction);
rmsField=zeros(A(1,1),A(1,2));
for i=1:A(1,1)
    for j=1:A(1,2)
        prediction_series=smooth(prediction(i,j,:),11,'sgolay')';
        reference_series=smooth(reference(i,j,:),11,'sgolay')';
        rmsField(i,j)=rms(prediction_series-reference_series);
    end
end
end