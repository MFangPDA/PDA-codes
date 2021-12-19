function corField=subf_corField(prediction,reference)
A=size(prediction);
corField=zeros(A(1,1),A(1,2));
for i=1:A(1,1)
    for j=1:A(1,2)
        prediction_series=reshape(prediction(i,j,:),1,A(1,3))';
        reference_series=reshape(reference(i,j,:),1,A(1,3))';
        Cor=corrcoef(prediction_series,reference_series);
        corField(i,j)=Cor(1,2);
    end
end
end
