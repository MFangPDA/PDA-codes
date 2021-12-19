function ceField = subf_ceField(prediction,reference)
A=size(prediction);
ceField=zeros(A(1,1),A(1,2));
for i=1:A(1,1)
    for j=1:A(1,2)        
        rec=prediction(i,j,:); ref=reference(i,j,:);
        ceValue=1-(sum((ref-rec).*(ref-rec)))/(sum((ref-mean(ref)).*(ref-mean(ref))));
        ceField(i,j)=ceValue;
    end
end
end

