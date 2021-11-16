function [Slope,Pvalue]=MKTrendDetection(series,name)
X=series;
N=length(X);
U=zeros(N-1,1);
for t=2:N
    x=X(1:t);
    S=0;
    n=length(x);
    for k=1:(n-1)
        for j=(k+1):n
            S=S+sign(x(j)-x(k));
        end
    end
    VarS=n*(n-1)*(2*n+5)/18;
    if S>0
        Z=(S+1)/sqrt(VarS);
    elseif S==0
        Z=0;
    else
        Z=(S-1)/sqrt(VarS);
    end
    U(t-1)=Z;
end

% significan level
Alpha=1-normcdf(U(end),0,1); 

% variation rate
Qi=zeros(N*(N-1)/2,1);
counter=1;
for k=1:(N-1)
    for j=(k+1):N
        Qi(counter)=(X(j)-X(k))/(j-k);
        counter=counter+1;
    end
end
% Sen's Slope:the trend
Q=median(Qi)*100;  

fprintf('%12s \n',strcat('*******',32,'Mann-Kendall is performing on:',32,name,32,'*******'));
fprintf('%12s \n',strcat('>> AA index trend =',num2str(mean(Q),'%.2f'),'/100 yrs')); 
fprintf('%12s \n',strcat('>> AA index trend is statistically significant at Alpha=',num2str(1-Alpha),32,'level')); 

Slope=Q;
Pvalue=1-Alpha;

% sign(x)£ºSignum function
% x<0,sign(x)=-1
% x=0,sign(x)=0
% x>0,sign(x)=1
end

