function [out]=slidingBaseline2(ndata,windowSize,quantileThresh)
%load(sprintf('%s/corrdata/%s.mat', Project_path, name));
%pre-allocate the vector and then map the data into the pieces between the
%window fragments, then map the first window/2 values into the first pad
%and the last window/2 values into the last pad.

data2 = ndata;
out=zeros(size(data2));

if mod(windowSize,2)
    windowSize=windowSize+1;
else
end

halfInd=fix(windowSize/2);



padD=zeros(size(data2,1),size(data2,2)+windowSize);
%disp(size(padD))

padD(:,halfInd+1:end-halfInd)=data2;
padD(:,1:halfInd)=data2(:,halfInd+1:(halfInd+1)+(halfInd-1));
padD(:,end-(halfInd-1):end)=data2(:,end-(halfInd-1):end);
%disp(size(padD))

startQuant=halfInd;
if size(data2,1)>1
    for n=1:size(data2,2)
        td=padD(:,(n+startQuant)-startQuant:(n+startQuant)+(startQuant-1));
        aa=quantile(td,quantileThresh,2);
        out(:,n)=diag(aa);
    end
else
    for n=1:size(data2,2)
        td=padD(:,(n+startQuant)-startQuant:(n+startQuant)+(startQuant-1));
        aa=quantile(td,quantileThresh);
        out(:,n)=aa;
    end
end

out = max(out,1);

out=out;

