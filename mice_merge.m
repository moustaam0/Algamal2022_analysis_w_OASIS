for i = 1:length (DataD);
    temp= DataD{1,i}(1:5);
    temp2{1,i}= temp(:);
end

for i = 1:length (temp2);
    s1 = temp2{1,i};
    for j = 1:length (temp2);
        s2 = temp2{1,j};
        tf3 = strcmp( s1,s2 );
        tf2{1,j}= tf3(:);
    end
    tf{1,i}= tf2(:);
end

tf5= [tf{:}];
tf6= cell2mat(tf5);
[~,B]=max(tf6,[],2);
B2= unique (B);
B2size = size(B2);
for n= 1:(B2size-1);
    one= B2(n,:);
    two= (B2(n+1,:))-1;
    data= DataF(:,one:two);
    fileName= temp2{1,one};
    save (fileName, 'data');
end
one= B2(B2size,:);
two= size (DataF,2);
data= DataF(:,one:two);
fileName= temp2{1,one};
save (fileName, 'data')

movefile P* APP_PV_Mice;