function clusters=label2kataegis(X,label)
if isempty(label)
    clusters=[];
    return;
end
catg=unique(label);
[r,~]=size(catg);
clusters=cell(r,1);
for i=1:r
    temp=X(find(label==catg(i)),:);
    clusters{i}=temp;
end