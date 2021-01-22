%% --------------------HEA-----------------------
% author:guohu
% email:1378917721@qq.com
% --------------------------------------------------------
% Input arguments:
%   Data:             Mutation data (amino acid position vector), each element represents the mutation position.
%   min_num_clust:    Minimum number of mutations in Kataegis 
% Output arguments:
%   new_Data:       Mutation data (remove mutations that do not belong to any kataegis)
%   new_label:      kataegis label
% --------------------------------------------------------
function [new_Data,new_label,clear_number]=HEA(Data,min_num_clust)
if nargin==1
    min_num_clust=10;
end
X=Data;
[r,c]=size(X);
if r<min_num_clust
    new_Data=X;
    new_label=[];
    clear_number=r;
    return;
end
    
    
p = pdist(X);
D= squareform(p);
mean_distance=sum(D,2)./r;
threshold=sum(mean_distance)./(r);
num_peak=get_num_peak(X);
threshold=threshold/num_peak;    %The smaller, the more kataegis are generated
threshold_step=abs((max(mean_distance)-threshold)/10);%
label=-1*ones(r,1);
centroids=[];
new_centroids=[];
index=9999999999;

while exsit_not_partitioned_points(label) || exsit_changed_centroid(centroids,new_centroids)
    centroids=new_centroids;
    i=min_index(mean_distance,label);
    if i~=-1
        %centroids=[centroids;X(i,:)];
        if can_add_centroids(i,centroids,X,index,threshold)
            centroids=[centroids;X(i,:)];
    
        else
            threshold=threshold+threshold_step;
        end
    end
    label=partition(X,centroids,threshold);
    if i~=-1
        index=evaluate_index(X,label);
    end

    new_centroids=update_centroids(X,label);
end
[X1,label1,clear_number]=clear_small_clusters(X,label,min_num_clust);
new_Data=X1;
new_label=label1;
end
%% subfunction definition------------------------------------------------------
function num_peak=get_num_peak(X)
s=std(X);
n=length(X);
IQR=iqr(X);
h=0.9*min(s,IQR/1.34)*n^(-1/5);
[y,xi,bw] = ksdensity(X,1:0.1:max(X),'Bandwidth',h);
extrMaxValue = y(find(diff(sign(diff(y)))==-2)+1);
N1=length(extrMaxValue);
len=length(y);
N2=0;
for i=2:len-2
    if y(i)>y(i-1)
        j=i+1;
        while j<len
            if y(j)==y(i)
                j=j+1;
            else
                break;
            end
        end
        if j>len
            continue;
        end
         if j>i+1
             if y(j)<y(i)
                 N2=N2+1;
             end
         end
    end
end
num_peak=N1+N2;
if y(1)>=y(2)
    num_peak=num_peak+1;
end
if y(len)>=y(len-1)
    num_peak=num_peak+1;
end
if num_peak==0
    num_peak=1;
end
% figure
% plot(xi,y);
% disp(num_peak);
end
function [new_label,Number_mark]=mark_small_clusters(X,label)
new_X=X;
new_label=label;
Number_mark=0;
for i=1:max(label)
   N=length(find(label==i));
   if N<15
       new_label(find(new_label==i))=-1;
       Number_mark=Number_mark+1;
   end      
end
end
function [new_X,new_label,Num_clear]=clear_small_clusters(X,label,min_num_clust)
new_X=X;
new_label=label;
Num_clear=0;
for i=1:max(label)
   N=length(find(label==i));
   if N<min_num_clust
       new_X(find(new_label==i),:)=[];
       new_label(find(new_label==i))=[];
       Num_clear=Num_clear+1;
   end      
end
end
function MSE=MSE_index(X,label)
D=X(find(label~=-1),:);
label=label(find(label~=-1));
if max(label)==1
    MSE=9999999999;
    return
end
s=[];
for i=1:max(label)
    cluster=D(find(label==i),:);
    center=mean(cluster);
    temp=cluster-center;
    temp2=mean(vecnorm(temp,2,2).^2);
    s=[s;temp2];
end
MSE=mean(s);
end
function DBI=DBI_index(X,label)
D=X(find(label~=-1),:);
label=label(find(label~=-1));
if max(label)==1
    DBI=9999999999;
    return
end
eva_DBI = evalclusters(D,label,'DaviesBouldin');
DBI=eva_DBI.CriterionValues;
end
function mark=can_add_centroids2(i,centroids,X,DBI,threshold)
mark=0; 
 centroids=[centroids;X(i,:)];
 label=partition(X,centroids,threshold);
 D=X(find(label~=-1));
 label=label(find(label~=-1));
 if max(label)==1
     mark=1;
     return
 end
 new_DBI=MSE_index(D,label);
 if new_DBI<DBI
     mark=1;
 end
end
function mark=can_add_centroids(i,centroids,X,DBI,threshold)
mark=0; 
 centroids=[centroids;X(i,:)];
 label=partition(X,centroids,threshold);
 D=X(find(label~=-1));
 label=label(find(label~=-1));
 if max(label)==1
     mark=1;
     return
 end
 new_DBI=evaluate_index(D,label);
 if new_DBI<=DBI
     mark=1;
 end
end
function new_index=evaluate_index(D,label)%***************************evaluate index
%better index is less
    D=D(find(label~=-1));
    label=label(find(label~=-1));
    if max(label)==1
        new_index=9999999999;
        return
    end
        
    %     new_MSE=MSE_index(D,label);
    
    eva_DBI = evalclusters(D,label,'DaviesBouldin');
    new_DBI=eva_DBI.CriterionValues;
    
    new_index=new_DBI;
end
% function ch=check_new_centroid(X,centroids,old_DBI)
%     if isempty(centroids)
%         ch=1
%         return
%     end
%     [r,~]=size(centroids);
%     label=kmeans(X,r);
%     eva1_DBI = evalclusters(X,label,'DaviesBouldin')
%     new_DBI=eva1_DBI.CriterionValues;
%     if new_DBI>=old_DBI
%         ch=0
%     else
%         ch=1
%     end
% end
function partition_mark=exsit_not_partitioned_points(label)
N=length(label);
for i=1:N
    if label(i)==-1
        partition_mark=1;
        return 
    end
end
partition_mark=0;
return
end
function variation_centroids_mark=exsit_changed_centroid(centroids,new_centroids)
if length(centroids)~=length(new_centroids)
    variation_centroids_mark=1;
    return
end
[r,c]=find(centroids~=new_centroids);
if isempty(r)
    variation_centroids_mark=0;
    return 
else
    variation_centroids_mark=1;
    return
end
end
function ii=min_index(mean_distance,label)
N=length(label);
m=9999999999999;
ii=-1;
for i=1:N
    if label(i)==-1
        if m>mean_distance(i)
            m=mean_distance(i);
            ii=i;
        end
    end
end
end

function label=partition(X,centroids,threshold)
[r,c]=size(X);
label=-1*ones(r,1);
for i=1:r
    [d,index]=min(vecnorm(X(i,:)-centroids,2,2));
    if d<=threshold
        label(i)=index;  
    else
        label(i)=-1;
    end
    
end
% reduce some centroids  no points belong to 
p_label=label(label~=-1);
only_label=unique(p_label);
only_label=sort(only_label);
for i=1:length(only_label)
    p_label(p_label==only_label(i))=i;
end
label(label~=-1)=p_label;
end
function new_centroids=update_centroids(X,label)
clusters_num=max(label);
new_centroids=[];
for i=1:clusters_num
    c=X(find(label==i),:);
    if size(c,1)==1
        new_centroids=[new_centroids;c];
    else
        new_centroids=[new_centroids;mean(c)];
    end
    
end
end


    

