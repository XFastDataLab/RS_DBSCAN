function [Center,MinDis,MeanDis,CosDis] = Localrep(data,rnnidx)
label_idx = rnnidx;
%分类个数
Center = []; 
MinDis = [];
MeanDis = [];
CosDis = [];  
label = unique(label_idx);
len = length(label);
sample_cluster_point = [];
    for j = 1:len
        dis1=[];        
        dis2=[];
        labelname=label(j);   
        if labelname == -1
            continue
        else
        sample_cluster_point = [data(find(label_idx==labelname),:)]; 
        center = mean(sample_cluster_point,1);
        Center=[Center;center];
        [row,col] = size(sample_cluster_point);
        for p=1:row
            D=DS(sample_cluster_point(p,:),center(1,:));
            dis1 = [dis1,D];                                    %Store the distance from the point in the category to the center point
            A=sample_cluster_point(p,:)-center(1,:);
            B=center(1,:);
            cos=dot(A',B)/norm(A)/norm(B);
            thetaPos=acos(dot(A,B)/(norm(A)*norm(B)))*180/pi;    %Angle
            dis2 = [dis2, abs(D*cos)];
        end 
        dis2=dis2(1,1:row-1);
        min_dis = min(dis1');               % minimum distance                         
        mean_dis = mean(dis1',1);           % average distance
        cos_dis =sum(dis2);                 % sum of angle*distance
        end
        MinDis = [MinDis, min_dis];
        MeanDis=[MeanDis,mean_dis];
        CosDis=[CosDis,cos_dis];  
    end
        
    function D = DS(X,Y)
        D = X.^2*ones(size(Y'))+ones(size(X))*(Y').^2-2*X*Y';
    end
end
