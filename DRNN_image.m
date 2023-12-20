close all;
%I = imread('2.png'); 
%I = imread('moon.jpg'); 
I = imread('fruit.jpg'); 
% figure, imshow(I);
% xlim([0 1280])
% ylim([0 800])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.24,0.29])    % 0.2,0.2,0.125,0.145     0.1   0.116
% set(gca,'FontSize',18)   %18  14.4

%报错可能要改投影图像画图精度
%Center_accuracy = 0.01;
[originalRAW,originalCOL,R_double, G_double, B_double, A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l,A_row_gradient,B_row_gradient,C_row_gradient,A_line_gradient,B_line_gradient,C_line_gradient] = Characteristics(I);
[x,y]=meshgrid( (1 : originalRAW) ,fliplr(1:originalCOL));
ori_data =[R_double, G_double, B_double];


%采样
Sample_value = 4000;   %15440
data = [reshape(x, originalRAW*originalCOL,1),reshape(y, originalRAW*originalCOL,1), R_double, G_double, B_double,A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l];
%kmeans
% [kmeans_cluster_idx,Centers]=kmeans(ori_data, 4);
% sample_x = data(:,1);
% sample_y = data(:,2);
% figure,gscatter(sample_x,sample_y,kmeans_cluster_idx);
% xlim([0 1200])        % moon    mountain (1280,800)   moon(1460 991)  friut(1200 1200)
% ylim([0 1200])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.24,0.29])    % 0.2,0.2,0.125,0.145     0.1   0.116
% set(gca,'FontSize',30)   %18  14.4
% xlabel(''); % 隐藏 x 轴的标签
% ylabel(''); % 隐藏 y 轴的标签
% legend('off');


sample_nNeighbors = 90;        % nNeighborsIndex is how many neighbors used to create the knn index, and must be >= nNeighbors + 1
sample_nNeighborsIndex = sample_nNeighbors + 1; % Because the index includes self-edges (each point is its own nearest neighbor)

rep_nNeighbors = 20;
rep_nNeighborsIndex =rep_nNeighbors+1;

[epsilonDi, minptsDi] = deal(0.08,70);    %0.4 2            Person_1 0.21 14 1200 500
n = 50;
Allcenter = [];                            % Center Points
MinDistance = [];                          % The minimum distance from the point in the category to the center point distance
MeanDistance = [];                         % The mean value of the distance from the point in the category to the center point
MaxDistance = [];                          % Angle between the point of the category and the center point*Dis
points_in_class =[];
cluster_idx=cell(1,n);                     % Sample clustering labels
ori_sample_cell = cell(1, n);
d=0.6;
k=5;
sample_centers = [];

for i = 1:n
    %fprintf('第%d次采样\n',i);
    %数据随机抽取样本
    ori_sample_cell{1, i} = datasample(data, Sample_value);
    ori_sample = ori_sample_cell{1, i};
    sample_x =  ori_sample_cell{1, i}(:,1);
    sample_y =  ori_sample_cell{1, i}(:,2);
    
    ori_sample = [ori_sample_cell{1, i}(:,3), ori_sample_cell{1, i}(:,4), ori_sample_cell{1, i}(:,5)];
%     %[cluster_idx{1, i},Centers]=kmeans(sample_data, k);
     cluster_idx{1,i}=dbscan(ori_sample,epsilonDi,minptsDi);
%     
%     rnndbscan = RnnDbscan(ori_sample , sample_nNeighbors, sample_nNeighborsIndex);
%     %rnndbscan = RnnDbscan(sample_data , nNeighbors, nNeighborsIndex,'Method', 'nndescent');
%     rnndbscan.cluster();
%     % Inspect clusters, outliers, and labels
%     cluster=rnndbscan.Clusters;
%     noise=rnndbscan.Outliers;
%     cluster_idx{1,i}=rnndbscan.Labels;
     %figure,gscatter(sample_x,sample_y,cluster_idx{1,i});
    
    label = unique(cluster_idx{1, i});
    len=length(label);              % the number of clusters
    %len is number of categories
    for j = 1:len
        dis1=[];        
        dis2=[];
        labelname=label(j);   
        if labelname == -1
        sample_cluster_point = [ori_sample(find(cluster_idx{1,i}==labelname),:)]; 
        %center = mean(sample_cluster_point,1);
        points_number = length(sample_cluster_point);
        
        %compute the min_dist, mean_dist and cosine
        [center, min_dist,max_dist, mean_dist] = calculateMetrics(sample_cluster_point);
        Allcenter=[Allcenter;center];
        MinDistance = [MinDistance, min_dist];
        MaxDistance=[MaxDistance,max_dist];
        MeanDistance=[MeanDistance,mean_dist];
            continue
        else
        %find all points belong to category  with labele = labelname  
        sample_cluster_point = [ori_sample(find(cluster_idx{1,i}==labelname),:)]; 
        %center = mean(sample_cluster_point,1);
        points_number = length(sample_cluster_point);
        
        %compute the min_dist, mean_dist and cosine
        [center, min_dist,max_dist, mean_dist] = calculateMetrics(sample_cluster_point);
        Allcenter=[Allcenter;center];
        MinDistance = [MinDistance, min_dist];
        MaxDistance=[MaxDistance,max_dist];
        MeanDistance=[MeanDistance,mean_dist];
    end  
    end
end
%% clustering of local representative
% all_representative = [Allcenter,MinDistance',MaxDistance',MeanDistance'];      % local representative
% [coeff, score] = pca(all_representative);     % dimensionality reduction
% res = score(:, 1:2);
% all_representative_DP = res;
% %[distance_idx,p,s,K,C,Klist] = Dpeak(all_representative_DP, d, k);      % clustering of local representative
% %[distance_idx,p,s,K,C,Klist] = Dpeak(all_representative, d, k);      % clustering of local representative
%  rnndbscan = RnnDbscan(all_representative_DP, rep_nNeighbors, rep_nNeighborsIndex);
% % rnndbscan = RnnDbscan(all_representative, nNeighbors, nNeighborsIndex);
%  rnndbscan.cluster
%  % Inspect clusters, outliers, and labels
% cluster=rnndbscan.Clusters;
% noise=rnndbscan.Outliers;
% rep_cluster_idx=rnndbscan.Labels;
% %final_representative_centers
% figure(),gscatter(all_representative_DP(:,1), all_representative_DP(:,2), rep_cluster_idx);
% 
% [rep_num,dim]=size(all_representative);
% %[shapes_id,colors_rgb,label_idx]=drawshapes(all_representative,rep_cluster_idx,rep_num);
% 
% rep_uni_label = unique(rep_cluster_idx);
% cluster_num=length(rep_uni_label);              % the number of clusters
% matrix_center_of_rep_within_same_category=zeros(cluster_num,dim);
% label_id_of_rep_within_same_category=zeros(cluster_num,1);
% for j = 1:cluster_num
%         cur_label_id=rep_uni_label(j);
%         if cur_label_id == -1
%             continue
%         else
%             %save the label id
%             label_id_of_rep_within_same_category(j)=cur_label_id;
%             
%             %find all representatives belong to category  with labele = cur_label_id  
%             rep_of_same_category = all_representative(find(rep_cluster_idx==cur_label_id),:); 
%             
%             %save the center of the current category
%             matrix_center_of_rep_within_same_category(j,:)=mean(rep_of_same_category);
%  
%         end  
% end
% zeroRows = any(label_id_of_rep_within_same_category == 0, 2);
% label_id_of_rep_within_same_category(zeroRows, :) = [];
% matrix_center_of_rep_within_same_category(zeroRows, :) = [];
% zeroIndices = find(zeroRows);
% %figure(),gscatter(matrix_center_of_rep_within_same_category(:,1),matrix_center_of_rep_within_same_category(:,2),label_id_of_rep_within_same_category);
% 
% [overall_label] = calculateClusterLabels(data(:,3:5), matrix_center_of_rep_within_same_category,label_id_of_rep_within_same_category);
% sample_x = data(:,1);
% sample_y = data(:,2);
% figure,gscatter(sample_x,sample_y,overall_label);
% xlim([0 1200])        % moon    mountain (1280,800)   moon(1460 991)  friut(1200 1200)
% ylim([0 1200])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.24,0.29])    % 0.2,0.2,0.125,0.145     0.1   0.116
% set(gca,'FontSize',30)   %18  14.4
% xlabel(''); % 隐藏 x 轴的标签
% ylabel(''); % 隐藏 y 轴的标签
% legend('off');


 %% center
all_representative = [Allcenter,MinDistance',MaxDistance',MeanDistance'];     % local representative

% all_representative_PCA =mapminmax(all_representative',1,10)';       % normalize
% [coeff, score] = pca(all_representative_PCA);     % dimensionality reduction
% res = score(:, 1:2);
% all_representative_PCA = res;
% 
% rnndbscan = RnnDbscan(all_representative_PCA , rep_nNeighbors, rep_nNeighborsIndex);

rnndbscan = RnnDbscan(Allcenter , rep_nNeighbors, rep_nNeighborsIndex);
rnndbscan.cluster
% Inspect clusters, outliers, and labels
cluster=rnndbscan.Clusters;
noise=rnndbscan.Outliers;
rep_cluster_idx=rnndbscan.Labels;
%final_representative_centers
%figure(),gscatter(all_representative_PCA(:,1), all_representative_PCA(:,2),rep_cluster_idx);
    
rep_uni_label = unique(rep_cluster_idx);
cluster_num=length(rep_uni_label);              % the number of clusters

[~,dim]=size(Allcenter);
matrix_center_of_rep_within_same_category=zeros(cluster_num,dim);
label_id_of_rep_within_same_category=zeros(cluster_num,1);
for j = 1:cluster_num
        cur_label_id=rep_uni_label(j);
        if cur_label_id == -1
            continue
        else
            %save the label id
            label_id_of_rep_within_same_category(j)=cur_label_id;
            
            %find all representatives belong to category  with labele = cur_label_id  
            rep_of_same_category = Allcenter(find(rep_uni_label==cur_label_id),:); 
            
            %save the center of the current category
            matrix_center_of_rep_within_same_category(j,:)=mean(rep_of_same_category);
        end  
end
    
[overall_label] = ClusterLabels_center(data(:,3:5), matrix_center_of_rep_within_same_category,label_id_of_rep_within_same_category);

sample_x = data(:,1);
sample_y = data(:,2);
figure,gscatter(sample_x,sample_y,overall_label);
xlim([0 1200])        % moon    mountain (1280,800)   moon(1460 991)  friut(1200 1200)
ylim([0 1200])
set(gcf,'unit','normalized','position',[0.2,0.2,0.24,0.29])    % 0.2,0.2,0.125,0.145     0.1   0.116
set(gca,'FontSize',30)   %18  14.4
xlabel(''); % 隐藏 x 轴的标签
ylabel(''); % 隐藏 y 轴的标签
legend('off');
%figure(),gscatter(dataall(:,1),dataall(:,2),all_label);
  

