%rnn dbscan 
%clear;
%dataall=load('kdd04_norm.txt');
%dataall=load('aps_norm.txt');
%dataall=load('FMA_norm.txt');
%dataall=load('mocap_norm.txt');
%dataall=load('reaction_norm.txt');
%dataall=load('cloth_new_2008.txt');
%dataall=load('house1.txt');

dataall = load ('Skin.txt');
%dataall = load ('shuttle.txt');
%dataall= load('letter.txt');
%dataall = load ('SensIT_acoustic.txt');
%dataall = load('bean.txt');
%dataall=load('K15.txt');
%dataall=load('accelerate.txt');
%dbscan_cluster_idx = dataall (:,1);
%dataall = points;
dataall=mapminmax(dataall',1,10)';

%bean
% data_idx=dataall(:,16);
% dataall =mapminmax(dataall(:,1:16)',1,10)';
%dataall=[dataall,data_idx];
% [coeff, score]= pca(dataall(:,1:16)); 
% res = score(:, 1:2);
% %res=[res(:,1).*10,res(:,2),res(:,3)];
% dataall=res;
%dbscan_cluster_idx=dbscan(dataall,0.4,80);
% figure(),gscatter(dataall(:,1),dataall(:,2),dbscan_cluster_idx);

%%accelerate
% data_idx=dataall(:,1);
% dataall =mapminmax(dataall',1,10)';
% [coeff, score]= pca(dataall(:,1:5)); 
% res = score(:, 1:3);
% dataall= res;
% dbscan_cluster_idx=dbscan(dataall,1,30);
% figure(),scatter(dataall(:,1),dataall(:,2), dbscan_cluster_idx);

%% parameters
sample_nNeighbors = 23;        % nNeighborsIndex is how many neighbors used to create the knn index, and must be >= nNeighbors + 1
sample_nNeighborsIndex = sample_nNeighbors + 1; % Because the index includes self-edges (each point is its own nearest neighbor)

rep_nNeighbors = 25;
rep_nNeighborsIndex =rep_nNeighbors+1;

n = 50;                            % Number of Distributed Sites
sample = 400;                      % Number of samples
% k = 6;                              % Number of clusters
% d = 10;                              % DPeak parameters d_c
array_size =120;
[~,dimension]=size(dataall);
x=1;

% %% Overall data clustering
%rnndbscan = RnnDbscan(data , nNeighbors, nNeighborsIndex);
%figure(),gscatter(dataall(:,1),dataall(:,2),cluster_idx);
% tic
% dbscan_cluster_idx=dbscan(dataall,0.4,80);
% toc
% label = unique(dbscan_cluster_idx);
% len=length(label);
% X= get_BestMatch(overall_label,dbscan_cluster_idx);
%figure(),gscatter(dataall(:,1),dataall(:,2),cluster_idx);

%% Sampling
for i = 1:n
    ori_sample_cell{1, i} = datasample(dataall, sample);       %数据随机抽取样本
end

% for i = 1:n
%     ori_sample_cell{1, i} = datasample(dataall, sample);       %数据随机抽取样本
%     sample_data = ori_sample_cell{1, i};
%     sample_idx_cell{1,i}=sample_data(:,17);
%     sample_cell{1, i}=sample_data(:,1:16);
%     [coeff, score]= pca(sample_data(:,1:16)); 
%     res = score(:, 1:2);
%     %res=[res(:,1).*10,res(:,2),res(:,3)];
%     ori_sample_cell{1, i}=res;
% end
 %figure(),gscatter(ori_sample_cell{1,1}(:,1),ori_sample_cell{1,1}(:,2),sample_idx_cell{1,1});


%% Clustering of samples
Allcenter = zeros(array_size,dimension);                    % Center Points
MinDistance = zeros(array_size,1);                          % The minimum distance from the point in the category to the center point distance
MeanDistance = zeros(array_size,1);                         % The mean value of the distance from the point in the category to the center point
MaxDistance = zeros(array_size,1);                          % Angle between the point of the category and the center point*Dis
cluster_idx=cell(1,n);                     % Sample clustering labels

tic
%n is the number of sample
for i = 1:n    
    ori_sample = ori_sample_cell{1, i};
    rnndbscan = RnnDbscan(ori_sample , sample_nNeighbors, sample_nNeighborsIndex);
    rnndbscan.cluster
    % Inspect clusters, outliers, and labels
    cluster=rnndbscan.Clusters;
    noise=rnndbscan.Outliers;
    cluster_idx{1,i}=rnndbscan.Labels;
    
    % [cluster_idx{1, i} ,~]= kmeans(ori_sample, 9 );
    % figure(),gscatter(ori_sample(:,1),ori_sample(:,2),cluster_idx{1,i});
    %figure(),drawshapes3d(ori_sample,cluster_idx{1,i},sample);
      
    label = unique(cluster_idx{1, i});
    len=length(label);              % the number of clusters
    %len is number of categories
    
    for j = 1:len
        dis1=[];        
        dis2=[];
        labelname=label(j);
        if labelname == -1
            continue
        else
            %find all points belong to category  with labele = labelname  
            sample_cluster_point = ori_sample(find(cluster_idx{1,i}==labelname),:); 
            %center = mean(sample_cluster_point,1);
            points_number = length(sample_cluster_point);
          
            %compute the min_dist, mean_dist and cosine
            [center, min_dist,max_dist, mean_dist] = calculateMetrics(sample_cluster_point);
            Allcenter(x,:) = center;
            MinDistance(x,:) = min_dist;
            MaxDistance(x,:) = max_dist;
            MeanDistance(x,:) = mean_dist;
            x = x+1;
        end  
    end
end

%% clustering of local representative
all_representative = [Allcenter,MinDistance,MaxDistance,MeanDistance];     % local representative

% all_representative_PCA =mapminmax(all_representative',1,10)';       % normalize
% [coeff, score] = pca(all_representative_PCA);     % dimensionality reduction
% res = score(:, 1:3);
% all_representative_PCA = res;
% 
% rnndbscan = RnnDbscan(all_representative_PCA , rep_nNeighbors, rep_nNeighborsIndex);

rnndbscan = RnnDbscan(all_representative , rep_nNeighbors, rep_nNeighborsIndex);
rnndbscan.cluster
% Inspect clusters, outliers, and labels
cluster=rnndbscan.Clusters;
noise=rnndbscan.Outliers;
rep_cluster_idx=rnndbscan.Labels;
%final_representative_centers
%figure(),gscatter(all_representative_PCA(:,1), all_representative_PCA(:,2),rep_cluster_idx);
%figure(),drawshapes3d(all_representative_PCA,rep_cluster_idx,n_rep);
 
rep_uni_label = unique(rep_cluster_idx);
cluster_num=length(rep_uni_label);              % the number of clusters

[n_rep,dim]=size(all_representative);
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
            rep_of_same_category = all_representative(find(rep_cluster_idx==cur_label_id),:); 
            
            %save the center of the current category
            matrix_center_of_rep_within_same_category(j,:)=mean(rep_of_same_category);
        end  
end

% zeroRows = any(label_id_of_rep_within_same_category == 0, 2);
% label_id_of_rep_within_same_category(zeroRows, :) = [];
% matrix_center_of_rep_within_same_category(zeroRows, :) = [];
% zeroIndices = find(zeroRows);
    
[overall_label] = calculateClusterLabels_Large(dataall, matrix_center_of_rep_within_same_category,label_id_of_rep_within_same_category);
toc

% zeroRows = any(dbscan_cluster_idx == -1, 2);
% dbscan_cluster_idx(zeroRows, :) = [];
% overall_label(zeroRows, :) = [];
% zeroIndices = find(zeroRows);
%figure(),gscatter(dataall(:,1),dataall(:,2),all_label);
  

