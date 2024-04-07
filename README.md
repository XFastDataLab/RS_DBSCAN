# RS_DBSCAN
RS-DBSCAN is based on two key assumptions: (1) The data of each sufficient sample should have similar data distribution, as well as category distribution, to the entire data set; (2) the representative of each category in all sufficient samples conform to Gaussian distribution. It processes data in two stages, one is to classify data in each local sample independently, and the other is to globally classify data by assigning each point to the category of its nearest representative category center.

***********************************************************************************
The RS-DBSCAN program was compiled under Windows using matlab R2016b.
***********************************************************************************

Files
===================================================================================
These program mainly containing:

-startup code named `Huge_speed.m`.

-one main functions of C4Y named `calculateClusterLabels_Large.m`.

-The main code for running the RNN-DBSCAN algorithm named `RnnDbscan.m` and `RNN-DBSCAN tests.ipynb`, main functions of RNN-DBSCAN named `knnIndexToGraphEdges.m`, `knngragh.m` and `knnindex.m`.

-some data sets.

Dataset Format
===================================================================================
The dataset should be saved in a text file in the following format:
-First, prepare the dataset without column numbers and row numbers. When using Matlab to read the given data in the mat or txt format in the project, you can use the load function. By default, the variable name of the data matrix is named points.
For example, the first 3 lines of the sample dataset "reaction_norm.txt" (with 33913 data points and 29 dimensions) are as follows:

1	5.60150038000000	6.98368582000000	10	5.79009655000000	4.36853458000000	5.97079504000000	5.56423456900947	2.71410562000000	4.71469267000000	1.29810532600000	1.34582442400000	6.92467651000000	6.62162131000000	8.60340988000000	5.89977262000000	5.64927517000000	7.66309897000000	5.70834001000000	8.05536334000000	8.74150606000000	5.35851118000000	4.56869890000000	3.57616396000000	4.93070779000000	4.77036496000000	5.68987534000000	4.40919037000000	1.73058051700000
4.12808995000000	1	4.44066121000000	4.13071768000000	7.10646355000000	1.67400523000000	5.45808196000000	10	5.03483437000000	3.92544424000000	3.01669768000000	1.82090620900000	4.27905766000000	3.09534103000000	5.74199875000000	4.64042908000000	2.97484219000000	5.43718540000000	8.01752500000000	7.16370994000000	8.40023083000000	4.30502185000000	8.20240246000000	2.83482325000000	6.47268724000000	3.29867506000000	3.57226876000000	6.52711816000000	6.64724161000000
4.36351015000000	1.59094349200000	2.68672285000000	5.02900111000000	7.08358735000000	1	6.00224257000000	7.35759135095414	5.55801958000000	4.36338514000000	3.46431619000000	1.64193322600000	6.56113006000000	6.00243004000000	7.02559900000000	4.69794277000000	3.70601758000000	6.51251485000000	10	5.67305893000000	8.58003598000000	2.90317924000000	6.99304303000000	5.23052299000000	4.56796891000000	8.42081806000000	2.39853088000000	6.21819910000000	5.23795906000000

```matlab
dataall=load('reaction_norm.txt');
```
An example of quick start
===================================================================================
Step1:
Open the startup code `Huge_speed.m`, and load the data sets we prepared already.
```matlab
dataall=load('reaction_norm.txt');
```

Step2:
Adjust the appropriate parameters such as n (number of samples), sample (number of samples), sample_nNeighbors (RNN parameter for sample clustering), rep_nNeighbors (RNN parameter for global clustering), and perform clustering.
```matlab
%% parameters
sample_nNeighbors = 23;        % nNeighborsIndex is how many neighbors used to create the knn index, and must be >= nNeighbors + 1
sample_nNeighborsIndex = sample_nNeighbors + 1; % Because the index includes self-edges (each point is its own nearest neighbor)
rep_nNeighbors = 25;
rep_nNeighborsIndex =rep_nNeighbors+1;
n = 50;                            % Number of Distributed Sites
sample = 400;                      % Number of samples
array_size =120;
[~,dimension]=size(dataall);
x=1;

%% Sampling
for i = 1:n
    ori_sample_cell{1, i} = datasample(dataall, sample);       %数据随机抽取样本
end

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
rnndbscan = RnnDbscan(all_representative , rep_nNeighbors, rep_nNeighborsIndex);
rnndbscan.cluster
% Inspect clusters, outliers, and labels
cluster=rnndbscan.Clusters;
noise=rnndbscan.Outliers;
rep_cluster_idx=rnndbscan.Labels;
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
    
[overall_label] = calculateClusterLabels_Large(dataall, matrix_center_of_rep_within_same_category,label_id_of_rep_within_same_category);
toc
```

Step3：
Use `calculateClusterLabels_Large.m` to compute the distances from all data points in the dataset to the centroids of representative points, and assign corresponding labels to them.
```matlab
%% reaction
function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
    
    [m,n]=size(P);
    ones_vec=ones(m,1);
    [number_representatives,~]= size(finally_representatives_center);
    overall_dists_for_P_to_representatives=zeros(m,number_representatives);
    for i=1: number_representatives
        %get the original cneter coordinate of representatives
        cur_CC=finally_representatives_center(i,1:29); 
        cur_min_dist=finally_representatives_center(i,30);
        cur_max_dist=finally_representatives_center(i,31);
        cur_mean_dist=finally_representatives_center(i,32);
        
        tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
        ori_dists=pdist2(cur_CC,P)';
        Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
        cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
        overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
    end
   [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
    overall_label = finally_representatives_center_label(col_indices, :);      
end
```
Step4: -Press the “Run” button in Matlab and the code will be run.

Output Format
===================================================================================
Output the runtime of using RS-DBSCAN for this dataset. As shown below:
>> Huge_speed
历时 0.534213 秒。

Reaction n=30,sample=400,sample_nNeighbors = 50,rep_nNeighbors = 35.
