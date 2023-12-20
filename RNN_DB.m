data=load('synthesis_2.txt');
% data=load('multiple_grids_noise.mat');
%data=load('moons_with_noise.mat');
nNeighbors = 28;
nNeighborsIndex = 29;


% % Use the NN Descent algorithm to create the knn index; this is much faster than an exhaustive search
 tic
 rnndbscan = RnnDbscan(data, nNeighbors, nNeighborsIndex, 'Method', 'nndescent');
 toc

 tic
 rnndbscan = RnnDbscan(data , nNeighbors, nNeighborsIndex);
 toc

rnndbscan.cluster();
% Or
%cluster(rnndbscan);
cluster_idx=[];    
% Inspect clusters, outliers, and labels
cluster_idx=rnndbscan.Labels;
figure(1),gscatter(data(:,1),data(:,2),cluster_idx);

tic
cluster_idx=dbscan(data,25,50);
toc
figure(2),gscatter(data(:,1),data(:,2),cluster_idx);

