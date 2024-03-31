% %% household
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:7); 
%         cur_min_dist=finally_representatives_center(i,8);
%         cur_max_dist=finally_representatives_center(i,9);
%         cur_mean_dist=finally_representatives_center(i,10);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end

% % APS 
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:170); 
%         cur_min_dist=finally_representatives_center(i,171);
%         cur_max_dist=finally_representatives_center(i,172);
%         cur_mean_dist=finally_representatives_center(i,173);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


% %% Kdd04
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:16); 
%         cur_min_dist=finally_representatives_center(i,17);
%         cur_max_dist=finally_representatives_center(i,18);
%         cur_mean_dist=finally_representatives_center(i,19);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


% %% FAM
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:519); 
%         cur_min_dist=finally_representatives_center(i,520);
%         cur_max_dist=finally_representatives_center(i,521);
%         cur_mean_dist=finally_representatives_center(i,522);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end

% %% MoCap
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:36); 
%         cur_min_dist=finally_representatives_center(i,37);
%         cur_max_dist=finally_representatives_center(i,38);
%         cur_mean_dist=finally_representatives_center(i,39);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end

% %% reaction
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:29); 
%         cur_min_dist=finally_representatives_center(i,30);
%         cur_max_dist=finally_representatives_center(i,31);
%         cur_mean_dist=finally_representatives_center(i,32);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end

% %% Click
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:12); 
%         cur_min_dist=finally_representatives_center(i,13);
%         cur_max_dist=finally_representatives_center(i,14);
%         cur_mean_dist=finally_representatives_center(i,15);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


%% Skin 
function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
    
    [m,n]=size(P);
    ones_vec=ones(m,1);
    [number_representatives,~]= size(finally_representatives_center);
    overall_dists_for_P_to_representatives=zeros(m,number_representatives);
    for i=1: number_representatives
        %get the original cneter coordinate of representatives
        cur_CC=finally_representatives_center(i,1:4); 
        cur_min_dist=finally_representatives_center(i,5);
        cur_max_dist=finally_representatives_center(i,6);
        cur_mean_dist=finally_representatives_center(i,7);
        
        tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
        ori_dists=pdist2(cur_CC,P)';
        Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
        cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
        overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
    end
   [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
    overall_label = finally_representatives_center_label(col_indices, :);      
end

% % letter 
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%    [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:16); 
%         cur_min_dist=finally_representatives_center(i,17);
%         cur_max_dist=finally_representatives_center(i,18);
%         cur_mean_dist=finally_representatives_center(i,19);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


% %%acoustic
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%    [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:51); 
%         cur_min_dist=finally_representatives_center(i,52);
%         cur_max_dist=finally_representatives_center(i,53);
%         cur_mean_dist=finally_representatives_center(i,54);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


%%Dry Bean
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%    [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:16); 
%         cur_min_dist=finally_representatives_center(i,17);
%         cur_max_dist=finally_representatives_center(i,18);
%         cur_mean_dist=finally_representatives_center(i,19);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end


% %% accelerate
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%    [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:5); 
%         cur_min_dist=finally_representatives_center(i,6);
%         cur_max_dist=finally_representatives_center(i,7);
%         cur_mean_dist=finally_representatives_center(i,8);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end

% %% 15K
% function [overall_label] = calculateClusterLabels(P, finally_representatives_center, finally_representatives_center_label)
%     
%    [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(finally_representatives_center);
%     overall_dists_for_P_to_representatives=zeros(m,number_representatives);
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=finally_representatives_center(i,1:2); 
%         cur_min_dist=finally_representatives_center(i,3);
%         cur_max_dist=finally_representatives_center(i,4);
%         cur_mean_dist=finally_representatives_center(i,5);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives(:,i) = cur_dists_for_P_to_representatives;
%     end
%    [~, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = finally_representatives_center_label(col_indices, :);      
% end