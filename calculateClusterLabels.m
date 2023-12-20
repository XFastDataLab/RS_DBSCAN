%二维
% Functionality: 
% 
% Input parameters:
%                 P: the total data set P
%   representatives: all reprentatives 
%         label_idx: the label index of representatives.
% 
% Ouput:
%  overall_label: 
function [overall_label] = calculateClusterLabels(P, representatives, label_idx)
    
    [m,n]=size(P);
    ones_vec=ones(m,1);
    [number_representatives,~]= size(representatives);
    [number_P,~]= size(P);
    overall_dists_for_P_to_representatives=[];
    for i=1: number_representatives
        %get the original cneter coordinate of representatives
        cur_CC=representatives(i,1:2); 
        cur_min_dist=representatives(i,3);
        cur_max_dist=representatives(i,4);
        cur_mean_dist=representatives(i,5);
        
        tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
        ori_dists=pdist2(cur_CC,P)';
        Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
        cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
        overall_dists_for_P_to_representatives = [overall_dists_for_P_to_representatives,cur_dists_for_P_to_representatives];
    end
   [min_values, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
    overall_label = label_idx(col_indices, :);      


%% image
% function [overall_label] = calculateClusterLabels(P, representatives, label_idx)
%     
%     [m,n]=size(P);
%     ones_vec=ones(m,1);
%     [number_representatives,~]= size(representatives);
%     [number_P,~]= size(P);
%     overall_dists_for_P_to_representatives=[];
%     for i=1: number_representatives
%         %get the original cneter coordinate of representatives
%         cur_CC=representatives(i,1:3); 
%         cur_min_dist=representatives(i,4);
%         cur_max_dist=representatives(i,5);
%         cur_mean_dist=representatives(i,6);
%         
%         tmp_shift_from_P_to_CC= P-ones_vec*cur_CC;
%         ori_dists=pdist2(cur_CC,P)';
%         Y=[tmp_shift_from_P_to_CC,ori_dists-cur_min_dist, ori_dists-cur_max_dist,ori_dists-cur_mean_dist];
%         cur_dists_for_P_to_representatives= vecnorm(Y,2,2);
%         overall_dists_for_P_to_representatives = [overall_dists_for_P_to_representatives,cur_dists_for_P_to_representatives];
%     end
%    [min_values, col_indices] = min(overall_dists_for_P_to_representatives, [], 2);
%     overall_label = label_idx(col_indices, :);      
% end



