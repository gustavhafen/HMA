%   File: 'Tree_Top_Down_phase.m'
%
%   Author(s):  Gustavo A. Puerto-Souza
%   Created on: 2012
%   
%   (C) Copyright 2015 Gustavo Armando Puerto-Souza and Gian-Luca Mariottini, 
%   All Rights Reserved.
% 
% --- begin license - do not edit ---
% 
% This software is provided "as is" under an open source license, with
% no warranty.  The complete license can be found in LICENSE.txt
% 
% --- end license ---

function [Tree, Matches_Data,Final_Clusters] = Tree_Top_Down_phase(Tree, Matches_Data, Initial_Matches, Params, ind_repetitions, Images, Params_Display, option_plot)

%% TOP-DOWN phase
option_plot = 0;
%% 0. Processing Cluster's Data.
% training set
% indices of Inliers of all clusters
  ind_training = (Matches_Data.label == 1); % replaced.. correct and check later
% Inliers of all clusters (Query-keypoint positions)
training_set = Initial_Matches.Query(1:2, ind_training)';
% test set
ind_test = 1:Initial_Matches.num;
ind_test(ind_training) = [];
num_training = sum(ind_training);
num_test = Initial_Matches.num - num_training;
test_set = Initial_Matches.Query(1:2, ind_test)';
% computing labels of the clusters
labels = Matches_Data.ID_clusters(ind_training); % replaced.. correct and check later
if sum(labels) < 0
    Final_Clusters = [];
    display('No matches found!!!')
end
%% 1. Verification  phase:
%  %classes_test = classify(test_set, training_set, labels, 'diagquadratic');
% Look for the closest cluster
for i=1:num_test,
    candidate = test_set(i,:);
    [distance, min_distance] = min(sum((training_set - repmat(candidate, num_training, 1)).^2,2));
    classes_test(i) = labels(min_distance);
end

%% 2. Estimating the new ATs (from the whole cluster)
%Params.Spatial.min_pts_consensus = 4;
Final_Clusters = [];
ind_cluster2remove = [];
ind_solution = Tree.solution_nodes;
ind_classes = Tree.ind_IDs(ind_solution);
num_clusters = length(ind_solution);
for i=1:num_clusters,
    Final_Clusters.Mapping(i) = Tree.nodes(ind_solution(i)).Mapping;
    Final_Clusters.Inliers(i).Training = Initial_Matches.Training(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).Query = Initial_Matches.Query(:, Tree.nodes(ind_solution(i)).ind_inliers);
end

if option_plot == 1,   
    Final_Clusters.plot_type= 1;
    Plot_Multi_Model(Final_Clusters, Images,Params_Display, []); hold on;
    title('Before Recovering')
end
clear Final_Clusters;

for i=1:num_clusters,
    %indices belonging to the class i
    ind_node_i = ind_solution(i);
    class_i = ind_classes(i);
    ind_cluster_i = (classes_test == class_i);
    ind_test_i = ind_test(ind_cluster_i);
    matches = []; % structure for RANSAC
    matches.Training = Initial_Matches.Training(:, ind_test_i);
    matches.Query = Initial_Matches.Query(:, ind_test_i);
    [ind_consensus, aux, residuals_i] = Compute_Consesus_Model(Tree.nodes(ind_node_i).Mapping, matches, Params.RANSAC.error_bound, ind_repetitions(ind_cluster_i));
    %% 3.Recovery phase:
    if ~isempty(ind_consensus),
        num_new_elements = length(ind_consensus);
        ratio = ( Tree.nodes(ind_node_i).ratio*Tree.nodes(ind_node_i).num_matches + num_new_elements )/ (num_new_elements + Tree.nodes(ind_node_i).num_matches);
        Tree.nodes(ind_node_i).ind_matches = [Tree.nodes(ind_node_i).ind_matches, ind_test_i(ind_consensus)];
        Tree.nodes(ind_node_i).ind_inliers = [Tree.nodes(ind_node_i).ind_inliers, ind_test_i(ind_consensus)];
        Tree.nodes(ind_node_i).num_matches = Tree.nodes(ind_node_i).num_matches + num_new_elements;
        Tree.nodes(ind_node_i).residuals(end+1:end+num_new_elements) = residuals_i(ind_consensus);
        Tree.nodes(ind_node_i).ratio = ratio;
        %
        Matches_Data.ID_clusters(ind_test_i(ind_consensus)) = class_i;
        Matches_Data.rep_error(ind_test_i(ind_consensus)) = residuals_i(ind_consensus);
        Matches_Data.label(ind_test_i(ind_consensus)) = 1;    %inliers
    end
end

for i=1:num_clusters,
    Final_Clusters.Mapping(i) = Tree.nodes(ind_solution(i)).Mapping;
    Final_Clusters.Inliers(i).Training = Initial_Matches.Training(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).Query = Initial_Matches.Query(:, Tree.nodes(ind_solution(i)).ind_inliers);
end
Final_Clusters.plot_type= 1;

if option_plot == 1,
    Plot_Multi_Model(Final_Clusters, Images,Params_Display, []); hold on;
    title('After Recovering')
end