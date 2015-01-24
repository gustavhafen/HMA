%   File: 'HMA_Algorithm.m'
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

function [Final_Clusters] = HMA_Algorithm(Features, Parameters, Images, Initial_Matches, ratio_threshold, model)

% Hierarchical Multi-Affine algorithm (TMI 2012).
%%
bool_plot = 0;
%% 0. Loading data/Initializing variables
% a) Loading Parameters
%% this section is only for testing version.
if (nargin < 2) || (isempty(Parameters)),
    if (nargin == 7),
        Parameters = Set_HMA_parameters(Features.Type, model);
    else
        Parameters = Set_HMA_parameters(Features.Type);
    end
end
% b) Loading/Computing Initial Matches
if (nargin < 4)||(isempty(Initial_Matches)),
    Initial_Matches = Set_Initial_Matches(Features.Training, Features.Query, Parameters.Matching);
end
% c) Visual Options
if (nargin >= 3)&&(~isempty(Images)),
    bool_plot = 0;
    Params_Display = Set_Display_Params();
    %Params_Display.Matches.LineStyle = '-';
    
    % ADD options to plot.
end
% d) optional ratio threshold
if nargin == 5,
    Parameters.Clustering.ratio_threshold = ratio_threshold;
end
%%

%% 0) Initialization

factor_children = 1;

% MATCHES STRUCTURE:

% Matches_Data: Struct with the information of the matches in the tree.
% This structure contains the following fields:
% ind_clusters_Tree: match's current index of the cluster wrt ot the tree structure
% rep error: rep error of the match, Inf if they do not have AT nor cluster
% label: label indicating the class of the match: 0-outlier, 1-inlier, 2-limbo
% tree_level: depth level of the cluster in the tree (last-clustering)

Matches_Data = struct(...
    'ID_clusters', ones(Initial_Matches.num, 1), ...
    'rep_error',Inf*ones(Initial_Matches.num, 1), ...
    'label', zeros(Initial_Matches.num, 1) ...
    );
%'tree_level',  zeros(Initial_Matches.num, 1));

% TREE STRUCTURE

% Structure: inner_nodes (nodes to expand).
% This is the control structure for the expansion of inner nodes, and contains the following fields:
index_ID = 1;
tree_depth = 1;
% leaf: boolean var that specifies if the node is leaf (1) or inner (0)
% ind_matches: matches of the cluster
% solution: if the node is part of the final solution set.
root = struct(...
    'ID', index_ID,... %node ID
    'parent', [],... %parent node
    'children',[], ... %children nodes;
    'leaf_flag', 0, ... % leaf node flag
    'solution_flag', 0,... %solution flag
    'num_matches', Initial_Matches.num,... % number of elements;
    'ind_matches', 1:Initial_Matches.num,... % indices of the matches in the cluster
    'Mapping', [],... % Associated AT
    'residuals', Inf*ones(Initial_Matches.num, 1),... % reprojection errors of the matches in the cluster.
    'ind_inliers',[],... %indices of inliers
    'ind_outliers',[],... %indices of outliers
    'num_inliers', 0, ... %number of inliers
    'ratio',0, ...% ratio of inliers per cluster
    'depth', 1 ...
    );


Tree = struct(...
    'nodes', root, ... % vector of nodes
    'num_nodes', 3, ... % number of nodes in the tree
    'depth', tree_depth+1,... % maximal depth of the tree
    'ind_IDs', [1 2 3],... % indices of the ID's of each node in the tree.
    'ind_inner_nodes', 2, ... % indices of not expanded inner nodes in the Tree structure.
    'ind_leaf_nodes', 3, ... % indices of the leaf nodes in the Tree structure.
    'solution_nodes', [], ...
    'num_inner_nodes', 1 ...
    ); % Solution- node's indices in the Tree structure.

% Tree Construction
% deltas = [size_of_images(1)/4, size_of_images(2)/4, size_of_images(1)/2, size_of_images(2)/2, 2*2, pi/2];
% deltas = [size_of_images(1)/4, size_of_images(2)/4, size_of_images(1)/2*(2/3), size_of_images(2)/2*(2/3), 2*8, 2*pi];

%computed indices of repeated matches
ind_repetitions = compute_repeated_positions(Initial_Matches);
% HOUGH filtering
size_of_images = size(Images.Query);
deltas = [size_of_images(1)/4, size_of_images(2)/4, size_of_images(1)/4, size_of_images(2)/4, 2*8, 2*pi];
Hough_Clustering = HMA_Hough_Transform(Initial_Matches, size_of_images, 3, deltas, Images);
% re-initializing Tree
tree_depth = tree_depth + 1;
index_ID = index_ID + 1;
num_matches = length(Hough_Clustering.Clusters);
node1 = struct(...
    'ID', index_ID,... %node ID
    'parent', root.ID,... %parent node
    'children',[], ... %children nodes;
    'leaf_flag', 0, ... % leaf node flag
    'solution_flag', 0,... %solution flag
    'num_matches', num_matches,... % number of elements;
    'ind_matches', Hough_Clustering.Clusters,... % indices of the matches in the cluster
    'Mapping', [],... % Associated AT
    'residuals', Inf*ones(num_matches, 1),... % reprojection errors of the matches in the cluster.
    'ind_inliers',[],... %indices of inliers
    'ind_outliers',[],... %indices of outliers
    'num_inliers', 0, ... %number of inliers
    'ratio',0, ...% ratio of inliers per cluster
    'depth', tree_depth...
    );

index_ID = index_ID + 1;
num_matches = length(Hough_Clustering.removed);
node2 = struct(...
    'ID', index_ID,... %node ID
    'parent', root.ID,... %parent node
    'children',[], ... %children nodes;
    'leaf_flag', 1, ... % leaf node flag
    'solution_flag', 0,... %solution flag
    'num_matches', num_matches,... % number of elements;
    'ind_matches', Hough_Clustering.removed,... % indices of the matches in the cluster
    'Mapping', [],... % Associated AT
    'residuals', Inf*ones(num_matches, 1),... % reprojection errors of the matches in the cluster.
    'ind_inliers',[],... %indices of inliers
    'ind_outliers',[],... %indices of outliers
    'num_inliers', 0, ... %number of inliers
    'ratio',0, ...% ratio of inliers per cluster
    'depth', tree_depth ...
    );

Tree.nodes(2) = node1;
Tree.nodes(3) = node2;
clear node1 node2;

while Tree.num_inner_nodes, % while there are inner nodes
    %% 1) Clustering
    ind_inner_nodes_parents = Tree.ind_inner_nodes;
    inner_nodes = Tree.nodes(ind_inner_nodes_parents); % copying inner nodes
    Tree.ind_inner_nodes = [];
    num_inner_nodes = Tree.num_inner_nodes;
    Tree.num_inner_nodes = 0;
    for j = 1:num_inner_nodes,
        expanded_nodes = [];
        tree_depth = tree_depth + 1;
        leaves = [];
        IDs = [];
        children = [];
        num_expanded_nodes = 0; % counter that measures how many nodes are expanded per level in the tree
        % computing weighting function
        Weights = Weighting_Function(inner_nodes(j).residuals, Parameters.RANSAC.error_bound);
        % clustering with k-means
        [clusters, num_clusters] = Tree_Clustering(Initial_Matches, Weights, 0, Parameters.RANSAC.min_pts_consensus, inner_nodes(j).ind_matches);
        if bool_plot, % plotting
            Plot_Clustering(Initial_Matches.Training(1:2,inner_nodes(j).ind_matches), Initial_Matches.Query(1:2,inner_nodes(j).ind_matches), clusters, Images);
        end
        %% Affine Estimation
        k = 0;
        % index for the num of cluster per inner node
        while  k < num_clusters, % for each cluster
            k = k + 1;
            index_ID = index_ID + 1;
            num_expanded_nodes = num_expanded_nodes + 1;
            % a) new node construction
            cluster_k_elements = inner_nodes(j).ind_matches(clusters == k);
            num_cluster_k_elem = length(cluster_k_elements);
            expanded_nodes.nodes(num_expanded_nodes) = struct(...
                'ID', index_ID,... %node ID
                'parent', inner_nodes(j).ID,... %parent node
                'children',[], ... %children nodes;
                'leaf_flag', 0, ... % leaf node flag
                'solution_flag', 0,... %solution flag
                'num_matches', num_cluster_k_elem,... % number of elements;
                'ind_matches', cluster_k_elements,... % indices of the matches in the cluster
                'Mapping', [],... % Associated AT
                'residuals', [],... % reprojection errors of the matches in the cluster.
                'ind_inliers',[],... %indices of inliers
                'ind_outliers',[],... %indices of outliers
                'num_inliers', 0, ... %number of inliers
                'ratio', 0, ...% ratio of inliers per cluster
                'depth', tree_depth ...
                );
            leaves(num_expanded_nodes) = 0;
            IDs(num_expanded_nodes) = index_ID;
            if bool_plot, % plotting
                Plot_Initial_Matches_subset(Initial_Matches, cluster_k_elements, Images);
            end            
            % b) fitting a new model
            %    Case 1: The cluster does not have enough elements to fit a model
            %            then the elements are labeled as outliers.
            if num_cluster_k_elem < Parameters.RANSAC.min_pts_consensus,
                % The node is a leaf
                leaves(num_expanded_nodes) = 1;
                % Updating the Matches_Data struct
                Matches_Data.ID_clusters(cluster_k_elements) = index_ID;
                Matches_Data.rep_error(cluster_k_elements) = Inf; %we couldnt fit an AT
                Matches_Data.label(cluster_k_elements) = 0; % outliers
                expanded_nodes.nodes(num_expanded_nodes).leaf_flag = 1;
            else % enough points to fit a model
                % preprocessing (packing structure)
                matches_cluster_k.num = num_cluster_k_elem;
                matches_cluster_k.Query = Initial_Matches.Query(:,cluster_k_elements);
                matches_cluster_k.Training = Initial_Matches.Training(:,cluster_k_elements);
                % computing the AT
                [AT_cluster_k, ind_cluster_k_inliers, residuals_cluster_k] = RANSAC_Mapping_Model(matches_cluster_k, Parameters.RANSAC, Parameters.Model, ind_repetitions(cluster_k_elements));
                % CASE 2: The cluster does not satisfy the minimal consensus.
                % The elements are considered as outliers
                if isempty(AT_cluster_k), %it is not possible to fit an AT
                    % The node is a leaf
                    leaves(num_expanded_nodes) = 1;
                    % Updating the Matches_Data struct
                    Matches_Data.ID_clusters(cluster_k_elements) = index_ID;
                    Matches_Data.rep_error(cluster_k_elements) = Inf;
                    Matches_Data.label(cluster_k_elements) = 0; %outliers
                    expanded_nodes.nodes(num_expanded_nodes).leaf_flag = 1;
                    % CASE 3: The Alg. could estimate an AT for the cluster.
                else
                    if length(ind_cluster_k_inliers) < num_cluster_k_elem-4, %4 elements to count the quartiles
                    % Detecting large outliers
                    cluster_k_inliers = cluster_k_elements(ind_cluster_k_inliers); %indices wtr to the indices of initial matches. (absolute indices)
                    [cluster_k_outliers, ind_k_outliers] = setdiff(cluster_k_elements, cluster_k_inliers); %% earsier if a simple if is used (wrt ot the residuals)
                    p_100 = 100*(0.25:0.25:1);
                    residuals_outliers = residuals_cluster_k(ind_k_outliers);
                    quartiles = prctile(residuals_outliers, p_100);
                    max_outliers_thres = quartiles(3) + (quartiles(3)-quartiles(1)) * Parameters.Clustering.factor_IQR; %extended IQR (like in boxplots)
                    % determine the threshold to consider a data-point as outlier
%                     max_outliers_thres = min( max(max_outliers_thres,Parameters.Clustering.min_ERROR), Parameters.Clustering.max_ERROR);
                    % computing outliers (special case)
                    ind_outliers = residuals_outliers > max_outliers_thres; %checking if the residuals are larger than the threshold
                    outliers_to_remove = ind_k_outliers(ind_outliers); %(indices wrt to the residuals)
                    %storing matches to remove
%                     cluster_k_p_elements = cluster_k_outliers(ind_outliers);
                    num_cluster_k_p_elem = sum(ind_outliers);
%                     residuals_cluster_k_p = residuals_outliers(ind_outliers);
                    
                    % removing data from the outliers
                    residuals_inliers = residuals_cluster_k(ind_cluster_k_inliers);
%                     residuals_cluster_k(outliers_to_remove) = [];
                    num_cluster_k_elem = num_cluster_k_elem - num_cluster_k_p_elem;
                    % %                   % auxiliar for plotting
                    cluster_k_inliers = cluster_k_elements(ind_cluster_k_inliers);
                    cluster_k_outliers(ind_outliers) = [];
                    cluster_k_elements = [cluster_k_inliers, cluster_k_outliers];
                    residuals_outliers(ind_outliers) = [];
                    residuals_cluster_k = [residuals_inliers, residuals_outliers];
                    num_inliers = length(cluster_k_inliers);
                    else
                        cluster_k_inliers = cluster_k_elements(ind_cluster_k_inliers);
                        num_inliers = length(cluster_k_inliers);
                        [cluster_k_outliers, ind_k_outliers] = setdiff(cluster_k_elements, cluster_k_inliers); %% earsier if a simple if is used (wrt ot the residuals)
                        ind_outliers = 0;
                    end

                    %                     Inliers.Training = Initial_Matches.Training(:, cluster_k_inliers);
                    %                     Inliers.Query = Initial_Matches.Query(:, cluster_k_inliers);
                    % Storing the solution:
                    expanded_nodes.nodes(num_expanded_nodes).num_matches = num_cluster_k_elem;
                    expanded_nodes.nodes(num_expanded_nodes).ind_matches = cluster_k_elements;
                    expanded_nodes.nodes(num_expanded_nodes).Mapping = AT_cluster_k;
                    expanded_nodes.nodes(num_expanded_nodes).residuals = residuals_cluster_k;
                    expanded_nodes.nodes(num_expanded_nodes).ind_inliers = cluster_k_inliers;
                    expanded_nodes.nodes(num_expanded_nodes).ind_outliers = cluster_k_outliers;
                    expanded_nodes.nodes(num_expanded_nodes).num_inliers = num_inliers;
                    expanded_nodes.nodes(num_expanded_nodes).ratio = num_inliers/num_cluster_k_elem;
                    %
                    Matches_Data.ID_clusters(cluster_k_elements) = index_ID;
                    Matches_Data.rep_error(cluster_k_elements) = residuals_cluster_k;
                    Matches_Data.label(cluster_k_outliers) = 0;
                    Matches_Data.label(cluster_k_inliers)  = 1;
                    if bool_plot, % plotting
                        Inliers.Training = Initial_Matches.Training(:, cluster_k_inliers);
                        Inliers.Query = Initial_Matches.Query(:, cluster_k_inliers);
                        Sol.Mapping = AT_cluster_k;
                        Sol.Inliers = Inliers;
                        Sol.plot_type = 1;
                        %plotting inliers
                        Plot_Multi_Model(Sol, Images, Params_Display, []); hold on;
                        %plotting outliers
                        if ~isempty(cluster_k_outliers),
                            Plot_Initial_Matches_subset(Initial_Matches, cluster_k_outliers, [], 'r', 'x', size(Images.Training,2));
                        end
                        ind_aux = logical((residuals_cluster_k > 3*Parameters.RANSAC.error_bound).*(residuals_cluster_k < Parameters.Clustering.limbo_factor_threshold*3*Parameters.RANSAC.error_bound));
                        cluster_k_limbo = cluster_k_elements(ind_aux);
                        if ~isempty(cluster_k_limbo),
                            Plot_Initial_Matches_subset(Initial_Matches, cluster_k_limbo, [], 'g', 'o', size(Images.Training,2));
                        end
                        if sum(ind_outliers)>0,
                            Plot_Initial_Matches_subset(matches_cluster_k, outliers_to_remove , [], 'y', 's', size(Images.Training,2),1);
                        end
                        clear Sol;
                    end
                    % CASE 3-bis: If there are large outliers we isolate
                    % them in new nodes which are defined as outliers
                    if sum(ind_outliers) > 0,
                        num_clusters = num_clusters + 1; % augmenting number of clusters
                        % getting the real indices of outliers
                        ind_clust = find(clusters == k);
                        ind_clust = ind_clust(outliers_to_remove);
                        clusters(ind_clust) = num_clusters;
                    end
                end % end for AT estimation (Cases 2 and 3)
                %adding the children indices of the node
            end % end for AT estimation (Cases 1, 2 and 3)
        end % end for the cluster j (variable k)
        
        %% HERE isolate leaf nodes
        leaves = logical(leaves);
        num_leaves = sum(leaves);
        %% NEW (April 17)
        if num_leaves == length(leaves), % all expanded clusters are bad
            ind_parent = find(Tree.ind_IDs == inner_nodes(j).ID);
            Tree.nodes(ind_parent).leaf_flag = 1;
            Tree.ind_leaf_nodes = [Tree.ind_leaf_nodes, ind_parent];
            if ~isempty(Tree.nodes(ind_parent).ind_inliers), 
                Tree.nodes(ind_parent).solution_flag = 1;
                Tree.solution_nodes = [Tree.solution_nodes, ind_parent];
            %
                cluster_k_elements = Tree.nodes(ind_parent).ind_matches;
                cluster_k_outliers = Tree.nodes(ind_parent).ind_outliers;
                cluster_k_inliers = Tree.nodes(ind_parent).ind_inliers;
            %%
                    Matches_Data.ID_clusters(cluster_k_elements) = Tree.nodes(ind_parent).ID;
                    Matches_Data.rep_error(cluster_k_elements) =  Tree.nodes(ind_parent).residuals;
                    Matches_Data.label(cluster_k_outliers) = 0;
                    Matches_Data.label(cluster_k_inliers)  = 1;
            end
        else
            % copying the leaf nodes
            leaf_nodes.nodes = expanded_nodes.nodes(leaves);
            % removing them from the current nodes
            expanded_nodes.nodes(leaves) = [];
            %% 3) Recovery phase:
            % function that reassigns the matches according to spatial and
            % geometrical constraints:
            if ~isempty(expanded_nodes),
                [expanded_nodes, Matches_Data] = HMA_Local_Recovery_Phase(expanded_nodes, Initial_Matches, Matches_Data, Parameters, ind_repetitions, Images, Params_Display, bool_plot);
                
                %% 4. Stop Criteria
                num_clusters = length(expanded_nodes.nodes);
                leaf_flag = zeros(1, num_clusters); %vector with the indices of the clusters to remove
                solutions = zeros(1, num_clusters); %vector with the indices of the clusters to remove
                IDs_exp = zeros(1, num_clusters); %vector with the indices of the clusters to remove
                for i = 1:num_clusters, %for each cluster, check if it is possible to expand it
                    IDs_exp(i) = expanded_nodes.nodes(i).ID;
                    num_inliers_i = expanded_nodes.nodes(i).num_inliers; %number of inliers in cluster i
                    num_outliers_i = expanded_nodes.nodes(i).num_matches - num_inliers_i; %number of outliers in cluster i
                    %if (expanded_nodes.nodes(i).ratio >= Parameters.Clustering.ratio_threshold)||...
                     if (num_inliers_i>0)&&(num_outliers_i<Parameters.RANSAC.min_pts_consensus)&&(expanded_nodes.nodes(i).num_inliers < Parameters.Clustering.max_num_inliers),
                        leaf_flag(i) = 1;
                        solutions(i) = 1;
                        expanded_nodes.nodes(i).leaf_flag = 1;
                        expanded_nodes.nodes(i).solution_flag = 1;
                    elseif (expanded_nodes.nodes(i).num_matches < Parameters.RANSAC.min_pts_consensus),
                        leaf_flag(i) = 1;
                        expanded_nodes.nodes(i).leaf_flag = 1;
                    end
                end
                num_leaf_flag = sum(leaf_flag);
                if num_leaf_flag > 0,
                    leaf_flag = logical(leaf_flag);
                    Tree.nodes(Tree.num_nodes+1:Tree.num_nodes+num_leaf_flag) = expanded_nodes.nodes(leaf_flag); % adding to the tree
                    Tree.ind_leaf_nodes = [Tree.ind_leaf_nodes, Tree.num_nodes+1:Tree.num_nodes+num_leaf_flag];
                    ind_solution = find(solutions(leaf_flag));
                    Tree.solution_nodes = [Tree.solution_nodes, Tree.num_nodes+ind_solution];
                    Tree.num_nodes = Tree.num_nodes + num_leaf_flag;
                    Tree.nodes(ind_inner_nodes_parents(j)).children = [Tree.nodes(ind_inner_nodes_parents(j)).children, IDs_exp(leaf_flag)];
                    Tree.ind_IDs = [Tree.ind_IDs IDs_exp(leaf_flag)];
                    expanded_nodes.nodes(leaf_flag) = []; %removing the cluster that cannot be expanded.
                    IDs_exp(leaf_flag) = [];
                end
            end            
            %% 5. Updating the tree if ~isempty(expanded_nodes.nodes),
            num_inner = length(expanded_nodes.nodes);
            num_leaves = length(leaf_nodes.nodes);
            Tree.num_inner_nodes = Tree.num_inner_nodes + num_inner;
            % adding inner nodes
            Tree.nodes(Tree.num_nodes+1:Tree.num_nodes+num_inner) = expanded_nodes.nodes;
            Tree.ind_inner_nodes = [Tree.ind_inner_nodes, Tree.num_nodes+1:Tree.num_nodes+num_inner];
            Tree.nodes(ind_inner_nodes_parents(j)).children = [Tree.nodes(ind_inner_nodes_parents(j)).children, [expanded_nodes.nodes.ID]];
            Tree.ind_IDs = [Tree.ind_IDs [expanded_nodes.nodes.ID]]; % adding node's ID
            Tree.num_nodes = Tree.num_nodes + num_inner;
            Tree.depth = tree_depth;
            % adding leaf nodes
            Tree.nodes(Tree.num_nodes+1:Tree.num_nodes+num_leaves) = leaf_nodes.nodes;
            Tree.ind_leaf_nodes = [Tree.ind_leaf_nodes, Tree.num_nodes+1:Tree.num_nodes+num_leaves];
            Tree.num_nodes = Tree.num_nodes + num_leaves;
            Tree.ind_IDs = [Tree.ind_IDs [leaf_nodes.nodes.ID]];
            Tree.nodes(ind_inner_nodes_parents(j)).children = [Tree.nodes(ind_inner_nodes_parents(j)).children, [leaf_nodes.nodes.ID]];
        end % end for the inner nodes (variable j)
    end
end
%% 6. Top-Down phase
[Tree, Matches_Data,Final_Clusters] = Tree_Top_Down_phase(Tree, Matches_Data, Initial_Matches, Parameters, ind_repetitions, Images, Params_Display, bool_plot);
% %%
% Params_Display = Set_Display_Params();
if isempty(Final_Clusters),
   Solution = [];
   return;
end
% [Tree, Matches_Data, Final_Clusters] = RecoveryStrategy(Tree, Matches_Data, Initial_Matches, Parameters, ind_repetitions, Images, Params_Display, bool_plot);
%packing Solution
Final_Clusters.plot_type = 1;

%[Tree, Matches_Data, Final_Clusters] = HMA_Algorithm(Features, Parameters, Images, Initial_Matches, ratio_threshold, model)