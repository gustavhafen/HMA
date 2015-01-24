%   File: 'HMA_Local_Recovery_Phase.m'
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

function [corrected_nodes, Matches_Data] = HMA_Local_Recovery_Phase(expanded_nodes, Initial_Matches, Matches_Data, Params, ind_repetitions, Images, Params_Display, option_plot)
%Function that corrects the classes provided by the clustering. This is done
%by applying LDA to the estimated classes.

option_plot = 0;
corrected_nodes = expanded_nodes; %% copying the contents of the original nodes
if length(expanded_nodes.nodes) > 1,
    
    %% 0. Processing Cluster's Data.
    
    % a) Training set
    % indices of Inliers of all clusters in the current branch
    ind_inliers = [expanded_nodes.nodes.ind_inliers];
    % Inliers of all clusters (Query-keypoint positions)
    training_set = Initial_Matches.Query(1:2, ind_inliers)';
    %%
    % computing labels of the clusters
    label = Matches_Data.ID_clusters(ind_inliers); %[expanded_nodes.nodes(:).ID];
    num_expanded_nodes = length(expanded_nodes.nodes);
    % b) Test set
    
    % Indices of all points in the clusters (inliers and outliers)
    ind_clusters = [expanded_nodes.nodes(:).ind_matches];
    num_data_cluster = length(ind_clusters);
    % clusters' data
    matches_candidates.Query = Initial_Matches.Query(1:2, ind_clusters);
    matches_candidates.Training = Initial_Matches.Training(1:2, ind_clusters);
    %
    test_set = matches_candidates.Query';
    %%
    % clusters ID.
    ind_cluster_ID = [expanded_nodes.nodes(:).ID]; % previous labels
    num_cluster_classes = length(ind_cluster_ID); % num labels
    
    %% 1. Correction phase:
    classes_test = classify(test_set, training_set, label, 'diaglinear');
    %
    %% plotting
    if option_plot,
        hfig1 = figure;
        hfig2 = figure;
        width = Plot_Dual_Image(Images.Training, Images.Query,hfig1); hold on;
        Plot_Dual_Image(Images.Training, Images.Query,hfig2); hold on;
        for i=1:num_expanded_nodes,
            ind_aux = expanded_nodes.nodes(i).ind_matches;
            figure(hfig1);
            Plot_Initial_Matches_subset(Initial_Matches, ind_aux, [], list_of_colors(i), 'x',width);
            ind_aux = (classes_test == ind_cluster_ID(i));
            figure(hfig2);
            Plot_Initial_Matches_subset(matches_candidates, ind_aux, [], list_of_colors(i), 's',width);
        end
    end
    %% check
    ind_clusters_flaw = []; % wrong solutions
    ind_class_flaw = [];% indices of incorrect nodes
    ind_train_flaw = []; %indices of the matches in incorrect nodes
    ind_cluster2remove = zeros(1, num_cluster_classes);
    for i=1:num_cluster_classes, % The correction is just for the inner nodes
        %indices belonging to the class i
        class_i = ind_cluster_ID(i);
        ind_cluster_i = find(classes_test' == class_i); % matches of the new clustering
        if ~isempty(ind_cluster_i), % there are elements in the new cluster
            matches_cluster_i = []; % structure for RANSAC
            matches_cluster_i.Training = matches_candidates.Training(:, ind_cluster_i);
            matches_cluster_i.Query = matches_candidates.Query(:, ind_cluster_i);
            Mapping_i = corrected_nodes.nodes(i).Mapping;
            % Checking for matches agreeing with the previous model:
            [ind_consensus_clusters, aux, residuals] = Compute_Consesus_Model(Mapping_i, matches_cluster_i, Params.RANSAC.error_bound, ind_repetitions(ind_clusters(ind_cluster_i)));
            %
            if length(ind_consensus_clusters) < Params.RANSAC.min_pts_consensus,
                % The AT was wrong, reassigning the matches to closest clusters.
                ind_class_flaw = [ind_class_flaw, class_i];
                ind_clusters_flaw = [ind_clusters_flaw, ind_cluster_i];
                ind_train_flaw = [ind_train_flaw, find(label == class_i)'];
                ind_cluster2remove(i) = 1; % marked to remove
            else
                % Packing Solution
                cluster_i_elements = ind_clusters(ind_cluster_i);
                num_cluster_i_elements = length(cluster_i_elements);
                cluster_i_inliers = cluster_i_elements(ind_consensus_clusters);
                num_i_inliers = length(cluster_i_inliers);
                cluster_i_outliers = setdiff(cluster_i_elements, cluster_i_inliers);
                % updating the node's contents
                corrected_nodes.nodes(i).num_matches = num_cluster_i_elements;
                corrected_nodes.nodes(i).ind_matches = cluster_i_elements;
                corrected_nodes.nodes(i).residuals = residuals;
                corrected_nodes.nodes(i).ind_inliers = cluster_i_inliers;
                corrected_nodes.nodes(i).ind_outliers = cluster_i_outliers;
                corrected_nodes.nodes(i).num_inliers = num_i_inliers;
                corrected_nodes.nodes(i).ratio = num_i_inliers/num_cluster_i_elements;
                % Updating the Matches_Data matrix
                Matches_Data.ID_clusters(cluster_i_elements) = class_i;
                Matches_Data.rep_error(cluster_i_elements) = residuals;
                Matches_Data.label(cluster_i_outliers) = 0; %outliers
                Matches_Data.label(cluster_i_inliers)  = 1; %inliers
            end
        else
            ind_cluster2remove(i) = 1;
        end
    end
    % removing failure nodes
    corrected_nodes.nodes(logical(ind_cluster2remove)) = [];
    % updating the number of expanded nodes
    % removing the indices from the list of candidates (for the new clustering)
    ind_cluster_ID(logical(ind_cluster2remove)) = [];
    num_cluster_classes = length(ind_cluster_ID);
    % reassigning matches from defecting nodes
    %% Special cases
    if ~isempty(ind_class_flaw),
        test_set = test_set(ind_clusters_flaw, :); %inliers reassigned to this cluster
        % structure for the computation of the consensus
        matches_candidates.Training = matches_candidates.Training(:,ind_clusters_flaw);
        matches_candidates.Query = matches_candidates.Query(:,ind_clusters_flaw);
        %computing indices of elements to erase from training set
        ind_clusters_flaw_i = ind_clusters(ind_clusters_flaw);
        training_set(ind_train_flaw,:) = [];
        label(ind_train_flaw) = [];
        ind_cluster_classes = classify(test_set, training_set, label, 'diaglinear'); %classes of the flawed nodes
        %maybe add outliers
        % reassigning
        for i=1:num_cluster_classes,
            %indices belonging to the class i
            class_i = ind_cluster_ID(i);
            ind_cluster_i = find(ind_cluster_classes' == class_i );
            if ~isempty(ind_cluster_i), % there are elements in the new cluster
                % copying required elements
                matches_cluster_i = []; % structure for consensus
                matches_cluster_i.Training = matches_candidates.Training(:, ind_cluster_i);
                matches_cluster_i.Query = matches_candidates.Query(:, ind_cluster_i);
                Mapping_i = corrected_nodes.nodes(i).Mapping;
                % Checking for matches agreeing with the model:
                [ind_consensus_new, aux, residuals_new] = Compute_Consesus_Model(Mapping_i, matches_cluster_i, Params.RANSAC.error_bound, ind_repetitions(ind_clusters_flaw_i));
                % validating the new model
                cluster_new_elements = ind_clusters(ind_clusters_flaw(ind_cluster_i)); % elements to add to the cluster i
                if ~isempty(ind_consensus_new),
                    cluster_new_inliers = cluster_new_elements(ind_consensus_new);
                    cluster_new_outliers = setdiff(cluster_new_elements, cluster_new_inliers);
                    % Updating the node:
                    num_new_elements = length(cluster_new_elements);
                    num_new_inliers = length(cluster_new_inliers);
                    num_new_outliers = length(cluster_new_outliers);
                    corrected_nodes.nodes(i).ratio = (corrected_nodes.nodes(i).ratio*corrected_nodes.nodes(i).num_matches+ num_new_inliers)/(corrected_nodes.nodes(i).num_matches+ num_new_elements);
                    corrected_nodes.nodes(i).num_matches = corrected_nodes.nodes(i).num_matches+ num_new_elements;
                    corrected_nodes.nodes(i).ind_matches(end+1:end+num_new_elements) = cluster_new_elements;
                    corrected_nodes.nodes(i).residuals(end+1:end+num_new_elements) = residuals_new;
                    corrected_nodes.nodes(i).ind_inliers(end+1:end+num_new_inliers) = cluster_new_inliers;
                    corrected_nodes.nodes(i).ind_outliers(end+1:end+num_new_outliers) = cluster_new_outliers;
                    % Updating the Matches_Data matrix
                    Matches_Data.ID_clusters(cluster_new_elements) = class_i;
                    Matches_Data.rep_error(cluster_new_elements) = residuals_new;
                    Matches_Data.label(cluster_new_outliers) = 0; %outliers
                    Matches_Data.label(cluster_new_inliers)  = 1; %inliers
                else
                    cluster_new_outliers = cluster_new_elements;
                    Matches_Data.ID_clusters(cluster_new_elements) = class_i;
                    Matches_Data.rep_error(cluster_new_elements) = residuals_new;
                    Matches_Data.label(cluster_new_outliers) = 0; %outliers
                end
            end
        end
    end
    
    %% 4. Plotting
    if option_plot == 1,
        % plotting the leaf nodes
        Params_Display.Rep.Marker = '+';
        Params_Display.Rep.LineStyle = 'none';
        hfig =  figure;
        width = Plot_Dual_Image(Images.Training, Images.Query,hfig); hold on;        % plotting the cluster's matches
        for i=1:length(corrected_nodes.nodes),
            Plot_Initial_Matches_subset(Initial_Matches, expanded_nodes.nodes(i).ind_matches, [], list_of_colors(i), 'x',     width);
            if ~ind_cluster2remove(i),
                Plot_Initial_Matches_subset(Initial_Matches, corrected_nodes.nodes(i).ind_matches, [], list_of_colors(i), 'o', width);
                Plot_Initial_Matches_subset(Initial_Matches, corrected_nodes.nodes(i).ind_inliers, [], list_of_colors(i), '+', width);
            end
        end
    end
end