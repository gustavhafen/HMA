%   File: 'RecoveryStrategy.m'
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

function [Tree, Matches_Data, Final_Clusters] = RecoveryStrategy(Tree, Matches_Data, Initial_Matches, Params, ind_repetitions, Images, Params_Display, option_plot)

limbot = 25;
minimal_limbo = 2;
num_imported_inliers = 6;
option_plot = 0;

%% TOP-DOWN phase

%% 0. Processing Cluster's Data.
% training set
% indices of Inliers of all clusters
ind_training = (Matches_Data.label == 1);
% Inliers of all clusters (Query-keypoint positions)
training_set = Initial_Matches.Query(1:2, ind_training)';
% test set
ind_test = 1:Initial_Matches.num;
ind_test(ind_training) = [];
num_training = sum(ind_training);
num_test = Initial_Matches.num - num_training;
test_set = Initial_Matches.Query(1:2, ind_test)';
% computing labels of the clusters
labels = Matches_Data.ID_clusters(ind_training);

%% 1. Correction phase:
%classes_test = classify(test_set, training_set, labels, 'linear');
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
ind_classes = [Tree.nodes(ind_solution).ID];
num_clusters = length(ind_solution);
New_Clusters = [];
ind_new_cluster = 0;
for i=1:num_clusters,
    Final_Clusters.Mapping(i) = Tree.nodes(ind_solution(i)).Mapping;
    Final_Clusters.Inliers(i).Training = Initial_Matches.Training(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).Query = Initial_Matches.Query(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).ind = Tree.nodes(ind_solution(i)).ind_inliers;
    
    Matched_Points.Query = Final_Clusters.Inliers(i).Query(1:2,:);
    [pts_region, BB] = discard_repeated_rows_matrix(Matched_Points.Query');
    [region_polygon, BB] = Graham_Scan_batch(pts_region');
    Final_Clusters.Inliers(i).polygon = region_polygon;
    Final_Clusters.Inliers(i).centroid = mean(pts_region);
end
Final_Clusters.plot_type= 1;
%
if option_plot == 1,
    Plot_Multi_Model(Final_Clusters, Images,Params_Display, []); hold on;
    title('Before AT Recovering')
end

new_ind_inliers = []; %indices of the added inliers
for i=1:num_clusters,
    %indices belonging to the class i
    ind_node_i = ind_solution(i);
    class_i = ind_classes(i);
    ind_cluster_i = (classes_test == class_i);
    if sum(ind_cluster_i)>0,
        ind_test_i = ind_test(ind_cluster_i); %real indices
        matches = []; % structure for RANSAC
        matches.Training = Initial_Matches.Training(:, ind_test_i);
        matches.Query = Initial_Matches.Query(:, ind_test_i);
        [ind_consensus, aux, residuals_i] = Compute_Consesus_Model(Tree.nodes(ind_node_i).Mapping, matches, Params.RANSAC.error_bound, ind_repetitions(ind_test_i));
        %% 3.Recovery phase: (adding them!!!)
        % - points agreeing with the model
        if ~isempty(ind_consensus),
            num_new_elements = length(ind_consensus);
            ratio = ( Tree.nodes(ind_node_i).ratio*Tree.nodes(ind_node_i).num_matches + num_new_elements )/ (num_new_elements + Tree.nodes(ind_node_i).num_matches);
            Tree.nodes(ind_node_i).ind_matches = [Tree.nodes(ind_node_i).ind_matches, ind_test_i(ind_consensus)];
            Tree.nodes(ind_node_i).ind_inliers = [Tree.nodes(ind_node_i).ind_inliers, ind_test_i(ind_consensus)];
            Tree.nodes(ind_node_i).num_matches = Tree.nodes(ind_node_i).num_matches + num_new_elements;
            Tree.nodes(ind_node_i).residuals(end+1:end+num_new_elements) = residuals_i(ind_consensus);
            Tree.nodes(ind_node_i).ratio_inliers = ratio;
            %
            Matches_Data.ind_cluster(ind_test_i(ind_consensus)) = class_i;
            Matches_Data.rep_error(ind_test_i(ind_consensus)) = residuals_i(ind_consensus);
            Matches_Data.label(ind_test_i(ind_consensus)) = 1;    %inliers
            Matches_Data.tree_level(ind_test_i(ind_consensus)) = Tree.nodes(ind_node_i).depth;
            new_ind_inliers = [new_ind_inliers , ind_test_i(ind_consensus)];
        end
        
        %     %remaining points
        ind_limbo = find((residuals_i > Params.RANSAC.error_bound) & (residuals_i < limbot)); % limbo points\
        %points outside the polygon
        ind_limbo_disc = inpolygon(matches.Query(1,ind_limbo),matches.Query(2,ind_limbo),Final_Clusters.Inliers(i).polygon(1,:),Final_Clusters.Inliers(i).polygon(2,:)); %check if points fall inside the polygon
        % removing them
        ind_limbo(ind_limbo_disc) = []; %limbo indices (relative to the residuals (not the real indices))
        %% HERE PCA
        %     if ~isempty(ind_limbo)&&(length(ind_limbo)>=minimal_limbo),
        limbopts_all = matches.Query(1:2,ind_limbo);
        ind_real_limbo = ind_test_i(ind_limbo); %% new real ind_limbo indices (wrt to the init mat)
        % plotting
        if option_plot,
            fig_handle = figure;
            width = Plot_Dual_Image(Images.Training, Images.Query,fig_handle); hold on;%width of the query image (Image on the left)
            plot(Final_Clusters.Inliers(i).Training(1,:), Final_Clusters.Inliers(i).Training(2,:), 'cx');
            plot(matches.Training(1,ind_limbo), matches.Training(2,ind_limbo), 'y+');
            plot(width + Final_Clusters.Inliers(i).Query(1,:), Final_Clusters.Inliers(i).Query(2,:), 'cx');
            plot(width + matches.Query(1,ind_limbo), matches.Query(2,ind_limbo), 'y+');
            pol_i = Final_Clusters.Inliers(i).polygon;
            pol_i(1,:) = width + pol_i(1,:);
            draw_polygon(pol_i);
        end
        symb = 'ro';
        %% Here  a  while over the possible directions!!
        if (length(limbopts_all)>=minimal_limbo),
            %        %
            if option_plot,
                plot(matches.Training(1,ind_limbo), matches.Training(2,ind_limbo), symb);
                plot(width + matches.Query(1,ind_limbo), matches.Query(2,ind_limbo), symb);
            end
            % looking for closest neighbors in query image
            
            centroid = mean(limbopts_all, 2);
            inliers = Final_Clusters.Inliers(i).Query(1:2,:);
            distances = dist([centroid, inliers]);
            [mute, closest_neighbors] = sort(distances(1,2:end)); % computing the distances from center to the inliers
            imported_inliers = closest_neighbors(1:num_imported_inliers);
            %% new
            ind_real_imported_inliers = Final_Clusters.Inliers(i).ind(closest_neighbors(1:num_imported_inliers)); %%
            if option_plot,
                plot(Final_Clusters.Inliers(i).Training(1,imported_inliers), Final_Clusters.Inliers(i).Training(2,imported_inliers), 'co');
                plot(width + Final_Clusters.Inliers(i).Query(1,imported_inliers), Final_Clusters.Inliers(i).Query(2,imported_inliers), 'co');
            end
            % fitting an AT
            matches_rec = []; % structure for RANSAC
            matches_rec.num = length(ind_limbo) + num_imported_inliers;
            matches_rec.Training = [Final_Clusters.Inliers(i).Training(:, imported_inliers), matches.Training(:,ind_limbo)];% adding the imported inliers to the limbo points
            matches_rec.Query = [Final_Clusters.Inliers(i).Query(:, imported_inliers), matches.Query(:,ind_limbo)];% adding the imported inliers to the limbo points
            Mapping_i = Final_Clusters.Mapping(i);
            % Checking for matches agreeing with the model:
            aux = [ ind_repetitions(Tree.nodes(ind_solution(i)).ind_inliers(imported_inliers)) ,ind_repetitions(ind_test_i(ind_limbo)) ];
            [AT_limbo, ind_limbo_inl, reslimbo] = RANSAC_Mapping_Model(matches_rec, Params.RANSAC, Params.Model, aux);
            if ~isempty(AT_limbo)&&(sum(ind_limbo_inl > num_imported_inliers)) %matches that haven't been detected.
                % ind_limbo_inl = indices relative to ind_l_p, or to matches_rec
                if option_plot,
                    ind_support = ind_limbo_inl(ind_limbo_inl <= num_imported_inliers);
                    plot(matches_rec.Training(1,ind_support), matches_rec.Training(2,ind_support), 'm*');
                    plot(width + matches_rec.Query(1,ind_support), matches_rec.Query(2,ind_support), 'm*');
                end
                %
                aux_ind = ind_limbo_inl(ind_limbo_inl > num_imported_inliers) - num_imported_inliers;
                ind_real_recovered = ind_real_limbo(aux_ind);
                ind_supporting = ind_real_imported_inliers(ind_limbo_inl < num_imported_inliers);
                ind_recovered = ind_limbo_inl;
                ind_new_cluster = ind_new_cluster + 1;                
                %
                if option_plot,
                    plot(matches_rec.Training(1,ind_recovered), matches_rec.Training(2,ind_recovered), 'g*');
                    plot(width + matches_rec.Query(1,ind_recovered), matches_rec.Query(2,ind_recovered), 'g*');
                end
                New_Clusters.Mapping(ind_new_cluster) = AT_limbo;
                New_Clusters.Inliers(ind_new_cluster).Training = matches_rec.Training(:,ind_recovered); %Initial_Matches.Training(:, );
                New_Clusters.Inliers(ind_new_cluster).Query = matches_rec.Query(:, ind_recovered);
                %                 new_ind_inliers = [new_ind_inliers , aux(ind_recovered)];
                % imported inliers
                New_Clusters.Inliers(ind_new_cluster).Query = matches_rec.Query(:, ind_recovered);
                % new indices:
                New_Clusters.Inliers(ind_new_cluster).ind = ind_real_recovered;
                New_Clusters.Inliers(ind_new_cluster).supporting = ind_supporting;
                New_Clusters.cluster_index(ind_new_cluster) = i;
                % for now not implemented for Tree and Matches Structure!!
                %             num_new_elements = length(ind_consensus);
                %             ratio = ( Tree.nodes(ind_node_i).ratio_inliers*Tree.nodes(ind_node_i).num + num_new_elements )/ (num_new_elements + Tree.nodes(ind_node_i).num);
                %             Tree.nodes(ind_node_i).ind_matches = [Tree.nodes(ind_node_i).ind_matches, ind_test_i(ind_consensus)];
                %             Tree.nodes(ind_node_i).ind_inliers = [Tree.nodes(ind_node_i).ind_inliers, ind_test_i(ind_consensus)];
                %             Tree.nodes(ind_node_i).num = Tree.nodes(ind_node_i).num + num_new_elements;
                %             Tree.nodes(ind_node_i).residuals(end+1:end+num_new_elements) = residuals_i(ind_consensus);
                %             Tree.nodes(ind_node_i).ratio_inliers = ratio;
                %             %
                %             Matches_Data.ind_cluster(ind_test_i(ind_consensus)) = class_i;
                %             Matches_Data.rep_error(ind_test_i(ind_consensus)) = residuals_i(ind_consensus);
                %             Matches_Data.label(ind_test_i(ind_consensus)) = 1;    %inliers
                %             Matches_Data.tree_level(ind_test_i(ind_consensus)) = Tree.nodes(ind_node_i).depth;
            end
            %             end
        end
    end
end

clear Final_Clusters;
% old clusters
for i=1:num_clusters,
    Final_Clusters.Mapping(i) = Tree.nodes(ind_solution(i)).Mapping;
    Final_Clusters.Inliers(i).Training = Initial_Matches.Training(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).Query = Initial_Matches.Query(:, Tree.nodes(ind_solution(i)).ind_inliers);
    Final_Clusters.Inliers(i).ind = Tree.nodes(ind_solution(i)).ind_inliers;
    Final_Clusters.Inliers(i+num_clusters).supporting = [];
end

%
for i = 1:ind_new_cluster,
    if New_Clusters.cluster_index(i) ~= 0,
        Final_Clusters.Inliers(New_Clusters.cluster_index(i)).Training =  [Final_Clusters.Inliers(New_Clusters.cluster_index(i)).Training, New_Clusters.Inliers(i).Training];
        Final_Clusters.Inliers(New_Clusters.cluster_index(i)).Query =  [Final_Clusters.Inliers(New_Clusters.cluster_index(i)).Query, New_Clusters.Inliers(i).Query];
    end
end

Final_Clusters.plot_type= 1;
if option_plot == 1,
    Plot_Multi_Model(Final_Clusters, Images,Params_Display, []); hold on;
    title('After Recovering')
end