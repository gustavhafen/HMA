%   File: 'Set_Initial_Matches.m'
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

function Initial_Matches = Set_Initial_Matches(feat_t, feat_q, Matching_options)
%function that finds the initial matches between the query image and the 
%training image by using the second-closest-neighbor ratio. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@parameters:
%feat_t := cell containing {Keypoints, Descriptors, dimensions}
%                      of the database images
%feat_q    := cell containing {Keypoints, Descriptors, dimensions}
%                      of the query image.
%Matching_options.ratio_threshold    := parameter to discriminate a false match with high
%                      probability. [Lowe04].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output: 
%initial_matching := cell array containing the indexes of the matched 
%                    features.  Every column correspond a pair:
%                    [query feature;training feature].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In particular this version is slower, but is designed to support very
%large images.
%%%NOTES: Originally was called neighbohrs_ratiov1_2_1 (and previously
%neighbors_ratio, version that support multiples 'textureless' images)%%%.

%% PREPROCESSING 
%all the descriptors, each column i is a descriptor vector of feature.
D_training = feat_t.Descriptors;
%descriptors, each column i is a descriptor vector of feature.
D_query = feat_q.Descriptors;
%number of features.
number_of_features_in_training = size(D_training,2);
number_of_features_in_query = size(D_query,2);

%Normalizing Descriptors (to make them robust to illumination)
for feature_i=1:number_of_features_in_training,
    D_training(:,feature_i) = D_training(:,feature_i)/norm(D_training(:,feature_i));
end


%% COMPUTING DISTANCES
%Compute the distances is equivalent to compute: d' = 1-<dq,d_ij> (for each
%query feature).
minimal_indexes = zeros(number_of_features_in_query,1);
condition_indexes = zeros(number_of_features_in_query,1);
for feature_i=1:number_of_features_in_query,
    %computing dot product for feature i.
    D_query_dot_proD_training =  D_query(:,feature_i)'/norm(D_query(:,feature_i))  * D_training;
    %computing the minimal in d' is equivalent to find the maximal <dquery_i,dtraining_ij>, it
    %means, we want to search for the maximal.
    [minimal_distance, minimal_idx] = max(D_query_dot_proD_training); %max works by column
    minimal_indexes(feature_i) = minimal_idx;
    D_query_dot_proD_training(1, minimal_idx) = -Inf;   %removing the minimal
    second_minimal_distance = max(D_query_dot_proD_training); %finding max's among the 'remaining' distances
    condition_indexes(feature_i,1) = sqrt(1-minimal_distance) < Matching_options.ratio_threshold .* sqrt(1 - second_minimal_distance); %formula for d'.
end

%% Discarding the unimpotrtant keypoints
query_indexes = 1:number_of_features_in_query;
%Packing the matches that fullfil the criteria.
initial_matching_indexes(1,:) = query_indexes(condition_indexes == 1); %first row to query 
initial_matching_indexes(2,:) = minimal_indexes(condition_indexes == 1)'; %second row to dabaseimage i
Initial_Matches.Query = feat_q.Keypoints(:,initial_matching_indexes(1,:));
Initial_Matches.Training = feat_t.Keypoints(:,initial_matching_indexes(2,:));
Initial_Matches.Original_indices = initial_matching_indexes;
Initial_Matches.num = size(Initial_Matches.Training,2);