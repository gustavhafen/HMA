%   File: 'Compute_Consesus_Model.m'
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

function [inliers_idx, consensus, errors] = Compute_Consesus_Model(Model, matched_keypoints, error_bound, ind_repetitions)
%compute the number of inliers inside error bounds. Error is defined as
%scale_bounds * pixel (=how many pixels of error tolerance) thresholds
    
num_kpts = size(matched_keypoints.Query,2);
%%estimations
Rep_Query_pts = Model.Transformation * [matched_keypoints.Training(1:2,:); ones(1,num_kpts)];
Rep_Training_pts = Model.Transformation_inv * [matched_keypoints.Query(1:2,:); ones(1,num_kpts)];
Rep_Query_pts = Rep_Query_pts./(ones(3,1)*Rep_Query_pts(3,:));
Rep_Training_pts = Rep_Training_pts./(ones(3,1)*Rep_Training_pts(3,:));

%% sym error
Diff_Rep_Query = matched_keypoints.Query(1:2,:) - Rep_Query_pts(1:2,:);
Diff_Rep_Training = matched_keypoints.Training(1:2,:) - Rep_Training_pts(1:2,:);
%error_vectors = [Diff_Rep_Training + Diff_Rep_Query];
errors = 0.5*sqrt(sum(Diff_Rep_Query .* Diff_Rep_Query))+0.5*sqrt(sum(Diff_Rep_Training .* Diff_Rep_Training));
%%indexes
inliers_idx = (errors < 3*error_bound*ones(1,num_kpts));

vector = ind_repetitions(inliers_idx);
vector = sort(vector);
consensus = sum(vector(2:end) - vector(1:end-1) > 0)+1;
inliers_idx = find(inliers_idx); % with repetitions