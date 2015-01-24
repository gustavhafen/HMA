%   File: 'RANSAC_Mapping_Model.m'
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

function [Max_Model, max_consensus_idxs, max_residuals] = RANSAC_Mapping_Model(matched_keypoints, parameters, Model, ind_repetitions)
%function that finds the affine transforms between dabatase image and
%query avoiding outliers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% I packed all the constrains in the following vector:
% RANSAC_parameters := [N_model_minimal, N_minimal_inliers, N_max_iterations, N_enought_inliers, error_threshold]
%   N_model_minimal   := minimal points needed for calculate the model
%   N_minimal_inliers := minimal points needed to accept the model
%   N_max_iterations  := maximal RANSAC's iterations
%   N_enought_inliers := minimal points needed to consider the model as a
%                      good estimation.
%   error_threshold   := error bound for RANSAC (real error is
%                        3*error_threshold).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@return:
%max_consensus_idxs := indexes of inliers of the best a.t. calculated in
%                      the space of images.
%Max_Model := the best Affine transformations (M = [A|T]).
%statistics := vector that contains in the first column the quadratic error
%              of the approximation(of the aff trans M). the iteration when
%              the solution was find is in the second column and the third
%              has the total number of iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_size = size(matched_keypoints.Query,2);%number of features
acceptable_consensus = max(parameters.min_pts_consensus, floor((1-parameters.eps)*sample_size));
[A_full,b_full] = Model.Build_Matrix(matched_keypoints);%complete matrix A and b of Ax=b problem to find best a.t.
%% generating combinations
combinations = my_random_s(sample_size,parameters.min_pts_model, parameters.N_max_it);
% removing combinations with repetitions
 comb = ind_repetitions(combinations)';
 ind_rep_comb = (comb(:,1) == comb(:,2))|(comb(:,1) == comb(:,3))|(comb(:,2) == comb(:,3));
 combinations(:,ind_rep_comb) = [];
%
max_consensus = parameters.min_pts_consensus - 1; %parameters.min_pts_consensus - 1;%flag to know if there is a combination with at least threshold # of points
Max_Model = [];
max_consensus_idxs = [];
max_residuals = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~isempty(combinations)
    sample = combinations(:,1);  %we have our sample
    combinations(:,1) = []; %removing the sampled combination
    %for each we calculate the AT
    A(1:parameters.min_pts_model,:) = A_full(sample, :);
    A(parameters.min_pts_model+1:2*parameters.min_pts_model,:) = A_full(sample_size+sample,:);
    b(1:parameters.min_pts_model,:) = b_full(sample, :);
    b(parameters.min_pts_model+1:2*parameters.min_pts_model,:) = b_full(sample_size+sample,:);
    Model_Solution = Model.Estimation(A, b);
    if isempty(Model_Solution), %
        continue;
    end
    %consensus
    [idx, consensus, res] = Compute_Consesus_Model(Model_Solution, matched_keypoints, parameters.error_bound,ind_repetitions);
         if consensus > max_consensus,
            max_consensus = consensus;
            max_consensus_idxs = idx;
            Max_Model = Model_Solution;
            max_residuals = res;
            if (consensus > acceptable_consensus) %termination condition
                break;
            end
         end
end
%if the model_solution satisfy the minimal number of points:
if max_consensus > parameters.min_pts_consensus - 1,    
    %re-estimating the best Trasformation with all points
    max_Sample = length(max_consensus_idxs);
    A(1:max_Sample,:) = A_full(max_consensus_idxs, :);
    A(max_Sample+1:2*max_Sample,:) = A_full(sample_size+max_consensus_idxs,:);
    b(1:max_Sample,:) = b_full(max_consensus_idxs, :);
    b(max_Sample+1:2*max_Sample,:) = b_full(sample_size+max_consensus_idxs,:);
%
    Model_Solution = Model.Estimation(A, b);
    [idx, consensus, res] = Compute_Consesus_Model(Model_Solution, matched_keypoints, parameters.error_bound,ind_repetitions);
%
    if consensus >= max_consensus,
        Max_Model = Model_Solution;
        max_consensus_idxs = idx;
        max_residuals = res;
    end
end