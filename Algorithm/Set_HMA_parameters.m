%   File: 'Set_HMA_parameters.m'
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

function Params = Set_HMA_parameters(Feature_type, transformation_model, ratio_threshold)
% Function that sets HMA options.
%%
%@ params:
% Feature Type: 'SIFT', 'SURF', or 'ASIFT'
% transformation model: 'A' (Affine)
% ratio threshold: Initial matching ratio threshold (In case of the initial matches are not provided).

if nargin < 1,
    Params = [];
    display('error');
end
if nargin < 3,
    transformation_model = 'A';
end
if nargin < 4,
    if transformation_model == 'A',
        Model.Estimation = @Estimation_Affine;
        Model.Build_Matrix = @Build_Matrix_Affine;
    else %homography
        Model.Estimation = @Estimation_Homography;
        Model.Build_Matrix = @Build_Matrix_Homography;
    end
    %
    Matching.ratio_threshold = 0.8;
end
if nargin < 5,
    % RANSAC's parameters
    RANSAC.prob = 0.999; %confidence of results
    RANSAC.eps = 0.4; %probability of a point is an outlier
    if strcmp(transformation_model, 'A'),
        RANSAC.min_pts_model = 3; %minimal # of points to compute the model
    else
        RANSAC.min_pts_model = 4; %minimal # of points to compute the model
    end
    RANSAC.Model = transformation_model; %Homography (H) or Affine (A)
    
    %minimal # of points to accept the model
    if strcmp(Feature_type, 'ASIFT'),
        RANSAC.min_pts_consensus = 10;
        RANSAC.prob = 0.999999999; %confidence of results
    else
        RANSAC.min_pts_consensus = 6;
    end
        RANSAC.N_max_it = compute_N_threshold_RANSAC(RANSAC.prob, RANSAC.eps, RANSAC.min_pts_model); %maximal number of iterations of RANSAC
    % RANSAC's thresholds
    % indicates "t" error bound for the RANSAC in the image (euclidean
    % distance) - RANSAC will use "3*t" (pixel Euclidean distance)
    RANSAC.error_bound = 5/3;
    %RANSAC.N_enought_inliers = compute_N_threshold_RANSAC(RANSAC.prob, RANSAC.eps, RANSAC.min_pts_consensus);
    
end
if nargin < 6,
    Clustering.ratio_threshold = 0.9;
    Clustering.limbo_factor_threshold = 10;
%     Clustering.min_ERROR = 8*3*RANSAC.error_bound;
    Clustering.min_ERROR = Inf;%7*3*RANSAC.error_bound; 
    Clustering.max_ERROR = 0;%2.5*Clustering.min_ERROR;
    Clustering.factor_IQR = 1.5;
    if strcmp(Feature_type, 'ASIFT'),
        Clustering.max_num_inliers = 10*RANSAC.min_pts_consensus;
    else
        Clustering.max_num_inliers = 3*RANSAC.min_pts_consensus;
    end
    %     HKC.recovery_ratio_threshold = 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Packing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Matching = Matching;
Params.Model = Model;
Params.RANSAC = RANSAC;
Params.Clustering = Clustering;