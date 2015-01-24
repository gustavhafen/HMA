%   File: 'README.m'
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


README

Hierachical Multi-Affine (HMA) Toolbox v1.0 
Submitted to Transactions of Medical Imaing (TMI)
May 30th, 2012.
G. Puerto and G.L. Mariottini
ASTRA Lab
University of Texas at Arlington

1. Requirements:

The toolbox requires the VL_Feat Library, available at: 
http://www.vlfeat.org/
The SURF features were extracted by using the Dirk-Jan Kroon's implementation, available at:
http://www.mathworks.com/matlabcentral/fileexchange/28300

This release of the toolbox also includes a demo which shows HMA's feature matching with a non-planar object.

The demo.mat file includes an image pair (Before, After), the set of extracted SIFT features (Features), and the set of initial matches (Initial_Matches) computed with nearest neighbor distance ratio technique.

2. The MIS Database:

The MIS database is a collection of MAT files (each mat file represents an image pair) which includes:

*An image pair (Before, After). Note the images are named Before and After because they were extracted Before and After an occlusion, or other MIS event.
*Type of matched features (Feat_Type).
*A structure containing the keypoints and descriptors of the extracted features for both images (Features).
*The mapping ground-truth: the manually selected pairs of corresponding points per each image (Ground_truth_Mapping).
*The matching ground-truth (SIFT only): the manually selected labels of the matches divided into two fields: true and false (these fields contain the indices of the true/false matches wrt the Initial_Matches)
*The set of initial matches computed with nearest neighbor distance ratio technique (Initial_Matches).
*The delimiting polygon of the kidney in the Training image.
 
***

To visualize the set of corresponding points run the following code:

	fig_handle = 1;
	width = Plot_Dual_Image(Before, After, fig_handle);
	X = [width + Ground_truth_Mapping.Query(1,:)', Ground_truth_Mapping.Training(1,:)']';
	Y = [Ground_truth_Mapping.Query(2,:)', Ground_truth_Mapping.Training(2,:)']';
	line(X,Y,'Marker','+','LineStyle',':','Color','y', 'LineWidth', 0.5, 'MarkerSize', 8);

***

To visualize the set of corresponding points run the following code:

	fig_handle = 2;
	width = Plot_Dual_Image(Before, After, fig_handle);
	% true (correct) matches
	X_T = [width + Initial_Matches.Query(1,Ground_truth_Matches.True)', Initial_Matches.Training(1,Ground_truth_Matches.True)']';
	Y_T = [Initial_Matches.Query(2,Ground_truth_Matches.True)', Initial_Matches.Training(2,Ground_truth_Matches.True)']';
	line(X_T,Y_T,'Marker','+','LineStyle',':','Color','g', 'LineWidth', 0.5, 'MarkerSize', 8);
	% false (wrong) matches
	X_F = [width + Initial_Matches.Query(1,Ground_truth_Matches.False)', Initial_Matches.Training(1,Ground_truth_Matches.False)']';
	Y_F = [Initial_Matches.Query(2,Ground_truth_Matches.False)', Initial_Matches.Training(2,Ground_truth_Matches.False)']';
	line(X_F,Y_F,'Marker','+','LineStyle',':','Color','r', 'LineWidth', 0.5, 'MarkerSize', 8);

3. Running an instance of the database (e.g., MIS_SIFT_1.mat)

Run the following code:

    %% preprocessing:
    load('MIS_SIFT_1.mat'); % loading database
    % parameters:
    Images = prepare_Images(Before, After, Polygons)
    Params = Set_HMA_parameters('SIFT');
    % Stop condition (% oof inliers) at each node level: 
    ratio_threshold = 0.9;
    % model (affine only)
    model = 'A';
    % HMA
    Solution = HMA_Algorithm(Features, Params, Images, Initial_Matches, ratio_threshold, model);
    % Visualization
    Display_Params = Set_Display_Params();
    Plot_Multi_Model(Solution, Images, Display_Params, []);