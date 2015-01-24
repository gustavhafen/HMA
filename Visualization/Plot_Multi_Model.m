%   File: 'Plot_Multi_Model.m'
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

function fig_handle = Plot_Multi_Model(Solution, Images,Visual_Options, fig_handle)

%function that visualize (ONE Image) with the matched_keypoints(spicified
%in indexes vector), the polygons and the projected points of the
%corresponding database keypoints, and the projected db polygon , using the
%colors specified in colots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M := The Affine transform (M = [A|T]).
%indexes := indexes of projected points in the image. (could be the points
%           used to obtain the A.T. or the inliers in the error bounds).
%matched_keypoints := matched Keypoints of query and database
%Images   := cell with the images. {I_l I_r}, where I_l correspond to the
%            query image and I_r is the database image.
%polygons := cell containing the vertexes of both polygons {pol_l, pol_r}
%            pol_l correspond to the polygon of the ground truth(query
%            image).
%            pol_r corresponds to the polygon from the database image(will
%            be mapped to the query image).
%colors := cell containing the colors of the plots on the following order:
%          {col_pol_query col_proj_pol col_matches_line col_error_lines},
%          where col refers to the color, style line, and marker of each
%          line, for example col1 = 'r-.' produces a red segment of a line
%          with the begining point and ending point  marked with a point.
%          col_pol_query = color of ground truth polygon.
%          col_pol_db = color of projected db polygon.
%          col_matches_line = color of lines showing the matched
%                             keypoints.
%          col_error_line = color of lines that unite the projected db
%                           keypoint with the corresponding keypoint.
%sbool, = 1 color, 0 gray.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No return:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(fig_handle),
    if Visual_Options.Style == 1,
        fig_handle=figure;
    else
        fig_handle(1)=figure;
        fig_handle(2)=figure;
    end
end
%%
if Solution.plot_type == 1,
    num = length(Solution.Mapping);
    for i=1:num,
        %% Selecting points: %%fixed
        %% Generating the convex-hull polygon for the matches points
        Matched_Points.Query = Solution.Inliers(i).Query(1:2, :);
        Matched_Points.Training = Solution.Inliers(i).Training(1:2,:);        
        [pts_region, ind_red] = discard_repeated_rows_matrix(Matched_Points.Training');
        [region_polygon, ind_ch] = Graham_Scan_batch(pts_region');
        %% NEW %%
        
        ind_query = ind_red(ind_ch);
        polygon_query = Matched_Points.Query(:, ind_query);
        %% %%
        Images.polygon_training = region_polygon';
        Images.polygon_query = polygon_query;
        Visual_Options.Pol_Training.Color = list_of_colors(i);
        Visual_Options.Pol_Rep.Color = list_of_colors(i);
        Visual_Options.Matches.Color = list_of_colors(i);
        
        if i==1, %the first time we plot the gt.
            Plot_Model(Solution.Mapping(i), [], Matched_Points, Images, Visual_Options,fig_handle);
        else
            Images.polygon_query_gt = []; % removing ground truth
            Plot_Model(Solution.Mapping(i), [], Matched_Points, Images, Visual_Options,fig_handle,1);
        end
    end
else
    %%style 2
    polygon_training = Images.polygon_training;
    Images.polygon_training = [];
    Images.polygon_query_gt = [];
    Matched_Points.Query = Solution.Inliers.Query(1:2, :);
    Matched_Points.Training = Solution.Inliers.Training(1:2,:);
    Visual_Options.Pol_Training.Color = 'r';
    Visual_Options.Pol_Rep.Color = 'r';
    Visual_Options.Matches.Color = 'y';
    Plot_Model([], [], Matched_Points, Images, Visual_Options,fig_handle);
    Visual_Options.Pol_Training.Color = 'r';
    draw_polygon(polygon_training', Visual_Options.Pol_Training);
end