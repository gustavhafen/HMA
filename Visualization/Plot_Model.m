%   File: 'Plot_Model.m'
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

function fig_handle = Plot_Model(Model, indices, matched_keypoints, Images, Visual_Options, fig_handle, aux)
%function that visualize (ONE Image) with the matched_keypoints(spicified 
%in indices vector), the polygons and the projected points of the 
%corresponding database keypoints, and the projected db polygon , using the
%colors specified in colots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M := The Affine transform (M = [A|T]).
%indices := indices of projected points in the image. (could be the points
%           used to obtain the A.T. or the inliers in the error bounds).
%matched_keypoints := matched Keypoints of query and database
%Images   := cell with the images. {I_l I_r}, where I_l correspond to the
%            query image and I_r is the database image.
%polygons := cell containing the vertices of both polygons {pol_l, pol_r}
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

%% Plotting Images

if nargin == 7, 
    if Visual_Options.Style == 1, %One image mode
        width = size(Images.Training, 2);
    else %% Two images mode
        width = 0;
    end
else
    if Visual_Options.Style == 1, %One image mode
        figure(fig_handle);
        width = Plot_Dual_Image(Images.Training, Images.Query,fig_handle);%width of the query image (Image on the left)
        hold on;
    else %% Two images mode
        figure(fig_handle(1));
        imshow(Images.Training); hold on;    
        figure(fig_handle(2));
        imshow(Images.Query); hold on;
        width = 0;
    end
end

%% Reading points
if ~isempty(matched_keypoints),
    %number of lines plotted.
    if isempty(indices),
        n = size(matched_keypoints.Query,2);
        indices = 1:n;
    else
        n = length(indices);
    end
    points_Query    = matched_keypoints.Query(1:2,indices);
    points_Training = matched_keypoints.Training(1:2,indices);
else
    n = 0; %no points to plot
end
    
%% Drawing lines
if (n > 0) %&& (~isempty(Model.Transformation)), %there are lines to plot
    if ~isempty(Model)&&~isempty(Model.Transformation),%modify later version with reprojections
        if Visual_Options.Style == 1,
            Plot_Reprojected_Matches_Model(Model.Transformation,points_Query,points_Training, width, Visual_Options.Matches,Visual_Options.Rep);%plotting the image and the lines(matching lines and projected error lines)
        else
            Plot_Reprojected_Matches_Model(Model.Transformation,points_Query,points_Training, width, Visual_Options.Matches,Visual_Options.Rep,fig_handle);%plotting the image and the lines(matching lines and projected error lines)
        end
    else 
        if Visual_Options.Style == 1,
            Plot_Matches_Model(points_Query,points_Training, width, Visual_Options.Matches,Visual_Options.Rep);%plotting the image and the lines(matching lines and projected error lines)
        else
            Plot_Matches_Model(points_Query,points_Training, width, Visual_Options.Matches,Visual_Options.Rep,fig_handle);%plotting the image and the lines(matching lines and projected error lines)
        end
    end
end

if ~isempty(Images.polygon_training),
%drawing training-projected polygon:
    training_pol = Images.polygon_training';
    if (~isempty(training_pol))&&(~isempty(Model.Transformation)), 
        %number of vertices
        proj_pol = Images.polygon_query;
        proj_pol(1,:) = width + proj_pol(1,:);
        %% %%
        %drawing the projected polygon
        if Visual_Options.Style == 1, %One image mode        
            draw_polygon(proj_pol, Visual_Options.Pol_Rep);
            draw_polygon(training_pol, Visual_Options.Pol_Training);
        else
            figure(fig_handle(2));
            draw_polygon(proj_pol, Visual_Options.Pol_Rep);
            figure(fig_handle(1));
            draw_polygon(training_pol, Visual_Options.Pol_Training);
        end
    end
end