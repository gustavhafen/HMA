%   File: 'Plot_Matches_Model.m'
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

function Plot_Matches_Model(pts_query, pts_training, width, Style_Matches, Style_Rep, fig_handle)
%function that plot the image (if both cell values are not empty), plot the
%line showing the matched keypoints, and the projected-error lines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A := The Affine matrix transform.
% T := The Affine translation vector.
%points_query := locations of the query keypoints.
%points_db := locations of the (correswponding)  db keypoints.
%Images   := cell with the images. {I_l I_r}, where I_l correspond to the
%            query image and I_r is the database image. If Images{2} is
%            empty the images are not plotted neither.
%colors := cell containing the colors of the plots on the following order: 
%          {col_matches_line col_error_lines}, 
%          where col refers to the color, style line, and marker of each 
%          line, for example col1 = 'r-.' produces a red segment of a line
%          with the begining point and ending point  marked with a point.
%          col_matches_line = color of lines showing the matched
%                             keypoints.
%          col_error_line = color of lines that unite the projected db
%                           keypoint with the corresponding keypoint.
%sbool, = 1 color, 0 gray.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No return:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the images
plot_style_points = [Style_Rep.Color, Style_Rep.Marker];
plot_style_lines = [Style_Rep.Color, Style_Rep.LineStyle];
%% plotting
%number of points to be plotted
num_pts = size(pts_training,2);
%%Paper Style
num_images = 1;
if nargin == 7,
    num_images = length(fig_handle);
end
if num_images == 1,
    line([width + pts_query(1,:), pts_training(1,:)],[pts_query(2,:), pts_training(2,:)],'Marker','o','LineStyle',Style_Matches.LineStyle,'Color','k', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
%     line([width + pts_query(1,:); pts_training(1,:)],[pts_query(2,:); pts_training(2,:)],'Marker','o','LineStyle',':','Color','y', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
else
    figure(fig_handle(1)),
    plot(pts_training(1,:), pts_training(2,:),'Marker','o','LineStyle',Style_Matches.LineStyle,'Color',':', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
    figure(fig_handle(2)),
    plot(width + pts_query(1,:), pts_query(2,:),'Marker','o','LineStyle',Style_Matches.LineStyle,'Color','k', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
end
