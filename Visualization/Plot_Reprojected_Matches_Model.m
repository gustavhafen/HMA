%   File: 'Plot_Reprojected_Matches_Model.m'
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

function Plot_Reprojected_Matches_Model(Transformation, pts_query, pts_training, width, Style_Matches, Style_Rep, fig_handle)
%function that plot the image (if both cell values are not empty), plot the
%line showing the matched keypoints, and the projected-error lines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No return:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the images
plot_style_points = [Style_Rep.Color, Style_Rep.Marker];
plot_style_lines = [Style_Rep.Color, Style_Rep.LineStyle];
%% plotting
%number of points to be plotted
num_pts = size(pts_training,2);
projected_pts = Transformation*[ pts_training; ones(1,num_pts)];
projected_pts = projected_pts./(ones(3,1)*projected_pts(3,:));
%%
projected_t_pts = Transformation\[ pts_query; ones(1,num_pts)];
projected_t_pts = projected_t_pts ./(ones(3,1)*projected_t_pts (3,:));

 %%Paper Style
num_images = 1;
if nargin == 7,
    num_images = length(fig_handle);
end

if num_images == 1,
% style 2 (multiple lines)
    X = [width + pts_query(1,:)', pts_training(1,:)']';
    Y = [pts_query(2,:)', pts_training(2,:)']';
   line(X,Y,'Marker','o','LineStyle',Style_Matches.LineStyle,'Color','k', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
else
    figure(fig_handle(1)),
    plot(pts_training(1,:), pts_training(2,:),'Marker','o','LineStyle',Style_Matches.LineStyle,'Color','k', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
    figure(fig_handle(2)),
    plot(width + pts_query(1,:), pts_query(2,:),'Marker','o','LineStyle',Style_Matches.LineStyle,'Color','k', 'LineWidth', 0.5, 'MarkerFaceColor', Style_Matches.Color, 'MarkerSize', 8);
end
