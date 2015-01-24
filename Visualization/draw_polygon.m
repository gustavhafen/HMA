%   File: 'draw_polygon.m'
%
%   Author(s):  Gustavo A. Puerto-Souza
%   Created on: 2011
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

function draw_polygon(limits, style_params, color)
%function that draws a polygon "limits" in the current image with lines
%with stile defined at style_param
if nargin == 1, 
   style_params.Color = 'r';
   style_params.LineStyle = '-';
   style_params.Marker = '.';
   style_params.LW = 2;
   style_params.LW_ext = 3;
   style_params.Color_ext = 'k';
end
if isempty(style_params)
       style_params.Color = 'r';
       style_params.LineStyle = '-';
       style_params.Marker = '.';
       style_params.LW = 2;
       style_params.LW_ext = 3;
       style_params.Color_ext = 'k';
end
if nargin == 3, 
       style_params.Color = color;
end

line([limits(1,:), limits(1,:)],[limits(2,:), limits(2,:)],'Marker','none','LineStyle',style_params.LineStyle,'Color',style_params.Color_ext, 'LineWidth', style_params.LW_ext);
line([limits(1,:), limits(1,:)],[limits(2,:), limits(2,:)],'Marker',style_params.Marker,'LineStyle',style_params.LineStyle,'Color',style_params.Color, 'LineWidth', style_params.LW);


