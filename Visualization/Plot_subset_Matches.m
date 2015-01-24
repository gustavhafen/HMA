%   File: 'Plot_subset_Matches.m'
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

function Plot_subset_Matches(pts_t, ptsq, Images, indices, hfig, colors, symbol, linest, texto)

width = Plot_Dual_Image(Images.Training, Images.Query, hfig); hold on;
if nargin < 5,
    colors = 'r';
end
if nargin < 6,
    symbol = '.';
end
if nargin < 7,
    linest= 'none';
end
if nargin < 8,
    texto = 0;
end

if isempty(indices),
    indices = 1:size(pts_t,2);
end

line(  [pts_t(1, indices); width + ptsq(1, indices)], ...
    [pts_t(2, indices); ptsq(2, indices)],...
    'color', colors, 'marker', symbol, 'LineStyle', linest);
if texto == 1
    for i =1:length(indices),
        text(pts_t(1, indices(i)),pts_t(2, indices(i)), num2str(indices(i)));
        text(width+ptsq(1, indices(i)), ptsq(2, indices(i)), num2str(indices(i)));
    end
end