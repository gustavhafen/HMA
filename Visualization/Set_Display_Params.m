%   File: 'Set_Display_Params.m'
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

function Params_Display = Set_Display_Params()

%% General options
Params_Display.Image_Channels = 1; %image in gray (1) or color (3), if is available.
Params_Display.Style = 1; % 0 - no plots, 1 in one window, 2 in separate windows.
% Plot Style
%Matches
Matches.LW = 1; %LW=line width
Matches.Color = 'g'; %Color for matches and line betw. matches
Matches.Marker = 's';%'square'
Matches.LineStyle = 'none'; % between training and query
Params_Display.Matches = Matches;
%Reprojection error
Rep.LineStyle = '-';
Rep.Color = 'y';
Rep.Marker = 'o';
Rep.LW = 1;
Params_Display.Rep = Rep;
% Trainign image: Polygons Style
Pol_Training.Marker = '.';
Pol_Training.Color = 'r';
Pol_Training.LW = 2;
Pol_Training.LineStyle = '-';
Pol_Training.LW_ext = 3;
Pol_Training.Color_ext = 'k';
Params_Display.Pol_Training = Pol_Training;
%
Pol_Rep.Marker = '.';
Pol_Rep.Color = 'r';
Pol_Rep.LW = 2;
Pol_Rep.LineStyle = '-';
Pol_Rep.LW_ext = 3;
Pol_Rep.Color_ext = 'k';
Params_Display.Pol_Rep = Pol_Rep;
%
Pol_GT.Marker = '.';
Pol_GT.Color = 'c';
Pol_GT.LW = 2;
Pol_GT.LineStyle = ':';
Pol_GT.LW_ext = 3;
Pol_GT.Color_ext = 'k';
Params_Display.Pol_GT = Pol_GT;