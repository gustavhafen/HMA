%   File: 'Set_Param_Features.m'
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

function Params_Features = Set_Param_Features(Feat_Type)
%Setting parameters for the feature extraction.

if strcmp(Feat_Type,'SIFT'),
    Params_Features.Type = 0;
    Params_Features.peak = 0.5;
    Params_Features.edge = 800;
else if strcmp(Feat_Type,'SURF'),
        Params_Features.Type = 1;
        Params_Features.upright = true;
        Params_Features.thresh = 0.0001;
     else
         Params_Features.Type = 2;
         fpath = fileparts(which(mfilename));
         Params_Features.exec_str = ['"' fpath '\demo_ASIFT"'];
    end
end