%   File: 'prepare_Images.m'
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

function Images = prepare_Images(Before, After, Polygons)
    
    Images.crop_query = [];
    Images.polygon_query_gt = [];
    gray_Img1 = Adapt_Image_Channels(Before,1);
    gray_Img2 = Adapt_Image_Channels(After,1);
    Images.Training = gray_Img1;
    Images.Query = gray_Img2;
    Images.polygon_training = Polygons.Training;