%   File: 'Adapt_Image_Channels.m'
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

function Image = Adapt_Image_Channels(Src, Num_channels)
%function that converts the image to RGB or graymap

channels = size(Src,3);
if Num_channels == 1, %graymap
    if channels == 3, 
        Image = rgb2gray(Src);
    else
        Image = Src;
    end
else
    if channels == 3, 
        Image = Src;
    else
        Image(:,:,1) = Src;
        Image(:,:,2) = Src;
        Image(:,:,3) = Src;
    end
end