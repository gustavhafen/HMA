%   File: 'compute_neighbohrs_linear_case.m'
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

function  neighbohrs = compute_neighbohrs_linear_case(real_value, bin_value, bin_size)
%function to compute if the value falls on the left or right side of the
%bin. The name linear case Is referred to the transformation to obtain the
%bin should be linear(for example for scale doesn't work).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real_value := real value
%bin_value  := value of the bin (index of the bin)
%bin_size   := value of the bin, or size of the bin delta_bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@return:
%neighbohrs = vector of the same size of real value, with the values of the
%             closest neighbohring bins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%does fall in left?  (retoric question)
left  = ( bin_value * bin_size - real_value > bin_size/2 );
%does fall in right? (retoric question)
right = ( bin_value * bin_size - real_value <= bin_size/2 );
%then assign
neighbohrs(left)  = bin_value(left)  - 1;
neighbohrs(right) = bin_value(right) + 1;