%   File: 'Build_Matrix_Affine.m'
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

function [A,b] = Build_Matrix_Affine(matched_keypoints)

n = size(matched_keypoints.Query,2);
Xq = matched_keypoints.Query(1:2,:)';
Xt = matched_keypoints.Training(1:2,:)';

R1 = [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 0 0 1 0];
  
R2 = [0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 0 1];

A = [ [Xt ones(n,1)] * R1;  [Xt ones(n,1)] * R2 ];
b = [Xq(:,1); Xq(:,2)];