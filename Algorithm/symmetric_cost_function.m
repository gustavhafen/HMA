%   File: 'symmetric_cost_function.m'
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
function cost = symmetric_cost_function(X)

%% reading components of the function
% Transformation
a = X(1:6);  
M = reshape(a,2,3);
M = [M; 0 0 1];
A = M(1:2,1:2);
T = M(1:2,3);
M_inv = inv(M);  
A_inv = M_inv(1:2,1:2);
T_inv = M_inv(1:2,3);
% Keypoints
remaining = X(7:end);
num_points = length(remaining)/4; %2 sets each one of 2 dimensions
query_keypoints = reshape(remaining(1:2*num_points), 2, num_points);
db_keypoints = reshape(remaining(2*num_points+1:end), 2,num_points);

%% projecting db points
proj_db_keypoints = A*db_keypoints + T*ones(1, num_points);
proj_q_keypoints = A_inv*query_keypoints + T_inv*ones(1, num_points);
%% computing error
diff_in_query = query_keypoints - proj_db_keypoints;
diff_in_db = db_keypoints - proj_q_keypoints;
cost = sqrt(sum(diff_in_query.*diff_in_query) + sum(diff_in_db.*diff_in_db));