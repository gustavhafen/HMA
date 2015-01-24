%   File: 'Tree_Clustering.m'
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

function  [clusters,K] = Tree_Clustering(Initial_Matches, W, polar, min_pts, indices)

coeff = 1;
offset = 2;
if nargin < 5,
    indices = 1:Initial_Matches.num;
end
% selecting the subset of data to partition
Initial_Matches = struct('Training', Initial_Matches.Training(:,indices),...
                         'Query', Initial_Matches.Query(:,indices), ...
                         'Original_indices', Initial_Matches.Original_indices(indices), ...
                         'num', length(indices));
%% Rotation
%t = Initial_Matches.Query(1:2, :) - Initial_Matches.Training(1:2, :);
sigma = Initial_Matches.Query(3,:) ./ Initial_Matches.Training(3,:);
query_theta = mod(Initial_Matches.Query(4,:), 2*pi);
training_theta = mod(Initial_Matches.Training(4,:), 2*pi);

% a+b mod p = a mod p + b mod p
theta = mod(query_theta - training_theta, 2*pi);
%% translation
translation = zeros(2,Initial_Matches.num);
for i=1:Initial_Matches.num,
    Rotation = sigma(i)*[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
    translation(:,i) = Initial_Matches.Query(1:2, i) - Rotation*Initial_Matches.Training(1:2, i);    
end

if polar == 1,
    translation_ro = sqrt(translation(1,:).^2 + translation(2,:).^2);
    translation_theta = mod(atan2(translation(2,:),translation(1,:)),2*pi);
    translation = [translation_ro; translation_theta];
end

%% packing
data_points = [Initial_Matches.Query(1:2, :); translation; sigma; theta;];
%% normalizing
for i=5:6,
    data_points(i,:) = data_points(i,:)/norm(data_points(i,:));
end
%% Weighting
wp=W(1);
wt=W(2);
ws=W(3);

data_points(1,:) = wp*data_points(1,:);
data_points(2,:) = wp*data_points(2,:);
data_points(3,:) = wt*data_points(3,:);
data_points(4,:) = wt*data_points(4,:);
data_points(5,:) = ws*data_points(5,:);
data_points(6,:) = ws*data_points(6,:);

% Adaptive selection of the clusters
  if Initial_Matches.num > 4*(coeff*min_pts+offset),
     K = 4;
  else
     K = max(2,floor(Initial_Matches.num/(coeff*min_pts+offset))); %% check the factor 2!!
  end
[mute, clusters] = vl_kmeans(data_points, K, 'initialization', 'plusplus', 'algorithm', 'elkan', 'numrepetitions', 3);          
