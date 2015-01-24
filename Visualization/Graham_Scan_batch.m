%   File: 'Graham_Scan_batch.m'
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

function [convex_hull indexes] = Graham_Scan_batch(points)
%Algorith to find the convex hull of a cloud of points
%% inputs:
% points: 2xn matrix with the (x,y) coordinates of the points.
%% outputs:
% convex_hull: ordered points belonging to the convex hull in clockwise order)
% indexes: index of the respective point in the convex hull in the cloud of
%          points

%% special case
if size(points,2)<2,
    convex_hull = points;
    indexes = 1;
    return 
end

%% selecting anchor point p0 (lowest point)
min_y = min(points(2,:)); 
min_idxs = find(points(2,:) == min_y); %In the general case the lowest point can be multiple.
if length(min_idxs) > 1
     min_x= min(points(1, min_idxs));
     p0_idx = min_idxs(points(1,min_idxs) == min_x);
     p0 = points(:, p0_idx); %hyphotesis that no are 2 equal points (relaxing the non-colinearity constrain)
else
    p0_idx = min_idxs;
    p0 = points(:, p0_idx);    
end

%% constructing non-convex polygon
points(3,:) = 1:size(points,2); %indexes
points(:,p0_idx) = []; %discarding lowest point
n = size(points,2);
dif = points(1:2,:) - p0*ones(1,n);
cosines = (points(1,:) - p0(1))./sum(sqrt(dif.*dif));
[AA,cos_idx] = sort(cosines, 'descend');
points = [[p0; p0_idx]  points(:,cos_idx)];


%% Constructing convex polygon
M = 2;
i=3;
x_y = []; y_y = [];
while i<n+2,
    %checking the 'direction' of the turn (left or right)
    while ccw(points(:,M-1), points(:,M), points(:,i)) <= 0,
          %left turn
          if M == 2,
                  aux = points(:,M);
                  points(:,M) = points(:,i);
                  points(:,i) = aux;
                  i = i + 1;
          else
                  M = M - 1;
          end
    end
    %right turn
    M = M + 1;
    aux = points(:,M);
    points(:,M) = points(:,i);
    points(:,i) = aux;
    i=i+1;
end

%% Solution
convex_hull = points(1:2,1:M);
indexes = points(3,1:M);