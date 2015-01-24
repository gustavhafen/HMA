%   File: 'HMA_Space_Quantization.m'
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

function [bins, parameters]= HMA_Space_Quantization(matched_keypoints, size_of_images, deltas)
%function auxiliar in Hough Transformation. Quantize the space of
%parameters (location[x,y], scale, orientation) and return the closest(and
%closest neighbohr) bin for each parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matched_keypoints := cell containing matched Keypoints.
%size_of_images   := matrix with the sizes of each image(in each col)
%criterion        := to determinate the size of bins to axis x and y.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@return:
%bins := a cell with the votes of the keypoint(of all images) %% extended
%(includes kpt query position)
%parameters :=  the size of bins for location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Hough, first part
%%%quantization%%%

%obtaining number of keypoints for each image
query_keypoints = matched_keypoints.Query;
training_keypoints = matched_keypoints.Training;
total_keypoints = matched_keypoints.num;
% number of bins

if nargin == 3,
    if deltas(1) == 0,
        delta_pos_x = 0.25*max(size_of_images(1));
    else
        delta_pos_x = deltas(1);
    end
    if deltas(2) == 0,
        delta_pos_y = 0.25*max(size_of_images(2));
    else
        delta_pos_y = deltas(2);
    end
    if deltas(3) == 0,
        delta_x = 0.25*max(size_of_images(1));
    else
        delta_x = deltas(3);
    end
    if deltas(4) == 0,
        delta_y = 0.25*max(size_of_images(2));
    else
        delta_y = deltas(4);
    end
    sigma_factor = deltas(5);
    delta_theta = deltas(6);
else
    delta_pos_x = 0.25*max(size_of_images(1));
    delta_pos_y = 0.25*max(size_of_images(1));
    delta_x = 0.25*max(size_of_images(1));
    delta_y = 0.25*max(size_of_images(2));
    sigma_factor = 2;
    delta_theta = (pi/6);
end
%% Hough obtaining projections, compute bins
%initializing
sigma_bin = zeros(2,total_keypoints);
theta_bin = zeros(2,total_keypoints);
x_bin = zeros(2,total_keypoints);
y_bin = zeros(2,total_keypoints);
%% new
pos_x_bin = zeros(2,total_keypoints);
pos_y_bin = zeros(2,total_keypoints);
%%
%the factor of scale, sigma = training_query(variable)/query_sigma
sigma = query_keypoints(3,:) ./ training_keypoints(3,:);
%values of the indexes bins
sigma_bin(1,:) = ceil(log2(sigma)); %closest point
%If the scale falls into the left side of the bin, the neighbohr is on the
%left, on the contrary falls on the right
left_neighbohrs = (sigma > 1.5 * sigma_factor.^(sigma_bin(1,:) - 1));
right_neighbohrs = logical(ones(1,total_keypoints) - left_neighbohrs);
%assigning values of the neighbohrs
sigma_bin(2,left_neighbohrs)  = sigma_bin(1,left_neighbohrs)  - 1;
sigma_bin(2,right_neighbohrs) = sigma_bin(1,right_neighbohrs) + 1;
%%Theta (Orientation)
%to reduce error, we apply the mod 2pi to all angles (because maybe the
%numbers are too different in value but very closer in angle, like 0.0001
%and 2pi*20, the values of orientation are not "normalized")
query_theta = mod(query_keypoints(4,:), 2*pi);
training_theta = mod(training_keypoints(4,:), 2*pi);
% a+b mod p = a mod p + b mod p
theta = mod(query_theta - training_theta, 2*pi);
theta_bin(1,:) = ceil(theta./delta_theta);
theta_bin(2,:) = compute_neighbohrs_linear_case(theta, theta_bin(1,:), delta_theta);
%boundary constrains: [12] right neighbohr  is [0],  and [1] left neighbohr
%is [12]
max_theta_bins = ceil(2*pi/delta_theta);

theta_bin(2, theta_bin(2,:) == 0) = max_theta_bins;
theta_bin(2, theta_bin(2,:) == max_theta_bins+1) = 1;

%% translation
%Xq = sigma*RotMat*Xt + T => T = Xq - sigma*RotMat*Xt
% tx = query_keypoints(1,:) - sigma .* (cos(theta) .* training_keypoints(1,:) - sin(theta) .* training_keypoints(2,:));
% ty = query_keypoints(2,:) - sigma .* (sin(theta) .* training_keypoints(1,:) + cos(theta) .* training_keypoints(2,:));

tx = query_keypoints(1,:) - training_keypoints(1,:);
ty = query_keypoints(2,:) - training_keypoints(2,:);

%closest bin
x_bin(1,:) = ceil(tx./delta_x);
y_bin(1,:) = ceil(ty./delta_y);
%neighbohrs
x_bin(2,:) = compute_neighbohrs_linear_case(tx, x_bin(1,:), delta_x);
y_bin(2,:) = compute_neighbohrs_linear_case(ty, y_bin(1,:), delta_y);

%% position
%closest bin
pos_x_bin(1,:) = ceil(query_keypoints(1,:)./delta_pos_x);
pos_y_bin(1,:) = ceil(query_keypoints(2,:)./delta_pos_y);
%neighbohrs
pos_x_bin(2,:) = compute_neighbohrs_linear_case(tx, pos_x_bin(1,:), delta_pos_x);
pos_y_bin(2,:) = compute_neighbohrs_linear_case(ty, pos_y_bin(1,:), delta_pos_y);

%packing bins
bins.pos_X = pos_x_bin;
bins.pos_Y = pos_y_bin;
bins.X = x_bin;
bins.Y = y_bin;
bins.Sigma = sigma_bin;
bins.Theta = theta_bin;
%parameters required later!
parameters.pos_X = delta_pos_x;
parameters.pos_Y = delta_pos_y;
parameters.X = delta_x;
parameters.Y = delta_y;
parameters.Sigma = sigma_factor;
parameters.Theta = delta_theta;
