%   File: 'RUN_SCRIPT.m'
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

%% Demo with SIFT features
%Required Libraries

%VL_Feat, available at:
%http://www.vlfeat.org/
%ADD library to the path, and then run this script.

clear
clc
close all
run('vl_setup'); % run once
cpath = pwd;
addpath(genpath(cpath));
% addpath('.\Visualization','.\Algorithm');

% Demo with the two images (Images), extracted SIFT features (Features), 
% and an intial matching (Initial_Matches), as  well as the region of 
% interest in the training image (Polygons).
load Demo;

%% parameters:
% basic parameters (RANSAC and matching parameters)
Params = Set_HMA_parameters('SIFT');

% Stop condition (% oof inliers) at each node level: 
ratio_threshold = 0.9;
% model (affine only)
model = 'A';

Solution = HMA_Algorithm(Features, Params, Images, Initial_Matches, ratio_threshold, model);

% Visualization
Display_Params = Set_Display_Params();
Plot_Multi_Model(Solution, Images, Display_Params, []);