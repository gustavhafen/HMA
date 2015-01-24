%   File: 'HMA_Hough_Transform.m'
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

function Hough_Clustering = HMA_Hough_Transform(initial_matching, size_of_images, min_votes, deltas, Images)

%% Adaptation to -1 vote -

% Hough Transformation: For all the initial_matching (in all the images)
% finds indexes of keypoints with a least 3 different keypoints (clusters) (in DB) and no repeated(clusters).
%also we can add bool to add plots in each stage.

[bins, param_hough] = HMA_Space_Quantization(initial_matching, size_of_images,deltas);
bplot = 0;
%% Making votes
% clusters = Votes_16_Neighbors(bins, min_votes); %creating clusters
[clusters, removed1] = Votes_0_Neighbors(bins, min_votes); %creating clusters
%if ~isempty(clusters),
% no more repeated clusters
%     clusters = Discard_Repeated_Clusters(clusters, min_votes); %remove repeated clusters
%[clusters, removed2] = Discard_Ill_Conditioned_Clusters(clusters, initial_matching, min_votes);%remove clusters with less than 3 different points
%end
%
if bplot,
    if ~isempty([clusters.ind])
        Plot_Initial_Matches_subset2(initial_matching,[clusters.ind], Images, 'y');
    end
    if ~isempty(removed1)
        Plot_Initial_Matches_subset2(initial_matching,removed1, Images, 'r');
    end
    %         Plot_Initial_Matches_subset2(initial_matching,removed2, Images, 'b');
end
Hough_Clustering.Clusters = [clusters.ind];
Hough_Clustering.Params = param_hough;
Hough_Clustering.removed = removed1; %[removed1, removed2]