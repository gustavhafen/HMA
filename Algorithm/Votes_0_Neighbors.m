%   File: 'Votes_0_Neighbors.m'
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

function [clusters, removed] = Votes_0_Neighbors(bins, min_votes)
% N number of votes
%function that computes the votes and discard clusters with minus than
%'min_votes' votes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bins := cell with all the unidimensional votes for that image(and only that image)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@return:
%clusters := cell with the arrays of indexes of each cluster. Each column 
%          of the cell contains a vector of indexes equivalent to one cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of matches
n = size(bins.X,2);
%% Generate all the votes(and the fifth column have the corresponding match
%to that vote)
    votes = [];
    votes(1:n, :) = [bins.pos_X(1,:)', bins.pos_Y(1,:)', bins.X(1,:)', bins.Y(1,:)', bins.Sigma(1,:)', bins.Theta(1,:)', [1:n]'];

%     for i=1:16
%         % Example:
%         % votes = [2 0 3 4 1;
%         %          .. .. .. ..]; 
%         % where 2 is the index of the bin for 'x', (0 for 'y', etc.) 
%         % and '1' is the index corresp. to which match was associated to that bin.
%         votes(n*(i-1)+1 : i*n, :) = [bins.X(ceil(i/8),:)', bins.Y(mod(ceil(i/4),2)+1,:)', bins.Sigma(mod(ceil(i/2),2)+1,:)',bins.Theta(mod(i,2)+1,:)' [1:n]'];
%     end
%% For each bin in hough generate the list of associated matches    
population = [];
clusters = [];
i=1; % scan over all the votes
while (size(votes,1) > 0) %
       candidate = votes(1,:); %vote to be contabilized
       rows = size(votes,1); 
       votes_to_remove  = logical(floor( sum((votes(:,1:6) == ones(rows,1) * candidate(1:6))')/6 )); %indexes of rows with same bins values
       clusters(i).ind = votes(votes_to_remove,7)';%Extract the index of matches associated with the candidate
       %counting how many votes have(size of cluster)
       population = [population sum(votes_to_remove)];
       %removing contabilized votes
       votes(logical(votes_to_remove),:) = [];
       i = i+1;
end
%removing clusters with less than 'min_votes' votes
removed = [clusters(population < min_votes).ind];
clusters(population < min_votes)  = [];