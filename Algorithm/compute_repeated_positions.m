%   File: 'compute_repeated_positions.m'
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

function ind_repetitions = compute_repeated_positions(Initial_Matches, tol)

if nargin < 2,
    tol = 0.001;
end
matrix_matches = [Initial_Matches.Training(1,:); Initial_Matches.Training(2,:); Initial_Matches.Query(1,:); Initial_Matches.Query(2,:)]';
num_matches = Initial_Matches.num;
indices = 1:num_matches;
ind_repetitions = zeros(1, num_matches);
class = 0;
while ~isempty(indices);
    class = class + 1;
    candidate = matrix_matches(1,:);
    ind_rep = ( (abs(candidate(1)-matrix_matches(:,1))<tol) & (abs(candidate(2)-matrix_matches(:,2))<tol) ) | ...
              ( (abs(candidate(3)-matrix_matches(:,3))<tol) & (abs(candidate(4)-matrix_matches(:,4))<tol) ); % trainig or query positions repeated        
    ind_repetitions(indices(ind_rep)) = class;
    indices(ind_rep) = [];
    matrix_matches(ind_rep,:) = [];
end