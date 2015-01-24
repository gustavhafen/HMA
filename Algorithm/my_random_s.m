%   File: 'my_random_s.m'
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

function sample = my_random_s(n_array, k_array, max_samples)
%n_array size of all the set
%k_array size of the subsets
%max_samples maximal number of combinations

%% internal parameters. For now fixed
MAX_ITERATIONS = 100000; 
MAX_SET = 10;
%%
if n_array<MAX_SET, 
    C = nchoosek(1:n_array,k_array);
    n = size(C,1);
    idxs = randperm(n);
    if ~isempty(max_samples),
        num_samples = min(max_samples, n);
        sample = C(idxs(1:num_samples),:)';
    else
        sample = C';
    end
else 
    if ~isempty(max_samples),
        sample = random('unid',n_array,[k_array max_samples]);
    else
        sample = random('unid',n_array,[k_array MAX_ITERATIONS]);
    end
end