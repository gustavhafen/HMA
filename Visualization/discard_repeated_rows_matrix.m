%   File: 'discard_repeated_rows_matrix.m'
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

function [reduced_matrix reduced_indexes] = discard_repeated_rows_matrix(matrix)
%% function that returns the reduced(without repetead rows) matrix.
%INPUT:
% matrix : nxm matrix. 
%OUTPUT:
% reduced_matrix : n'xm matrix without repeated rows (according to the
%                  tolerance threshold = 1e-5).
    
    %% ALGORITHM
    
    reduced_matrix = [];
    i = 0;  % row counter
    n = size(matrix,1);
    indexes = 1:n;
    reduced_indexes = [];
    while ~isempty(matrix),
        i = i+1;
        row = matrix(1,:); %retrieving first row
        reduced_matrix(i,:) = row;
        reduced_indexes(i) = indexes(1);
        ind_m = find_vector_f(matrix, row); %looking for repeated elem.
        matrix(ind_m,:) = []; %removing repeated elements.
        indexes(ind_m) = [];
    end