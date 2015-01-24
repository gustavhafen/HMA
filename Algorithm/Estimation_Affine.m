%   File: 'Estimation_Affine.m'
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

function M = Estimation_Affine(A, b)

%     if nargin < 4,
tol = 1e-6;%10;
%     end
AA = A'*A;
if rank(AA, tol) == 6, %ill conditioned
    Affine = AA\(A'*b); %(A'A)^{-1} * (A'b)
    R = reshape(Affine(1:4),2,2)';
    T = Affine(5:6);
    M.Transformation = [[R T]; 0 0 1];
    Affine_inv = R\[eye(2), -T]; % M = [A]b] -> M^{-1} = [A^{-1}| -A^{-1} b]
    M.Transformation_inv = [Affine_inv; 0 0 1];
    M.num_p = 6;
else
    M =[];
end