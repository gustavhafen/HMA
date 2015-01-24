%   File: 'Weighting_Function.m'
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

function W = Weighting_Function(residuals, tau)
% function that computes the weights (influence) for each dimension
% of the array constructed from the parameters (pos query x, pos query y,
% trans x, trans y, delta theta, delta sigma) given by the matches of the cluster. 
% The weights are computed from the quality of the current affine model,
% i.e., by measuring "how close" are the residual errors 'r' to the desired
% threshold 'tau'

%sigma = 10;
%r_bar = median(r);

% The vector W is formed by 3 values: 
%[weights position (query image), weights translation, weights 'similarity' (delta theta and delta sigma)]

% We observe that weight 'Wo' offers good performance in the initial
% clustering.
Wo = [0.25, 0.2, 0.05];
% When the number of outliers is reduced, i.e., the residual errors are
% small, we preffer to ensure the spatial contiguity of the clusters.
 Wf = [0.5, 0, 0];
% Wf = [0.4, 0.1 0];
% if we are in the first clustering, use the initial weights. 
if sum(isinf(residuals))>0,
    W = Wo;
else
    % The weights are estimated by interpolating between Wo and Wf:
    % W(t) = alpha*Wf+(1-alpha)*Wo.
    % if residual < 3tau => alpha = 1 => W=Wf, 
    % if residual ~ 3tau => alpha->1 => W->Wf,
    % if residual ~ 3sigma => alpha->0 => W->Wo.
    % We can convert this interpolation as non-linear by choosing alpha as a
    % non-linear continuous function with image in the interval [0,1].
    % We chose to model alpha as 1 if the residual errors are smaller than the threshold 'tau', 
    % otherwise, we use 2 times the half of a Gaussian cdf 
    % centered in the error threshold and with all the points 
    % (matches) inside the 3*sigma interval, i.e., 'sigma = max(residual)/3'
    
    sigma = max(residuals)/3;
    r_hat = median(residuals);
    
    if r_hat < 3*tau, % if the majority of the matches lie inside the threshold, then ensure spatial contiguity.
        alpha = 1;
    else
        %esimating the 1-area of the cdf of the Gaussian from the mean(error thresold) to the estimator (median or mean of the residuals).
        den = sqrt(2)*sigma;
        num = (r_hat -3*tau);
        alpha = 1 - erf(num/den);
    end
    
    % computing the interpolation
    W = Wo + alpha*(Wf-Wo);
end