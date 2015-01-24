%   File: 'Plot_Dual_Image.m'
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

function xl = Plot_Dual_Image(Im_l, Im_r, fig_handle) 
%function plot two images side by side. Assuming same number of channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@parameters:
%Im_l    := left image
%Im_r    := right image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return:
%x1 :=  the with of Im_l %could be extended to the dimensioms of teh images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% dimensions of the images
if isempty(Im_l), %for now  with non-sense
    xl = 0; return; %for now  with non-sense    
else
    [yl,xl,zl] = size(Im_l);    
end

if isempty(Im_r),
    yr = 0; xr = 0; zr = 0;
    return;
else
	[yr,xr,AA] = size(Im_r);
end

I = 255*ones(max(yl,yr),xl+xr,zl);

I(1:yl,1:xl,:) = Im_l;
if ~isempty(Im_r),
    I(1:yr,1+xl:xl+xr,:) = Im_r;
end
%% plotting
figure(fig_handle);
imshow(uint8(I));