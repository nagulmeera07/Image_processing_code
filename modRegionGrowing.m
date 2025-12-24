function [J] = modRegionGrowing(cIM, initPos, thresVal, maxLowDist, tfFillHoles, mask)
% REGIONGROWING Region growing algorithm for 2D grayscale images
%
% Syntax:
%   P = modRegionGrowing();
%   P = modRegionGrowing(cIM);
%   P = modRegionGrowing(cIM, initPos)
%   P = modRegionGrowing(..., thresVal, maxDist, tfFillHoles)
%   [P, J] = modRegionGrowing(...);
%
% Inputs:
%          cIM: 2D grayscale matrix                      {current image}
%      initPos: Coordinates for initial seed position     {ginput position}
%     thresVal: Absolute threshold level to be included     {5% of max-min}
%      maxDist: Maximum distance to the initial position in [px]      {Inf}
%  tfFillHoles: Fills enclosed holes in the binary mask              {true}
%
% Outputs:
%   J: Binary mask (with the same size as the input image) indicating
%      1 (true) for associated pixel/voxel and 0 (false) for outside
%
% Examples:
%   % 2D Example
%   load example
%   figure, imshow(cIM, [0 1500]), hold all
%   poly = modRegionGrowing(cIM, [], 300); % click somewhere inside the lungs
%
% Requirements:
%   TheMathWorks Image Processing Toolbox for imfill()
%   
% Remarks:
%   The queue is not preallocated. Seems to run faster this way.
%   A queue counter coud be tried 
%
% Author:
%   Heavily modified the version by Daniel Kellner. Traded robustness for
%   speed. Some functionality removed, others added. Accepts only 2D images now to 
%   run faster. Lines removed, and code vectorized to some extent to
%   to increase speed by by about 70%
%   Anand T.N.C anand@iitm.ac.in 2014/05/10
%   Copyright (c) 2014 Indian Institute of Technology Madras

% Original File: Copyright (c) 2011, Daniel released under BSD license
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% error checking on input arguments
if nargin > 6
    error('Wrong number of input arguments!')
end

if ~exist('cIM', 'var')
    himage = findobj('Type', 'image');
    if isempty(himage) || length(himage) > 1
        error('Please define one of the current images!')
    end
    
    cIM = get(himage, 'CData');
end

if ~exist('initPos', 'var') || isempty(initPos)
    himage = findobj('Type', 'image');
    if isempty(himage)
        himage = imshow(cIM, []);
    end
    
    % graphical user input for the initial position
    p = ginput(1);
    
    % get the pixel position concerning to the current axes coordinates
    initPos(1) = round(axes2pix(size(cIM, 2), get(himage, 'XData'), p(2)));
    initPos(2) = round(axes2pix(size(cIM, 1), get(himage, 'YData'), p(1)));
end

if ~exist('thresVal', 'var') || isempty(thresVal)
    thresVal = double((max(cIM(:)) - min(cIM(:)))) * 0.05;
end

if ~exist('maxDistSq', 'var') || isempty(maxDistSq)
    maxDistSq = Inf;
end

if ~exist('tfFillHoles', 'var')
    tfFillHoles = true;
end

if isequal(ndims(cIM), 2)
    
elseif isequal(ndims(cIM),1) || ndims(cIM) > 2
    error('Only 2D images are allowed!')
end

[nRow, nCol] = size(cIM);

if initPos(1) < 1 || initPos(2) < 1 ||...
        initPos(1) > nRow || initPos(2) > nCol
    error('Initial position is outside the image, please try again!');
end

if thresVal < 0 || maxDistSq < 0
    error('Threshold and maximum distance values must be positive!');
end

posX=initPos(1)-maxLowDist   ;
posY=initPos(2)-maxLowDist   ;

mapMask=false(nRow, nCol);
for i=1:nRow
    for j=1:nCol
        if(i-posX)> 0  &&  (i-posX) <=size(mask,1)
            if(j-posY)> 0  &&  (j-posY) <=size(mask,2)
                mapMask(i,j)=mask(i-posX,j-posY);
            end
        end
    end
end










% preallocate array
J = false(nRow, nCol);

% Initialising queue and adding the initial pixel to the queue
queue(1,:) = [initPos(1), initPos(2)];

%%% START OF REGION GROWING ALGORITHM
while any(queue(:))
    xv = queue(1,1);
    yv = queue(1,2);   % .. and remove the values in that row of queue
    queue(1,:) = [];
    % check the neighbors for the current position
    comb=[-1 -1; -1 0; -1 1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
    for i = 1:8
        if xv+comb(i,1) > 0  &&  xv+comb(i,1) <= nRow &&...          % within the x-bounds of image?
                yv+comb(i,2) > 0  &&  yv+comb(i,2)<= nCol &&...          % within the y-bounds of image?
                ~J(xv+comb(i,1), yv+comb(i,2)) &&...                     % pixel position already set?
                cIM(xv+comb(i,1), yv+comb(i,2)) >= thresVal &&... % greater than set value of the threshold?
                mapMask(xv+comb(i,1),yv+comb(i,2))  % within distance cutoff?
            % current pixel is true, if all properties are fullfilled
            J(xv+comb(i,1), yv+comb(i,2)) = true;
            
            queue(end+1,:) = [xv+comb(i,1), yv+comb(i,2)];
        end
    end
end
%%% END OF REGION GROWING ALGORITHM

% Fill holes unless the flag is set to 0 in the arguments passed
if tfFillHoles
    if any(J(:,:))
        % fill the holes inside the mask
        J(:,:) = imfill(J(:,:), 'holes');
        %B = bwboundaries(J(:,:,cSli), 8, 'noholes');
    end
end


% text output with final number of vertices
% disp(['RegionGrowing Ending: Found' num2str(size(P, 1)) ' polygon vertices)!'])
% if length(find(J))==0
%     disp(['RegionGrowing Ending: Found ' num2str(length(find(J)))...
%         ' pixels within the threshold range'])
end

