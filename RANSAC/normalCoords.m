function [ nPts T ] = normalCoords( pts )
% This function implements a normalization to the input coordinates such
% that their scale is at 1, their centroid is at the origin (0,0) and their 
% mean distance from that origin is sqrt(2). It is an essential 
% preprocessing step before any estimation of a transform due to limited 
% arithmetic accuracy
%
% Input elements
% -->PTS: The set of points to be normalized 
%
% Output elements
% -->NPTS: The normalized set of points
% -->T: The implemented normalization as a transform matrix. It is useful
% for denormalization

% Scale-normalised coordinates
pts(1,:)=pts(1,:)./pts(3,:);
pts(2,:)=pts(2,:)./pts(3,:);
pts(3,:)=1;

% Calculate the centroid of the coordinates and subtract it to shift them
% to the origin (0,0)
center=mean(pts(1:2,:)')';
tmpPts=pts;
tmpPts(1,:)=tmpPts(1,:)-center(1);
tmpPts(2,:)=tmpPts(2,:)-center(2);

% Calculate the mean distance from their new centroid
distances=sqrt(tmpPts(1,:).^2+tmpPts(2,:).^2);
mDist=mean(distances(:));

% The new scaling factor
s=sqrt(2)/mDist;

% The normalization process as a transformation matrix
T=[s 0 -s*center(1);0 s -s*center(2);0 0 1];
nPts=T*pts;

end

