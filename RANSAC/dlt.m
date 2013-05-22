function [ H A ] = dlt( ptsA, ptsB )
% This function computes a homography transformation from input point
% correspondences, using the dlt algorithm (Direct Linear Transform)
%
% Input elements
% ptsA: The initial set of points
% ptsB: The final set of points ( ptsA<--> ptsB )
%
% Output elements
% H: The extracted homography


% The structure of both point inputs must be of same size
if ~all(size(ptsA)==size(ptsB))
    error('Points structure must be of same size');
end

% Number of input points
points=size(ptsA,2);

% Normalizing the input coordinates
[ptsA TA]=normalCoords(ptsA);
[ptsB TB]=normalCoords(ptsB);

% Assembling A matrix which consists of all the necessery equations to
% obtain the final homography
A=zeros(3*points,9);
for i=1:points
    % Declaring variables for convenience
    xVec=ptsA(:,i)';
    x=ptsB(1,i);
    y=ptsB(2,i);
    z=ptsB(3,i);
    
    % Calculating Ai matrix and putting it to the final A matrix
    A(3*i-2,:)=[zeros(1,3) -z*xVec y*xVec];
    A(3*i-1,:)=[z*xVec zeros(1,3) -x*xVec];
    A(3*i,:)  =[-y*xVec x*xVec zeros(1,3)];
end

% Calculating the SVD decomposition of the A matrix and extracting the
% homography as the final column of the V matrix
[U D V]=svd(A,0);
H=reshape(V(:,9),3,3)';
H=TB\H*TA;

end