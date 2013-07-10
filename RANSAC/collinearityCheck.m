function [ check ] = collinearityCheck( xA, xB )
% This function checks whether a set of points of up to 4 samples, contains
% 3 points that are collinear. It is used by the ransac.m to determine if
% an MSS that we picked is degenerated.
%
% Input elements
% -->xA: The initial set of points
% -->xB: The final set of points (The algorithm checks xA and xB
% seperately)
%
% Output elements
% -->check: A logical variable showing if at least one collinearity was
% detected

% Initially, the collinearity check is false
check=0;

% Check for valid structure in arrays xA and xB
if (~all(size(xA)==size(xB)))
    error('Input points must be of same size');
end

% Number of samples
samples=size(xA,2);

% Form the matrices A and B such that their determinants point out if there
% is any collinearity in arrays xA and xB
if (samples==4)
    K=4;
    
    A=zeros(4,3,3);
    B=zeros(4,3,3);

    A(1,:,:)=[xA(:,1)';xA(:,2)';xA(:,3)'];
    A(2,:,:)=[xA(:,1)';xA(:,2)';xA(:,4)'];
    A(3,:,:)=[xA(:,1)';xA(:,3)';xA(:,4)'];
    A(4,:,:)=[xA(:,2)';xA(:,3)';xA(:,4)'];

    B(1,:,:)=[xB(:,1)';xB(:,2)';xB(:,3)'];
    B(2,:,:)=[xB(:,1)';xB(:,2)';xB(:,4)'];
    B(3,:,:)=[xB(:,1)';xB(:,3)';xB(:,4)'];
    B(4,:,:)=[xB(:,2)';xB(:,3)';xB(:,4)'];
elseif(samples==3)
    K=1;
    
    A=zeros(1,3,3);
    B=zeros(1,3,3);
    
    A(1,:,:)=[xA(:,1)';xA(:,2)';xA(:,3)'];
    B(1,:,:)=[xB(:,1)';xB(:,2)';xB(:,3)'];
elseif(samples==2);
    K=0;
else
    error('Number of samples must be 2, 3 or 4');
end

% Compute the determinants of matrices A and B and warn if any collinearity
% shows up
for i=1:K
    if (~det(squeeze(A(i,:,:)))) | (~det(squeeze(B(i,:,:))))
        check=1;
        warning('Collinearity detected');
        break;
    end
end

