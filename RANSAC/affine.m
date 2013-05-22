function [ H ] = affine( ptsA, ptsB )
% This function retrieves an affine transform from given correspondences
% in the input, in the least squares sense
%
% Input elements
% -->PTSA: The initial set of homogeneous coordinates
% -->PTSB: The final set of homogeneous coordinates
%
% Output elements
% -->H: The estimated transform obtained from correspondences PTSA<-->PTSB
samples=size(ptsA,2);

% Normalizing input data points
[xA TA]=normalCoords(ptsA);
[xB TB]=normalCoords(ptsB);

% Initializing necessary arrays to solve the transform in the least squares
% sense
A=zeros(2*samples,6);
xp=zeros(2*samples,1);

if (samples>2)
    % Forming matrix A (A*h=xp)
    % h = [h11 h12 h13 h21 h22 h23];
    for k=1:samples
        A(2*k-1,:)=[xA(1,k) xA(2,k) 1 0 0 0];
        A(2*k,:)  =[0 0 0 xA(1,k) xA(2,k) 1];
        
        xp(2*k-1)=xB(1,k);
        xp(2*k)  =xB(2,k);
    end
    
    % Solving the system via SVD and denormalize  extracted transform
    [U D V]=svd(A);
    h=V*pinv(D)*U'*xp;
    H=[h(1) h(2) h(3);h(4) h(5) h(6);0 0 1];
    H=[reshape(h,3,2) [0 0 1]']';
    H=TB\H*TA;
else
   error('At least 3 samples are required to obtain an affine transform'); 
end

end