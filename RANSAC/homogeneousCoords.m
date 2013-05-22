function [ hx ] = homogeneousCoords( x )
% This function converts data points of different scales to a common
% scale of 1. Note that points with scale of 0 are left unchanged.
%
% Input elements
% -->x: input data points
%
% Output elements
% -->hx: output data points

% Check for valid data points array
if (size(x,1)~=3)
    error('Input structure must be 3xN');
end

% Initialize output
pts=size(x,2);
hx=zeros(size(x));

% Bring all data points to a common scale of 1
for i=1:pts
   if (x(3,i)==0)
       warning('zero scale point detected. Point left unchanged');
   else
       hx(:,i)=x(:,i)/x(3,i);
   end
end

end