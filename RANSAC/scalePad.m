function [ x ] = scalePad( x, scaleVec )
% This function pads a points array (2xN) with a third row that contains
% the scales of the points. By default, the padded values are 1's
%
% Input elements
% -->x: Input data points array
% -->scaleVec: Contains the desired scales and if given, it is padded itself
% to the data points array
%
% Output elements
% -->x: The scale padded data points

% Check data points array for valid structure(2xN)
if size(x,1)~=2
    error('In order to pad a points array with scales, it must be 2xN');
end

% Number of samples;
N=size(x,2);

% Pad data points array with 1's or with the scale vector if it's provided
if nargin<2
    % Pad the scales to the data points array ( 1's )
    x=[x;ones(1,N)];
else
    % Check validity of the scales vector
    if size(scaleVec,2)~=N
        error('Data points array and scale vector must have same number of samples');
    elseif ~isrow(scaleVec)
        error('Scale vector must be a row vector');
    end
    
    % Pad the scale vector to the data points array
    x=[x;scaleVec];
end

end