function [ mapping matches ] = matchDescriptors( descA, descB, ratio )
% This function implements a descriptor matching of the descriptors of two
% images according to their correlation coefficients.
%
% Input elements:
% descA: The descriptors of the first image (each column is handled as a
% descriptor)
%
% descB: The descriptors of the second image (each column is handled as a
% descriptor)
%
% ratio: The maximum ratio that is acceptable between the two best
% correlation coefficients for every attempted matching
%
% Output elements:
% mapping: A vector that holds the mapping between the two given sets of
% descriptors. The indexes of the vector correspond to the indexes of descA
% and the numbers correspond to the indexes of descB. Unsuccesful matchings
% are marked as 0's

% Number of descriptors for each image
NoDescA=size(descA,2);
NoDescB=size(descB,2);

% Vector that holds the mapping between the descriptors
mapping=zeros(1,NoDescA);

for i=1:NoDescA
    % Calculate dot product of each descriptor of the first image with
    % every descriptor in the second image and find their correlation
    % coefficients
    d=dot(descA(:,i),descB(:,i))/(norm(descA(:,i))*norm(descB(:,i)));
    [sd indexes]=sort(acos(d));
    
    % If the nearest descriptor is not relatively close to the next one,
    % then confirm a match
    if (sd(1)/sd(2) < ratio)
        mapping(i)=indexes(1);
    else
        mapping(i)=0;
    end
end

% Number of successful matches
matches=nnz(mapping);

end