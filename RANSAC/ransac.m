function [ Inliers Model ] = ransac( ptsA, ptsB, P, T, maxIterations, maxInvalidCount, errorThr )
% This function implements the ransac algorithm adapted to estimate a
% geometric transformation and exports the data that produce the estimated
% transformation and the transformation itself.
%
% Input elements
% -->PTSA: A 3xN array that holds the initial set of homogeneous
% coordinates
% -->PTSB: A 3xN array that holds the final set of homogeneous coordinates
% (in correspondence with PTSA)
% -->P: The desired probability to pick an MSS(Minimum Sample Set) with no 
% outliers
% -->T: The type of the exported transform. It must be one of the following
% strings: {'homography','affine','similarity','euclidean','translation'}
% -->MAXITERATIONS: The maximum number of iterations allowed
% -->MAXINVALIDCOUNT: The maximum number of serial degenerated set picks
% -->ERRORTHR: The error threshold that determines whether a data point is
% an inlier or an outlier
%
% Output elements
% -->Inliers: The data points that produce the exported transform
% -->Model: The estimated transform



% structure of arrays ptsA and ptsB
n1=size(ptsA);
n2=size(ptsB);

% Check transform validity
if ~(ischar(T))
   error('Input T must be a valid transform string'); 
else
    % Convert transform string to lower case for easier handling
    T=lower(T);
    
    % Is the transform string valid ?
    if ~(strcmp(T,'homography')...
            ||strcmp(T,'affine')...
            ||strcmp(T,'similarity')...
            ||strcmp(T,'euclidean')...
            ||strcmp(T,'translation'))
        
        error('Input T must be a valid string. Type "help ransac" for details');
    end     
end

% Check for structure validity of arrays ptsA and ptsB
if ~all(n1==n2)
    error('Structures ptsA and ptsB must be of same dimensions and sizes');
end

% Check for scale component in correspondences
if (n1~=3)
    error('Structures ptsA and ptsB must be 3xN');
end

% Check validity of probability P
if ~isnumeric(P)
    error('Input P is a probability and must be a number');
elseif ~(P>0 && P<1)
    error('Input P is a probability and must be in the interval (0,1)');
end


% Normalise the input points such that their centroid is the origin (0,0)
% and their mean distance from the centroid is sqrt(2)
[ptsA T1]=normalCoords(ptsA);
[ptsB T2]=normalCoords(ptsB);

% Number of given samples as input in the algorithm
samples=n1(2);

% Minimum number of correspondences needed to fully obtain a model. Depends
% on the input transform (e.g. 'homography', 'affine', 'similarity', 
% 'euclidean' or 'translation')
if (strcmp(T,'homography'))
    minS=4;
elseif (strcmp(T,'affine'))
    minS=3;
elseif (strcmp(T,'euclidean')||strcmp(T,'similarity'))
    minS=2;
elseif (strcmp(T,'translation'))
    minS=1;
   % BLANK---------------------- 
end

% Required samples check
if (samples<minS)
    error(sprintf('At least %d samples are required',minS));
end

% Function return variables: 
% 1) The estimated model
% 2) The inliers to this model
Model=[];
Inliers=[];


% Auxiliary variables
bestConsensus=0;
iterations=1;

% Iteration counter
i=0;

% A flag that shows that the algorithm exceeded the maximum number of
% iterations
failFlag=0;

while (i<iterations)
    % Check for maximum iterations violation
    if (i>maxIterations)
        warning('Maximum number of iterations reached. Exiting algorithm...');
        failFlag=1;
        break
    end
    
    % Degenerated configurations in a row count initialization
    iCount=0;
    breakFlag=1;
    
    while (iCount<maxInvalidCount && breakFlag)
        % Generate an MSS (Minimum Sample Set) - Indexes
        perm=randperm(samples);
        initSet=perm(1:minS);
        
        % MSS correspondences
        xA=ptsA(:,initSet);
        xB=ptsB(:,initSet);
        
        % Check for degeneracy
        breakFlag=collinearityCheck(xA, xB);
        
        iCount=iCount+1;
    end
    
    % Compute model from MSS correspondences
    
    if (strcmp(T,'homography'))
        tmpH=dlt(xA, xB);
    elseif (strcmp(T,'affine'))
        tmpH=affine(xA, xB);
    elseif (strcmp(T,'euclidean')||strcmp(T,'similarity'))
        tmpH=similarity(xA, xB);
    elseif (strcmp(T,'translation'))
        % ------------------ BLANK ---------------------
    end
    
    % Temporary consensus set. The inliers are marked with 1's
    tmpConSet=zeros(1,samples);
    
    % Measuring distances of datapoints from estimated model
    for j=minS+1:samples
        
        % j-th pair of correspondences
        initX  = homogeneousCoords(ptsA(:,perm(j)));
        finalX = homogeneousCoords(ptsB(:,perm(j)));
        
        
        % Scale normilization before error measurement
%         Hx=hnormalise(tmpH*initX);
%         iHxp=hnormalise(tmpH\finalX);
%         
%         initX=hnormalise(initX);
%         finalX=hnormalise(finalX);
        Hx=homogeneousCoords(tmpH*initX);
        iHxp=homogeneousCoords(tmpH\finalX);
        
        % Calculate error for current point (Symmetric Transfer Error)
        dF=sum((finalX-Hx).^2);
        dB=sum((initX-iHxp).^2);
        d=dF+dB;

        % Determine whether data point is an inlier or an outlier
        if (d<errorThr)
            tmpConSet(perm(j))=1;
        end
        
    end
    
    % Determine whether the current model is better than the previous one.
    % This is done by checking the cardinality of the new consensus set
    tmpConsensus=nnz(tmpConSet);
    if (tmpConsensus>bestConsensus)
        
        % In case of a better estimation, update the consensus of the best
        % set, keep track of the new inliers and the new candidate model
        bestConsensus=tmpConsensus;
        Inliers=find(tmpConSet==1);
        Model=tmpH;
        
        % BEWARE OF DIVISION WITH -INF AND 0 !!!!!!
        % SOME ESSENTIAL CHECKING IS NEEDED.
        % Estimate the new number of necessary iterations
        iterations=ceil(log(1-P)/log(1-(tmpConsensus/samples)^minS));
    end
    
    % Iteration counter update
    i=i+1;
    
end

if (length(Inliers)>minS-1)
    % Estimate final model using only the found inliers and denormalize model
    if (strcmp(T,'homography'))
        % Homography Case
        Model=T2\(dlt(ptsA(:,Inliers),ptsB(:,Inliers)))*T1;
    
    elseif(strcmp(T,'affine'))
        % Affine Case
        Model=T2\(affine(ptsA(:,Inliers),ptsB(:,Inliers)))*T1;

    elseif(strcmp(T,'similarity'))
        % Similarity Case
        Model=T2\(similarity(ptsA(:,Inliers),ptsB(:,Inliers)))*T1;
       
    elseif(strcmp(T,'euclidean'))
       % Euclidean Case (similar to Similarity Case)
       Model=T2\(similarity(ptsA(:,Inliers),ptsB(:,Inliers)))*T1;
       
       % Eliminate scaling in extracted transform
       [U S V]=svd(Model(1:2,1:2));
       Model(1:2,1:2)=Model(1:2,1:2)*pinv(S);
       
    end
else
    if ~failFlag
        fprintf(1,'Not enough inliers were foun.Try changing configuration of the algorithm\n');
    end
end

% Executed iterations
fprintf(1,'\nRANSAC: Done in %d Iterations\n',i);

end