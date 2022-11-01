function STA = simpleSTA(Stim,sp,nkt,maxsize,varargin)
% [STA] = simpleSTC(Stim,sp,nkt,maxsize)
%
% Computes mean and covariance of spike-triggered (or response weighted)
% stimuli and raw stimuli
%
% INPUT:
%    Stim [N x M]   - stimulus matrix; 1st dimension = time, 2nd dimension = space
%      sp [N x 1]   - column vector of spike count in each time bin (can be sparse)
%         [Nsp x 1] - or, vector of spike times
%     nkt [1 x 1]   - # time samples to consider to be part of the stimulus
% MaxSize [1 x 1] (optional)   - max # of floats to store while computing cov
%                                (smaller = slower, but less memory required)
%
% OUTPUT:
%    STA [nkt x M]       - spike-triggered average (reshaped as matrix)
%
%  Notes:  
%  (1) Ignores response before nkt bins
%  (2) Faster if only 2 output arguments (raw mean and covariance not computed)
%  (3) Reduce 'maxsize' if getting "out of memory" errors
%  (4) If inputs are not spike times or counts, response-weighted covariance may be more appropriate.
%      See  "simpleRWC.m". 
%
%  --------
%  Details:
%  --------
%   Let X = "valid" design matrix (from makeStimRows)
%       Y = sp (spike train)
%     nsp = sum(Y); % number of spikes (sum of responses)
%       N = size(Stim,1); % number of stimuli
%
%   then  STA = X'*Y / nsp
%         RawMu = sum(X)/N;
%         RawCov = X*X'/(N-1);
%
% Copyright 2010 Pillow Lab. All rights reserved.
% % Modified by Suva Roy, 2021. 

%-------- Parse inputs  ---------------------------
if nargin < 4
    maxsize = 1e9; % max chunk size; decrease if getting "out of memory"
end
[slen,swid] = size(Stim); % stimulus size (# time bins x # spatial bins).


p = inputParser;
addParameter(p,'weighted',false, @islogical);
parse(p,varargin{:});
params = p.Results;

% Convert list of spike times to spike-count vector, if necessary
nsp = length(sp);
if (nsp ~= slen)
    fprintf(1, 'simpleSTC: converting spike times to counts\n');
    sp = hist(sp,1:slen)';
end
sp(1:nkt-1) = 0;  % Remove spikes before time n
spsum = sum(sp);  % Sum of spikes
if (spsum == 0) 
    spsum = 2;
end
nspnds = sum(sp~=0);


% ---------------------------------------------------
% 1. Compute only the spike-triggered STA

Msz = nspnds*swid*nkt; % size of design matrix for spiking stimuli
iisp = find(sp); 
splen = length(iisp);

if Msz < maxsize  % Compute in one chunk if small enough
    SS = makeStimRows(Stim, nkt,iisp);
    STA = (SS'*sp(iisp))/spsum;
    
else % Compute in multiple chunks if too large
    nchunk = ceil(Msz/maxsize);
    chunksize = ceil(length(iisp)/nchunk);
    fprintf(1, 'simpleSTC: using %d chunks to compute STA/STC\n', nchunk);
    
    % Initialize on 1st chunk
    i0 = 1;
    imx = chunksize;
    SS = makeStimRows(Stim,nkt,iisp(i0:imx));
    STA = (sp(iisp(i0:imx))'*SS)';
    
    % Compute for remaining chunks
    for j = 2:nchunk
        i0 = chunksize*(j-1)+1;
        imx = min(chunksize*j, splen);
        SS = makeStimRows(Stim,nkt,iisp(i0:imx));
        STA = STA + (sp(iisp(i0:imx))'*SS)';
    end
    % normalize by number of samples
    STA = STA/spsum;

end
    
STA = reshape(STA,[],swid);  % reshape to have same width as stimulus
