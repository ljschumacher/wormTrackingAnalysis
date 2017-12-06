function crossCorrelation = vectorCrossCorrelation2D(u,v,normalise,center)
% computes the instantaneous cross correlation of vectors
% input:
% u,v: vectors eg velocity/displacement in two dimensions - N by 1
% normalise: flag whether to normalise vectors by magnitude
if nargin < 4
    center = false;
end
if normalise
    s = sqrt(u.^2 + v.^2);
    % normalise for vector magnitude
    u = u./s;
    v = v./s;
end
if center
    u = u - mean(u);
    v = v - mean(v);
end
crossCorrelation = squareform(tril(... % only keep each pair once as cross-corrrelation should be symmetric
    [u v]*[u v]'... % calculate the scalar product between displacement vectors
    ,-1)); % ignores auto correlation
