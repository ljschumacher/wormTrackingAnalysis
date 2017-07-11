function crossCorrelation = vectorCrossCorrelation2D(u,v,normalise)
% computes the instantaneous cross correlation of velocity/displacment
% vectors
% input:
% u,v: vectors of velocity/displacement in two dimensions - N by 1
% normalise: flag whether to normalise displacments by magnitude

if normalise
    s = sqrt(u.^2 + v.^2);
    % normalise for velocity magnitude
    u = u./s;
    v = v./s;
end
crossCorrelation = squareform(tril(... % only keep each pair once as cross-corrrelation should be symmetric
    [u v]*[u v]'... % calculate the scalar product between displacement vectors
    ,-1)); % ignores auto correlation
