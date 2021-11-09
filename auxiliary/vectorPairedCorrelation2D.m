function pairedCorrelation = vectorPairedCorrelation2D(u,v,x,y,normalise,center)
% computes the instantaneous cross correlation of vectors
% input:
% u,v: vectors eg velocity/displacement in two dimensions - N by 1
% normalise: flag whether to normalise vectors by magnitude
if nargin < 4
    center = false;
end
if normalise
    s = sqrt(u.^2 + v.^2);
    d = sqrt(x.^2 + y.^2);
    % normalise for vector magnitude
    u = u./s;
    v = v./s;
    x = x./d;
    y = y./d;
end
if center % NOTE: whether to normalise or center first depends on the application
    u = u - nanmean(u);
    v = v - nanmean(v);
    x = x - nanmean(x);
    y = y - nanmean(y);
end
pairedCorrelation = u.*x + v.*y;
