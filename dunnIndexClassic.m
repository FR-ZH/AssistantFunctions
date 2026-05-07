function dunn = dunnIndexClassic(X, labels)
% Classic Dunn index (larger is better).
% Inter-cluster distance: minimum pairwise distance between clusters.
% Intra-cluster diameter: maximum pairwise distance within clusters.

u = unique(labels(:));
k = numel(u);

% intra-cluster diameters
diam = zeros(k,1);
for ii = 1:k
    Xi = X(labels == u(ii), :);
    if size(Xi,1) > 1
        diam(ii) = max(pdist(Xi, 'euclidean'));
    else
        diam(ii) = 0;
    end
end
maxDiam = max(diam);

if maxDiam == 0
    dunn = NaN;
    return;
end

% inter-cluster min distance
minInter = Inf;
for i = 1:k-1
    Xi = X(labels == u(i), :);
    for j = i+1:k
        Xj = X(labels == u(j), :);
        if isempty(Xi) || isempty(Xj), continue; end
        minInter = min(minInter, min(pdist2(Xi, Xj, 'euclidean'), [], 'all'));
    end
end

dunn = minInter / maxDiam;
end