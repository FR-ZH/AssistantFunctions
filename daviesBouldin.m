function db = daviesBouldin(X, labels)
% DAVIESBOULDIN Compute Davies–Bouldin index (lower is better).
%   db = daviesBouldin(X, labels)
%   X: n×p observations, labels: n×1 cluster labels (can be non-consecutive)

u = unique(labels(:));   % actual label values
k = numel(u);

centers = zeros(k, size(X,2));
S       = zeros(k,1);

for ii = 1:k
    lab = u(ii);
    Xi  = X(labels == lab, :);

    if isempty(Xi)
        centers(ii,:) = NaN;
        S(ii) = NaN;
        continue;
    end

    centers(ii,:) = mean(Xi, 1);
    di = sqrt(sum((Xi - centers(ii,:)).^2, 2));  % Euclidean
    S(ii) = mean(di);
end

M = pdist2(centers, centers, 'euclidean');

% avoid division by zero on diagonal
M(1:k+1:end) = Inf;

R = zeros(k,k);
for i = 1:k
    for j = i+1:k
        if ~isfinite(S(i)) || ~isfinite(S(j)) || ~isfinite(M(i,j)) || M(i,j)==0
            Rij = NaN;
        else
            Rij = (S(i) + S(j)) / M(i,j);
        end
        R(i,j) = Rij;
        R(j,i) = Rij;
    end
end

db = mean(max(R, [], 2, 'omitnan'), 'omitnan');

end