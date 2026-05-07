function CH = calinskiHarabaszIndex(X, idx)
% Calinski–Harabasz index:
% CH = (B/(k-1)) / (W/(n-k))
% B: between-cluster dispersion, W: within-cluster dispersion

[n, ~] = size(X);
labels = unique(idx(:));
k = numel(labels);

% overall mean
mu = mean(X,1);

W = 0;  % within
B = 0;  % between

for i = 1:k
    lab = labels(i);
    Xi = X(idx==lab, :);
    ni = size(Xi,1);
    if ni == 0, continue; end

    mui = mean(Xi,1);

    % within dispersion: sum ||x - mui||^2
    difW = Xi - mui;
    W = W + sum(difW(:).^2);

    % between dispersion: ni * ||mui - mu||^2
    difB = mui - mu;
    B = B + ni * sum(difB.^2);
end

CH = (B/(k-1)) / (W/(n-k));
end