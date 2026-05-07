function [idxA, CA, map] = align_labels_to_centers(idx, C, Cref)
%ALIGN_LABELS_TO_CENTERS Align k-means labels to reference centroids (within same feature space).
%
% This function resolves label switching across repeated k-means runs by
% mapping each run's cluster indices to a reference labeling defined by
% reference centroids Cref (k×p). It assumes C (k×p) and Cref (k×p) are in the
% same feature space (same number of columns p).
%
% Inputs
%   idx   - n×1 cluster labels from one k-means run (values 1..k)
%   C     - k×p centroids from that run
%   Cref  - k×p reference centroids (fixed ordering)
%
% Outputs
%   idxA  - n×1 aligned labels, mapped into the reference label space (1..k)
%   CA    - k×p centroids reordered to match reference label order (row j corresponds to ref label j)
%   map   - 1×k mapping vector: map(i)=j means run label i is mapped to reference label j
%
% Notes
%   - Uses matchpairs (Hungarian assignment) if available; otherwise uses a greedy fallback.
%   - This function does NOT compare across different feature sets with different p.

k = size(Cref, 1);
D = pdist2(C, Cref);  % k×k distances between run centers and reference centers

% Optimal assignment if matchpairs exists
if exist('matchpairs', 'file') == 2
    pairs = matchpairs(D, 1e9);   % [i j] = run cluster i -> ref cluster j
    map = zeros(1,k);
    for t = 1:size(pairs,1)
        map(pairs(t,1)) = pairs(t,2);
    end
else
    % Greedy fallback
    map = zeros(1,k);
    used = false(1,k);
    for i = 1:k
        [~, ord] = sort(D(i,:), 'ascend');
        for j = ord
            if ~used(j)
                map(i) = j;
                used(j) = true;
                break;
            end
        end
    end
end

% Remap labels
idxA = map(idx);

% Reorder centroids to reference order
CA = nan(size(C));
for i = 1:k
    j = map(i);
    CA(j,:) = C(i,:);
end
end