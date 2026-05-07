clc;

B = 1000;                  % bootstrap iterations (outer loop)
nRuns = 100;               % kmeans reruns per bootstrap sample (inner loop for voting)
ks = [3 4 5];              % sensitivity across k
rng(1);                    % reproducible at the whole-script level

% ---------- Put ALL your feature versions here ----------
Xsets = { ...
    Data_Ready_Z, ...
    Data_Ready_Z_Spectra, ...
    Data_Ready_Z_PCMI, ...
    Data_Ready_Z_PE, ...
};

Xnames = { ...
    'AllFeatures', ...
    'Spectral_only', ...
    'PCMI_only', ...
    'PE_Only', ...
};

% ---------- Sanity checks ----------
assert(numel(Xsets) == numel(Xnames), 'Xsets and Xnames must have same length.');
nSets = numel(Xsets);
n = size(Xsets{1},1);
for s = 1:nSets
    assert(size(Xsets{s},1)==n, 'All Xsets must have the same number of rows (subjects).');
end

% ---------- Storage ----------
Results = struct();
for s = 1:nSets
    Results(s).name = Xnames{s};
    Results(s).n = size(Xsets{s},1);
    Results(s).p = size(Xsets{s},2);
    Results(s).B = B;
    Results(s).nRuns = nRuns;
    Results(s).ks = ks;

    for kk = 1:numel(ks)
        Results(s).k(kk).k = ks(kk);

        Results(s).k(kk).sil_mean = nan(B,1);
        Results(s).k(kk).dbi      = nan(B,1);
        Results(s).k(kk).dunn     = nan(B,1);
        Results(s).k(kk).stab_med = nan(B,1); % median voteFrac (stability)

        % optional: keep voteFrac distribution if you want
        % Results(s).k(kk).voteFrac_all = cell(B,1);
    end
end

%% ---------- Bootstrap outer loop ----------
tStart = tic;

for b = 1:B

    % ---- progress display ----
    if b == 1
        fprintf('Starting bootstrap (%d iterations)...\n', B);
    end

    if mod(b,10)==0 || b==1
        elapsed = toc(tStart);
        pct = 100 * b / B;
        estTotal = elapsed / b * B;
        estRemain = estTotal - elapsed;

        fprintf('Bootstrap %4d / %4d (%.1f%%) | Elapsed: %.1f min | Remaining: %.1f min\n', ...
            b, B, pct, elapsed/60, estRemain/60);
    end

    % ---- paired bootstrap ----
    rng(b);
    bootIdx = randsample(n, n, true);

    for s = 1:nSets
        X = Xsets{s};
        Xb = X(bootIdx, :);

        for kk = 1:numel(ks)
            k = ks(kk);

            out = kmeans_majorityVotingAligned(Xb, k, nRuns, ...
                'Start','plus', 'Distance','sqeuclidean', ...
                'RefRun','bestSumd');

            Results(s).k(kk).sil_mean(b) = mean(silhouette(Xb, out.idx_vote));
            Results(s).k(kk).dbi(b)      = daviesBouldin(Xb, out.idx_vote);
            Results(s).k(kk).dunn(b)     = dunnIndexClassic(Xb, out.idx_vote);
            Results(s).k(kk).stab_med(b) = median(out.voteFrac);
        end
    end
end

fprintf('Bootstrap finished in %.2f minutes.\n', toc(tStart)/60);

%% ---------- Summaries: median + 95% bootstrap CI ----------
fprintf('\n=== Summary (median [2.5,97.5]) ===\n');
for s = 1:nSets
    fprintf('\nDataset: %s (n=%d, p=%d)\n', Results(s).name, Results(s).n, Results(s).p);
    for kk = 1:numel(ks)
        k = ks(kk);

        sil = Results(s).k(kk).sil_mean;
        dbi = Results(s).k(kk).dbi;
        dun = Results(s).k(kk).dunn;
        stb = Results(s).k(kk).stab_med;

        fprintf('  k=%d | Sil=%.4f [%.4f, %.4f] | DBI=%.4f [%.4f, %.4f] | Dunn=%.4f [%.4f, %.4f] | Stab=%.3f [%.3f, %.3f]\n', ...
            k, median(sil), prctile(sil,2.5), prctile(sil,97.5), ...
               median(dbi), prctile(dbi,2.5), prctile(dbi,97.5), ...
               median(dun), prctile(dun,2.5), prctile(dun,97.5), ...
               median(stb), prctile(stb,2.5), prctile(stb,97.5));
    end
end
%% ---------- Paired tests with direction (k=4 vs others) ----------
fprintf('\n=== Paired signrank tests (bootstrap-paired) + direction + 95%% CI ===\n');
targetK = 4;
s = 1;
    fprintf('\nDataset: %s\n', Results(s).name);

    for kk = 1:numel(ks)
        if kk == idxT, continue; end
        k_other = ks(kk);

        fprintf('  k=%d vs k=%d:\n', targetK, k_other);

        % ========= DBI (lower is better) =========
        vT = Results(s).k(idxT).dbi;
        vO = Results(s).k(kk).dbi;
        d  = vT - vO;  % negative => target better
        p  = signrank(vT, vO);

        fprintf('    DBI  (lower better): Δ=%.4f [%.4f, %.4f], p=%.3g  | favors k=%d if Δ<0\n', ...
            median(d,'omitnan'), prctile(d,2.5), prctile(d,97.5), p, targetK);

        % ========= Silhouette (higher is better) =========
        vT = Results(s).k(idxT).sil_mean;
        vO = Results(s).k(kk).sil_mean;
        d  = vT - vO;  % positive => target better
        p  = signrank(vT, vO);

        fprintf('    Sil  (higher better): Δ=%.4f [%.4f, %.4f], p=%.3g  | favors k=%d if Δ>0\n', ...
            median(d,'omitnan'), prctile(d,2.5), prctile(d,97.5), p, targetK);

        % ========= Dunn (higher is better) =========
        vT = Results(s).k(idxT).dunn;
        vO = Results(s).k(kk).dunn;
        d  = vT - vO;  % positive => target better
        p  = signrank(vT, vO);

        fprintf('    Dunn (higher better): Δ=%.4f [%.4f, %.4f], p=%.3g  | favors k=%d if Δ>0\n', ...
            median(d,'omitnan'), prctile(d,2.5), prctile(d,97.5), p, targetK);

        % ========= Stability (higher is better) =========
        if isfield(Results(s).k(idxT), 'stab_med') && isfield(Results(s).k(kk), 'stab_med')
            vT = Results(s).k(idxT).stab_med;
            vO = Results(s).k(kk).stab_med;
            d  = vT - vO;  % positive => target more stable
            p  = signrank(vT, vO);

            fprintf('    Stab (higher better): Δ=%.4f [%.4f, %.4f], p=%.3g  | favors k=%d if Δ>0\n', ...
                median(d,'omitnan'), prctile(d,2.5), prctile(d,97.5), p, targetK);
        end
    end

%% ---------- Boxplots: compare k=3/4/5 together ----------
s=1;
    % ===== Silhouette =====
    figure('Name',['Dataset: ' Results(s).name ' | Silhouette']);
    silMat = nan(B, numel(ks));
    for kk = 1:numel(ks)
        silMat(:,kk) = Results(s).k(kk).sil_mean;
    end
    boxplot(silMat, 'Labels', string(ks));
    xlabel('k'); ylabel('Mean silhouette');
    title(['Mean silhouette across bootstrap (higher is better) | ' Results(s).name]);

    % ===== DBI =====
    figure('Name',['Dataset: ' Results(s).name ' | DBI']);
    dbiMat = nan(B, numel(ks));
    for kk = 1:numel(ks)
        dbiMat(:,kk) = Results(s).k(kk).dbi;
    end
    boxplot(dbiMat, 'Labels', string(ks));
    xlabel('k'); ylabel('Davies–Bouldin index ');
    title(['DBI across bootstrap (lower is better)| ' Results(s).name]);

    % ===== Dunn =====
    figure('Name',['Dataset: ' Results(s).name ' | Dunn']);
    dunMat = nan(B, numel(ks));
    for kk = 1:numel(ks)
        dunMat(:,kk) = Results(s).k(kk).dunn;
    end
    boxplot(dunMat, 'Labels', string(ks));
    xlabel('k'); ylabel('Dunn index ');
    title(['Dunn across bootstrap (higher is better)| ' Results(s).name]);

    % ===== Stability =====
    figure('Name',['Dataset: ' Results(s).name ' | Stability']);
    stabMat = nan(B, numel(ks));
    for kk = 1:numel(ks)
        stabMat(:,kk) = Results(s).k(kk).stab_med;
    end
    boxplot(stabMat, 'Labels', string(ks));
    xlabel('k');
    ylabel('Median subject-level stability (higher is better)');
    title(['Stability across bootstrap (higher is better)| ' Results(s).name]);
%% ==========================================================
% k=4 | Compare each feature set vs baseline (s=1)
% Use paired bootstrap difference Δ and its 95% CI (no need for all pairwise)
% Baseline = sBase (default 1)
% ==========================================================

kTarget = 4;
kkT = find(ks==kTarget,1);
assert(~isempty(kkT), 'kTarget=%d not found in ks.', kTarget);

sBase = 1;   % baseline feature set index
baseName = Xnames{sBase};

fprintf('\n==================================================\n');
fprintf('k=%d | Feature sets compared vs baseline: %s (s=%d)\n', kTarget, baseName, sBase);
fprintf('Δ defined as (Set - Baseline)\n');
fprintf('  Sil/Dunn/Stab: Δ>0 => better than baseline\n');
fprintf('  DBI          : Δ<0 => better than baseline (lower is better)\n');
fprintf('==================================================\n');

% Pull baseline vectors (B×1)
sil0  = Results(sBase).k(kkT).sil_mean(:);
dbi0  = Results(sBase).k(kkT).dbi(:);
dun0  = Results(sBase).k(kkT).dunn(:);
stab0 = Results(sBase).k(kkT).stab_med(:);

for s = 1:nSets
    if s == sBase, continue; end

    name = Xnames{s};

    sil  = Results(s).k(kkT).sil_mean(:);
    dbi  = Results(s).k(kkT).dbi(:);
    dun  = Results(s).k(kkT).dunn(:);
    stab = Results(s).k(kkT).stab_med(:);

    % Paired bootstrap differences Δ = (set - baseline)
    d_sil  = sil  - sil0;   % >0 better
    d_dbi  = dbi  - dbi0;   % <0 better
    d_dun  = dun  - dun0;   % >0 better
    d_stab = stab - stab0;  % >0 better

    fprintf('\n%s vs %s:\n', name, baseName);

    fprintf('  Sil  (higher better): Δ=%.4f [%.4f, %.4f]  | better if Δ>0\n', ...
        median(d_sil,'omitnan'), prctile(d_sil,2.5), prctile(d_sil,97.5));

    fprintf('  DBI  (lower  better): Δ=%.4f [%.4f, %.4f]  | better if Δ<0\n', ...
        median(d_dbi,'omitnan'), prctile(d_dbi,2.5), prctile(d_dbi,97.5));

    fprintf('  Dunn (higher better): Δ=%.4f [%.4f, %.4f]  | better if Δ>0\n', ...
        median(d_dun,'omitnan'), prctile(d_dun,2.5), prctile(d_dun,97.5));

    fprintf('  Stab (higher better): Δ=%.4f [%.4f, %.4f]  | better if Δ>0\n', ...
        median(d_stab,'omitnan'), prctile(d_stab,2.5), prctile(d_stab,97.5));

end