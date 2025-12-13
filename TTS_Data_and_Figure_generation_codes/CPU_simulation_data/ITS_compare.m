%% plot_SAT_TTS_three_algs.m
% Plot SAT TTS scaling (median + IQR error bars) for three algorithms:
%   1) HO_color_full_precision
%   2) HO_single_spin_flip
%   3) QUBO
%
% Folder structure (relative to this script's location):
%   HO_color_full_precision/coarse_summary.csv
%   HO_color_full_precision/coarse_instance_tts.csv
%   HO_single_spin_flip/coarse_summary.csv
%   HO_single_spin_flip/coarse_instance_tts.csv
%   QUBO/coarse_summary.csv
%   QUBO/coarse_instance_tts.csv

clear; clc;

%% Global parameters
scale_factor    = 1;   % iterations -> "iteration-to-solution" units
delta_tol       = 1e-12;  % tolerance for matching thld_delta
quantiles_q     = [0.25, 0.50, 0.75];

% How large should NaN/Inf TTS be treated as, relative to the largest finite TTS?
LARGE_MULTIPLIER = 10;    % set to 100, 1000, etc. if you like
LARGE_FALLBACK   = 1e12;  % used if NO finite TTS at all for that prefix

baseDir = pwd;

%% Analyze each algorithm
dummy1 = repmat(char(160), 1, 36);   % 15 "invisible" space
res_color = analyze_one_sat_alg( ...
    fullfile(baseDir, 'HO_color_full_precision'), ...
    dummy1, ...
    scale_factor, delta_tol, quantiles_q, ...
    LARGE_MULTIPLIER, LARGE_FALLBACK);

res_spin  = analyze_one_sat_alg( ...
    fullfile(baseDir, 'HO_single_spin_flip'), ...
    dummy1, ...
    scale_factor, delta_tol, quantiles_q, ...
    LARGE_MULTIPLIER, LARGE_FALLBACK);

res_qubo  = analyze_one_sat_alg( ...
    fullfile(baseDir, 'QUBO'), ...
    dummy1, ...
    scale_factor, delta_tol, quantiles_q, ...
    LARGE_MULTIPLIER, LARGE_FALLBACK);

%% Combined plot: median TTS vs N with IQR error bars for all three
figure; clf; hold on;
set(gcf, 'Color', 'w');

% Colors & markers for the three algorithms
colors  = [0 0.4470 0.7410;    % blue
           0.8500 0.3250 0.0980; % orange
           0.4660 0.6740 0.1880];% green
markers = {'o', 's', '^'};

allRes   = {res_color, res_spin, res_qubo};
nAlgs    = numel(allRes);
legendHandles = gobjects(nAlgs, 1);
legendLabels  = cell(nAlgs, 1);

for t = 1:nAlgs
    res = allRes{t};

    N_vals = res.N_vals;
    Q25    = res.Q25;
    Q50    = res.Q50;
    Q75    = res.Q75;

    if isempty(N_vals)
        continue;
    end

    errLower = Q50 - Q25;
    errUpper = Q75 - Q50;

    c  = colors(t, :);
    mk = markers{t};

    % Median with IQR error bars
    h = errorbar(N_vals, Q50, errLower, errUpper, [mk '-'], ...
        'LineWidth', 1.8, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', c, ...
        'MarkerEdgeColor', c, ...
        'Color', c);

    legendHandles(t) = h;
    legendLabels{t}  = res.algName;
end

set(gca, 'YScale', 'log', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'TickLabelInterpreter', 'latex', ...
    'Layer', 'top');

%xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('Iteration to solution ($\times 10^{-8}$)', ...
%       'Interpreter', 'latex', 'FontSize', 18);

box on;

% Overall x-limits from all algorithms
allN = [res_color.N_vals; res_spin.N_vals; res_qubo.N_vals];
allN = allN(~isnan(allN));
if ~isempty(allN)
    xlim([min(allN)-10, max(allN)+10]);
end

hLeg = legend(legendHandles, legendLabels, ...
       'Location', 'northwest', ...
       'Interpreter', 'none');

pos = get(hLeg, 'Position');

% Example: shrink width and height by 20%, keep same center
pos(4) = pos(4) * 1.4;   % height
set(hLeg, 'Position', pos);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% Optional:
% title('SAT iteration-to-solution scaling for three algorithms', ...
%       'Interpreter', 'latex', 'FontSize', 18);

fprintf('Done. Figure created with three algorithms (median + IQR error bars).\n');

%% ========================================================================
function res = analyze_one_sat_alg(folder, algName, ...
                                   scale_factor, delta_tol, quantiles_q, ...
                                   LARGE_MULTIPLIER, LARGE_FALLBACK)
    % Analyze SAT TTS for a single algorithm in the given folder.
    % Returns median & IQR vs N after scaling by scale_factor.

    % Load CSVs
    summary = readtable(fullfile(folder, 'coarse_summary.csv'));
    inst    = readtable(fullfile(folder, 'coarse_instance_tts.csv'));

    % Ensure prefix columns are cellstr
    summary.prefix = normalize_prefix_column(summary.prefix);
    inst.prefix    = normalize_prefix_column(inst.prefix);

    % Parse N from prefix like "uf100 " -> 100
    summary.N = cellfun(@parseN_from_uf_prefix, summary.prefix);
    inst.N    = cellfun(@parseN_from_uf_prefix,  inst.prefix);

    % Unique prefixes and their sizes, sorted by N
    [prefixList, ~, ic] = unique(summary.prefix, 'stable');
    N_per_prefix = accumarray(ic, summary.N, [], @(x) x(1));

    [sortedN, sortIdx] = sort(N_per_prefix);
    sortedPrefixes     = prefixList(sortIdx);

    nSizes = numel(sortedN);

    TTS_q_iter = nan(nSizes, numel(quantiles_q));   % rows: size index, cols: q25,q50,q75

    for k = 1:nSizes
        pref = sortedPrefixes{k};
        Nval = sortedN(k); %#ok<NASGU>

        % ----- pick thld_delta that minimizes median_TTS_val for this prefix -----
        maskSum = strcmp(summary.prefix, pref);
        subSum  = summary(maskSum, :);

        if isempty(subSum)
            warning('No summary rows found for prefix %s.', pref);
            continue;
        end

        [~, idxMin] = min(subSum.median_TTS_val);
        delta_star  = subSum.thld_delta(idxMin);

        % ----- gather instance-wise data at (prefix, delta_star) -----
        maskInst = strcmp(inst.prefix, pref) & ...
                   abs(inst.thld_delta - delta_star) < delta_tol;
        subInst  = inst(maskInst, :);

        if isempty(subInst)
            warning('No instance rows found for prefix %s with thld_delta=%.6g.', ...
                    pref, delta_star);
            continue;
        end

        TTS_i_iter = subInst.TTS_at_tau_star;   % iterations

        % ====== treat NaN/Inf as large numbers for quantiles ======
        finiteMask = isfinite(TTS_i_iter);
        if any(finiteMask)
            maxFin   = max(TTS_i_iter(finiteMask));
            largeVal = LARGE_MULTIPLIER * maxFin;
        else
            largeVal = LARGE_FALLBACK;
        end

        TTS_all = TTS_i_iter;
        TTS_all(~finiteMask) = largeVal;

        % ----- quantiles q = 0.25, 0.5, 0.75 (over ALL instances) -----
        for qi = 1:numel(quantiles_q)
            q = quantiles_q(qi);
            TTS_q_iter(k, qi) = quantile(TTS_all, q);
        end
    end

    % Convert to scaled units
    TTS_q_sec = TTS_q_iter * scale_factor;
    Q25_all   = TTS_q_sec(:, 1);
    Q50_all   = TTS_q_sec(:, 2);
    Q75_all   = TTS_q_sec(:, 3);

    % Valid data (positive and finite IQR bounds)
    valid = isfinite(Q25_all) & isfinite(Q50_all) & isfinite(Q75_all) & ...
            (Q25_all > 0) & (Q75_all > 0);

    N_vals = sortedN(valid);
    Q25    = Q25_all(valid);
    Q50    = Q50_all(valid);
    Q75    = Q75_all(valid);

    % Pack results
    res.algName = algName;
    res.N_vals  = N_vals;
    res.Q25     = Q25;
    res.Q50     = Q50;
    res.Q75     = Q75;
end

%% ------------------------------------------------------------------------
function prefOut = normalize_prefix_column(prefIn)
    % Normalize a prefix column from readtable into a cell array of char.
    if isstring(prefIn)
        prefOut = cellstr(prefIn);
    elseif iscell(prefIn)
        if ~isempty(prefIn) && ischar(prefIn{1})
            prefOut = prefIn;
        else
            prefOut = cellstr(string(prefIn));
        end
    elseif ischar(prefIn)
        prefOut = cellstr(prefIn);
    else
        prefOut = cellstr(string(prefIn));
    end
end

function N = parseN_from_uf_prefix(pref)
    % Extract N from prefixes like "uf100", "uf100 " or "uf100-01", etc.
    s = strtrim(pref);
    tok = regexp(s, 'uf(\d+)', 'tokens', 'once');
    if isempty(tok)
        error('Could not parse N from prefix "%s". Expected pattern "uf<digits>".', pref);
    end
    N = str2double(tok{1});
end
