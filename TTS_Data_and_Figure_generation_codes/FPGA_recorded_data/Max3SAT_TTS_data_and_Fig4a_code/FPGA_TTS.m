%% plot_SAT_TTS_IQR_errorbars.m
% SAT TTS analysis with interquartile range as error bars.
% Assumes the current folder contains:
%   - coarse_summary.csv
%   - coarse_instance_tts.csv
%
% Prefix examples: 'uf100 ' => N = 100.

clear; clc;

%% Parameters
scale_factor = 1e-8;   % iterations -> seconds
delta_tol    = 1e-12;  % tolerance for matching thld_delta
quantiles_q  = [0.25, 0.50, 0.75];

% How large should NaN/Inf TTS be treated as, relative to the largest finite TTS?
LARGE_MULTIPLIER = 10;   % set to 100, 1000, etc. if you like
LARGE_FALLBACK   = 1e12; % used if NO finite TTS at all for that prefix

%% Load CSVs
summary = readtable('coarse_summary.csv');
inst    = readtable('coarse_instance_tts.csv');

% Ensure prefix columns are cell arrays of char
summary.prefix = normalize_prefix_column(summary.prefix);
inst.prefix    = normalize_prefix_column(inst.prefix);

%% Parse N from prefix like "uf100 " -> 100
summary.N = cellfun(@parseN_from_uf_prefix, summary.prefix);
inst.N    = cellfun(@parseN_from_uf_prefix, inst.prefix);

%% Unique prefixes and their sizes, sorted by N
[prefixList, ~, ic] = unique(summary.prefix, 'stable');
N_per_prefix = accumarray(ic, summary.N, [], @(x) x(1));

[sortedN, sortIdx] = sort(N_per_prefix);
sortedPrefixes     = prefixList(sortIdx);

nSizes = numel(sortedN);

%% Containers for quantiles (iterations)
TTS_q_iter = nan(nSizes, numel(quantiles_q));   % columns: [q25, q50, q75]

%% Main loop over sizes/prefixes
for k = 1:nSizes
    pref = sortedPrefixes{k};
    Nval = sortedN(k);

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
        this_q = quantile(TTS_all, q);
        TTS_q_iter(k, qi) = this_q;

        % If q = 0.25 or 0.75 and the quantile is effectively largeVal -> print message
        if (abs(q - 0.25) < 1e-9 || abs(q - 0.75) < 1e-9)
            if abs(this_q - largeVal) <= 1e-9 * max(1, largeVal)
                fprintf(['[INFO] For prefix %s (N=%d), q=%.2f TTS equals largeVal=%.3g ', ...
                         '(NaN/Inf/unsolved dominate that quantile).\n'], ...
                        pref, Nval, q, largeVal);
            end
        end
    end
end

%% Convert to seconds
TTS_q_sec = TTS_q_iter * scale_factor;
TTS_q25   = TTS_q_sec(:, 1);
TTS_q50   = TTS_q_sec(:, 2);   % median
TTS_q75   = TTS_q_sec(:, 3);

%% Prepare data for plotting (remove non-finite / non-positive)
valid = isfinite(TTS_q25) & isfinite(TTS_q50) & isfinite(TTS_q75) & ...
        (TTS_q25 > 0) & (TTS_q75 > 0);

N_vals = sortedN(valid);
Q25    = TTS_q25(valid);
Q50    = TTS_q50(valid);
Q75    = TTS_q75(valid);

if isempty(N_vals)
    error('No valid data to plot.');
end

% Asymmetric error bars for IQR
errLower = Q50 - Q25;   % distance from median down to q25
errUpper = Q75 - Q50;   % distance from median up to q75

%% Create paper-ready figure: median + IQR error bars
figure; clf; hold on;
set(gcf, 'Color', 'w');

% Aesthetic color
c = [0 0.4470 0.7410];  % MATLAB default blue

% Median with IQR error bars
errorbar(N_vals, Q50, errLower, errUpper, 'o-', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', c, ...
    'MarkerEdgeColor', c, ...
    'Color', c, ...
    'DisplayName', 'Median TTS (IQR)');

% After plotting errorbar(...)
yMin = min(Q25);
yMax = max(Q75);

padFactor = 1.5;    % 1.5× padding above/below the data
ylim([yMin / padFactor, yMax * padFactor]);
set(gca, 'YScale', 'log', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'TickLabelInterpreter', 'latex', ...
    'Layer', 'top');

%xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('Time-to-solution TTS (s)', 'Interpreter', 'latex', 'FontSize', 18);

box on;

xlim([min(N_vals)-10, max(N_vals)+10]);

% Set desired x-ticks
set(gca, 'XTick', [20 100 175 250]);

% No legend
legend('off');

fprintf('Done. Figure created with median TTS and IQR error bars.\n');


%% ---- Helper functions ---------------------------------------------
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
