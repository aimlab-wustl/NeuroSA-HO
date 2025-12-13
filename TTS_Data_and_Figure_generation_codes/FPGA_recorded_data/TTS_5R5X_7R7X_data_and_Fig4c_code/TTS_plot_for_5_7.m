%% analyze_XORSAT_5R5X_7R7X.m
% Analyze 5R5X and 7R7X XORSAT TTS scaling and plot them together.
% Run this script from the parent folder that contains:
%   - XORSAT_5R5X/coarse_summary.csv
%   - XORSAT_5R5X/coarse_instance_tts.csv
%   - XORSAT_7R7X/coarse_summary.csv
%   - XORSAT_7R7X/coarse_instance_tts.csv

clear; clc;

%% Global parameters
scale_factor   = 1e-8;   % iterations -> seconds
fitIdx_for_all = [];     % [] => use all sizes N for fitting
q_target       = 0.5;    % only q = 0.5

baseDir = pwd;

%% Analyze 5R5X and 7R7X separately
res5 = analyze_one_XORSAT_type(fullfile(baseDir, 'XORSAT_5R5X'), '5R5X', ...
                               scale_factor, fitIdx_for_all, q_target);

res7 = analyze_one_XORSAT_type(fullfile(baseDir, 'XORSAT_7R7X'), '7R7X', ...
                               scale_factor, fitIdx_for_all, q_target);

%% Combined plot: median TTS vs N with error bars and fits for both 5R5X and 7R7X
figure; clf; hold on;
set(gcf, 'Color', 'w');

colors  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];  % 5R5X, 7R7X
markers = {'o', 's'};

allRes = {res5, res7};

for t = 1:numel(allRes)
    res = allRes{t};

    Nvals   = res.N;
    medTTS  = res.median_sec;
    medErr  = res.err_sec;
    c       = colors(t, :);
    mk      = markers{t};
    typeStr = res.type;

    % --- Median points with error bars (marker+errorbar in ONE object) ---
    errorbar(Nvals, medTTS, medErr, ...
        mk, ...                       % marker style (o or s)
        'LineWidth', 1.8, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', c, ...
        'Color', c, ...
        'DisplayName', sprintf('%s median TTS', typeStr));

    % --- Fitted scaling line: log10(TTS_q) = alpha * N + beta ---
    alpha = res.alpha;
    beta  = res.beta;

    if ~isnan(alpha)
        N_line  = linspace(min(Nvals), max(Nvals), 200);
        TTS_fit = 10.^(alpha * N_line + beta);

        % Legend label with formatted alpha, beta (with parentheses)
        fitLabel = sprintf('%s fit: $q=0.5, \\alpha=%s, \\beta=%s$', ...
                           typeStr, res.alpha_fmt, res.beta_fmt);

        plot(N_line, TTS_fit, '-', ...
            'Color', c, ...
            'LineWidth', 1.8, ...
            'DisplayName', fitLabel);
    end
end

set(gca, 'YScale', 'log', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'TickLabelInterpreter', 'latex');

%xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('Optimal median TTS (s)', 'Interpreter', 'latex', 'FontSize', 18);

box on;

allN = [res5.N; res7.N];
xlim([min(allN)-2, max(allN)+2]);
ylim([5e-5, 1.7e-1]);

%legend('Location', 'northwest', 'Interpreter', 'latex');
legend off;

% Optional title
% title('5R5X and 7R7X optimal median TTS vs problem size', ...
%     'Interpreter', 'latex', 'FontSize', 18);

%% Save alpha/beta for q=0.5 to CSV in parent folder
typeCol  = {'5R5X'; '7R7X'};
qCol     = [q_target; q_target];

alphaCol = [res5.alpha;     res7.alpha];
alphaStd = [res5.alpha_std; res7.alpha_std];
alphaFmt = string({res5.alpha_fmt; res7.alpha_fmt});

betaCol  = [res5.beta;     res7.beta];
betaStd  = [res5.beta_std; res7.beta_std];
betaFmt  = string({res5.beta_fmt; res7.beta_fmt});

T_out = table(typeCol, qCol, ...
              alphaCol, alphaStd, alphaFmt, ...
              betaCol,  betaStd,  betaFmt, ...
    'VariableNames', {'type', 'q', ...
                      'alpha', 'alpha_std', 'alpha_fmt', ...
                      'beta',  'beta_std',  'beta_fmt'});

writetable(T_out, 'alpha_beta_TTS_5R5X_7R7X.csv');

fprintf('Done. Figure created and alpha_beta_TTS_5R5X_7R7X.csv written.\n');

%% ========================================================================
function res = analyze_one_XORSAT_type(folder, typeName, scale_factor, fitIdx_for_all, q_target)
    % Analyze a single XORSAT type (5R5X or 7R7X) in the given folder.

    delta_tol = 1e-12;  % tolerance for matching thld_delta

    % Load CSVs
    summary = readtable(fullfile(folder, 'coarse_summary.csv'));
    inst    = readtable(fullfile(folder, 'coarse_instance_tts.csv'));

    % Ensure prefix columns are cellstr
    summary.prefix = normalize_prefix_column(summary.prefix);
    inst.prefix    = normalize_prefix_column(inst.prefix);

    % Parse N from prefix like "..._n112" -> 112
    summary.N = cellfun(@parseN_from_prefix, summary.prefix);
    inst.N    = cellfun(@parseN_from_prefix,    inst.prefix);

    % Unique prefixes and their sizes, sorted by N
    [prefixList, ~, ic] = unique(summary.prefix, 'stable');
    N_per_prefix = accumarray(ic, summary.N, [], @(x) x(1));

    [sortedN, sortIdx] = sort(N_per_prefix);
    sortedPrefixes     = prefixList(sortIdx);

    nSizes = numel(sortedN);

    opt_thld_delta      = zeros(nSizes, 1);
    opt_median_TTS_iter = zeros(nSizes, 1);
    median_err_iter     = zeros(nSizes, 1);
    TTS_q_iter          = nan(nSizes, 1);  % q=0.5 across instances

    for k = 1:nSizes
        pref = sortedPrefixes{k};

        % ----- pick thld_delta that minimizes median_TTS_val for this prefix -----
        maskSum = strcmp(summary.prefix, pref);
        subSum  = summary(maskSum, :);

        if isempty(subSum)
            warning('No summary rows found for prefix %s.', pref);
            continue;
        end

        [minMedian, idxMin] = min(subSum.median_TTS_val);
        delta_star          = subSum.thld_delta(idxMin);
        opt_thld_delta(k)   = delta_star;
        opt_median_TTS_iter(k) = minMedian;  % iterations

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
        Ps         = subInst.Ps_hat;
        Nrep       = subInst.total_trials;

        % ====== treat NaN/Inf as large number for median quantile ======
        finiteMask0 = isfinite(TTS_i_iter);
        if any(finiteMask0)
            maxFin   = max(TTS_i_iter(finiteMask0));
            largeVal = 10 * maxFin;  % “large” to represent unsolved
        else
            largeVal = 1e12;         % fallback
        end

        TTS_all_for_q = TTS_i_iter;
        TTS_all_for_q(~finiteMask0) = largeVal;

        % quantile q (here q_target = 0.5) over ALL instances
        TTS_q_iter(k) = quantile(TTS_all_for_q, q_target);

        % ----- error bar on median TTS (SciAdv-style S10/S11) -----
        finiteMaskEB = isfinite(TTS_i_iter) & (Ps > 0) & (Ps < 1);
        TTS_finite   = TTS_i_iter(finiteMaskEB);
        Ps_finite    = Ps(finiteMaskEB);
        Nrep_finite  = Nrep(finiteMaskEB);

        if isempty(TTS_finite)
            warning('No finite TTS with 0 < Ps < 1 for prefix %s.', pref);
            continue;
        end

        % ΔP_S = sqrt(P_S - P_S^2)/sqrt(N_rep)
        dPs = sqrt(Ps_finite - Ps_finite.^2) ./ sqrt(Nrep_finite);

        % ΔTTS = TTS * ΔP_S / ((1 - P_S)*ln(1 - P_S))
        denom  = (1 - Ps_finite) .* log(1 - Ps_finite);
        dTTS_i = TTS_finite .* (dPs ./ denom);

        % Approximate SE of the median
        nF             = numel(TTS_finite);
        se_median_iter = sqrt(sum(dTTS_i.^2)) / nF;
        median_err_iter(k) = se_median_iter;
    end

    % Convert to seconds
    median_sec = opt_median_TTS_iter * scale_factor;
    err_sec    = median_err_iter     * scale_factor;
    TTS_q_sec  = TTS_q_iter          * scale_factor;

    % Scaling fit: log10(TTS_q) = alpha * N + beta
    valid = isfinite(TTS_q_sec) & (TTS_q_sec > 0);
    x = sortedN(valid);
    y = log10(TTS_q_sec(valid));

    alpha = NaN;
    beta  = NaN;
    alpha_std = NaN;
    beta_std  = NaN;

    if numel(x) >= 3
        if ~isempty(fitIdx_for_all)
            validIdx = fitIdx_for_all;
            validIdx = validIdx(validIdx <= numel(x));
            x_fit = x(validIdx);
            y_fit = y(validIdx);
        else
            x_fit = x;
            y_fit = y;
        end

        nFit = numel(x_fit);
        p = polyfit(x_fit, y_fit, 1);
        alpha = p(1);
        beta  = p(2);

        % Standard errors of slope and intercept
        y_hat     = polyval(p, x_fit);
        residuals = y_fit - y_hat;

        dof = nFit - 2;
        s2  = sum(residuals.^2) / dof;
        xbar = mean(x_fit);
        Sxx  = sum((x_fit - xbar).^2);

        var_a = s2 / Sxx;
        var_b = s2 * (1/nFit + xbar^2 / Sxx);

        alpha_std = sqrt(var_a);
        beta_std  = sqrt(var_b);
    else
        warning('Not enough valid data to fit scaling for %s.', typeName);
    end

    % Format alpha, beta with parentheses, like 0.0355(2), -4.48(3)
    nd_alpha = 4;
    nd_beta  = 2;

    alpha_fmt = format_with_paren(alpha, alpha_std, nd_alpha);
    beta_fmt  = format_with_paren(beta,  beta_std,  nd_beta);

    % Pack results
    res.type        = typeName;
    res.N           = sortedN;
    res.median_sec  = median_sec;
    res.err_sec     = err_sec;
    res.TTS_q_sec   = TTS_q_sec;
    res.alpha       = alpha;
    res.alpha_std   = alpha_std;
    res.alpha_fmt   = alpha_fmt;
    res.beta        = beta;
    res.beta_std    = beta_std;
    res.beta_fmt    = beta_fmt;
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

function N = parseN_from_prefix(pref)
    % Extract N from prefix patterns like "5R5X_n112", "7R7X_n96", etc.
    tok = regexp(pref, '_n(\d+)', 'tokens', 'once');
    if isempty(tok)
        error('Could not parse N from prefix "%s". Expected "_nNNN".', pref);
    end
    N = str2double(tok{1});
end

function s = format_with_paren(value, sigma, ndigits)
    % Format "value +/- sigma" as "v(s)" with ndigits decimal places.
    % Example: value=0.03547, sigma=0.00019, ndigits=4 -> "0.0355(2)".
    if isnan(value) || isnan(sigma)
        s = "";
        return;
    end
    v_round = round(value, ndigits);
    scale   = 10^ndigits;
    sigma_scaled = round(sigma * scale);  % integer in last-digit units
    fmt = sprintf('%%0.%df(%%d)', ndigits);
    s = sprintf(fmt, v_round, sigma_scaled);
end
