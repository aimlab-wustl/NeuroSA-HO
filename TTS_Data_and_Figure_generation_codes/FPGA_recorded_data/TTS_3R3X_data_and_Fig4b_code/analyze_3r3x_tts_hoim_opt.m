%% analyze_TTS_3R3X.m
% Analyze 3R3X TTS data and reproduce SBM-style scaling plots and table.

clear; clc;

%% Parameters
scale_factor = 1e-8;         % iterations -> seconds
quantiles_q  = [0.25, 0.50, 0.75];

% If you want to restrict which sizes N to use in the fit, set these indices.
fitIdx_for_all_q = [];       % [] => use all valid N for each q

% Small tolerance for matching thld_delta between tables (avoid FP issues)
delta_tol = 1e-12;

%% Load CSVs
summary = readtable('coarse_summary.csv');
inst    = readtable('coarse_instance_tts.csv');

% Ensure prefix columns are cell arrays of char
summary.prefix = normalize_prefix_column(summary.prefix);
inst.prefix    = normalize_prefix_column(inst.prefix);

%% Parse N from prefix "3R3X_n112" -> N = 112
getN = @(s) sscanf(s, '3R3X_n%d');

summary.N = cellfun(getN, summary.prefix);
inst.N    = cellfun(getN, inst.prefix);

%% Unique prefixes and their sizes, sorted by N
[prefixList, ~, ic] = unique(summary.prefix, 'stable');
N_per_prefix = accumarray(ic, summary.N, [], @(x) x(1));

[sortedN, sortIdx] = sort(N_per_prefix);
sortedPrefixes     = prefixList(sortIdx);

nSizes = numel(sortedN);

%% Containers
opt_thld_delta      = zeros(nSizes, 1);
opt_median_TTS_iter = zeros(nSizes, 1);  % median in iterations
median_err_iter     = zeros(nSizes, 1);  % error bar on median (iterations)

TTS_quantiles_iter  = nan(nSizes, numel(quantiles_q)); % per size, per q (iterations)

%% Main loop over sizes
for k = 1:nSizes
    pref = sortedPrefixes{k};
    Nval = sortedN(k);

    % ----- pick thld_delta that minimizes median_TTS_val for this prefix -----
    maskSum = strcmp(summary.prefix, pref);
    subSum  = summary(maskSum, :);

    if isempty(subSum)
        warning('No summary rows found for prefix %s (N=%d).', pref, Nval);
        continue;
    end

    [minMedian, idxMin]  = min(subSum.median_TTS_val);
    delta_star           = subSum.thld_delta(idxMin);
    opt_thld_delta(k)    = delta_star;
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

    % ====== NEW: treat NaN/Inf as large number for quantiles ======
    % First find finite TTS values
    finiteMask0 = isfinite(TTS_i_iter);
    if any(finiteMask0)
        maxFin   = max(TTS_i_iter(finiteMask0));
        largeVal = 10 * maxFin;        % you can change factor 10 if desired
    else
        % pathological case: no finite TTS at all
        largeVal = 1e12;               % arbitrary very large number
    end

    TTS_all_for_q = TTS_i_iter;
    TTS_all_for_q(~finiteMask0) = largeVal;
    % ===============================================================

    % ----- quantiles q = 0.25, 0.5, 0.75 (over ALL instances) -----
    % NOTE: NaN/Inf have been replaced by a large value, so unsolved
    % instances contribute correctly to large quantiles.
    for qi = 1:numel(quantiles_q)
        q = quantiles_q(qi);
        TTS_quantiles_iter(k, qi) = quantile(TTS_all_for_q, q);
    end

    % ----- error bar on median TTS using SciAdv S10 & S11 -----
    % For error bars we still only use well-defined (finite) TTS with 0<Ps<1.
    finiteMaskEB = isfinite(TTS_i_iter) & (Ps > 0) & (Ps < 1);
    TTS_finite   = TTS_i_iter(finiteMaskEB);
    Ps_finite    = Ps(finiteMaskEB);
    Nrep_finite  = Nrep(finiteMaskEB);

    if isempty(TTS_finite)
        warning('No finite TTS with 0 < Ps < 1 for prefix %s (N=%d).', pref, Nval);
        continue;
    end

    % ΔP_S = sqrt(P_S - P_S^2)/sqrt(N_rep)
    dPs = sqrt(Ps_finite - Ps_finite.^2) ./ sqrt(Nrep_finite);

    % ΔTTS = TTS * ΔP_S / ((1 - P_S)*ln(1 - P_S))
    denom  = (1 - Ps_finite) .* log(1 - Ps_finite);
    dTTS_i = TTS_finite .* (dPs ./ denom);

    % Approximate SE of the median: treat contributions as independent,
    % with weight ~1/n for each instance.
    nF              = numel(TTS_finite);
    se_median_iter  = sqrt(sum(dTTS_i.^2)) / nF;
    median_err_iter(k) = se_median_iter;
end

%% Convert everything to seconds
opt_median_TTS_sec = opt_median_TTS_iter * scale_factor;
median_err_sec     = median_err_iter     * scale_factor;
TTS_quantiles_sec  = TTS_quantiles_iter  * scale_factor;

%% Scaling fits: log10(TTS_q) = alpha_q * N + beta_q
nQ = numel(quantiles_q);
alphas     = nan(nQ, 1);
betas      = nan(nQ, 1);
alpha_std  = nan(nQ, 1);
beta_std   = nan(nQ, 1);

for qi = 1:nQ
    yq = TTS_quantiles_sec(:, qi);
    valid = isfinite(yq) & (yq > 0);

    x = sortedN(valid);
    y = log10(yq(valid));  % log10(TTS_q)

    if isempty(x) || numel(x) < 3
        warning('Not enough valid data to fit q = %.2f.', quantiles_q(qi));
        continue;
    end

    % Optionally restrict to intermediate sizes
    if ~isempty(fitIdx_for_all_q)
        validIdx = fitIdx_for_all_q;
        validIdx = validIdx(validIdx <= numel(x));
        x_fit = x(validIdx);
        y_fit = y(validIdx);
    else
        x_fit = x;
        y_fit = y;
    end

    % Linear regression: y = alpha * N + beta
    nFit = numel(x_fit);
    p = polyfit(x_fit, y_fit, 1);
    a = p(1);
    b = p(2);

    % Standard errors of slope and intercept
    y_hat     = polyval(p, x_fit);
    residuals = y_fit - y_hat;

    dof = nFit - 2;
    s2  = sum(residuals.^2) / dof;
    xbar = mean(x_fit);
    Sxx  = sum((x_fit - xbar).^2);

    var_a = s2 / Sxx;
    var_b = s2 * (1/nFit + xbar^2 / Sxx);

    alphas(qi)    = a;
    betas(qi)     = b;
    alpha_std(qi) = sqrt(var_a);
    beta_std(qi)  = sqrt(var_b);
end

%% Plot: median TTS vs N with error bars and fitted scaling line (q = 0.5)
figure; clf; hold on;
dummy1 = repmat(char(160), 1, 36);
% Median points with error bars
errorbar(sortedN, opt_median_TTS_sec, median_err_sec, ...
    'o', 'LineWidth', 1.8, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'DisplayName', 'Median TTS');

% Fitted line for q = 0.5
idx_med = find(abs(quantiles_q - 0.50) < 1e-9, 1);
alpha_med = alphas(idx_med);
beta_med  = betas(idx_med);

if ~isnan(alpha_med)
    % Line to plot
    N_line  = linspace(min(sortedN), max(sortedN), 200);
    TTS_fit = 10.^(alpha_med * N_line + beta_med);

    % Use the formatted values with parentheses in the legend, like 0.0355(2)
    nd_alpha_plot = 4;  % same as nd_alpha
    nd_beta_plot  = 2;  % same as nd_beta

    alpha_med_str = format_with_paren(alpha_med, alpha_std(idx_med), nd_alpha_plot);
    beta_med_str  = format_with_paren(beta_med,  beta_std(idx_med),  nd_beta_plot);

    fitLabel = sprintf('Fit: $ \\alpha=%s, \\beta=%s$', ...
                       alpha_med_str, beta_med_str);

    plot(N_line, TTS_fit, '-', 'LineWidth', 1.8, ...
        'DisplayName', fitLabel);
end


set(gca, 'YScale', 'log', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...
    'TickLabelInterpreter', 'latex');

%xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
%ylabel('Optimal median TTS (s)', 'Interpreter', 'latex', 'FontSize', 18);

box on;
xlim([min(sortedN)-5, max(sortedN)+5]);

%legend('Location', 'northwest', 'Interpreter', 'latex');
legend('off');
%title('3R3X optimal median TTS vs problem size', ...
%    'Interpreter', 'latex', 'FontSize', 18);

%% Prepare table with alpha, beta for q = 0.25, 0.5, 0.75
nd_alpha = 4;   % digits after decimal for alpha
nd_beta  = 2;   % digits after decimal for beta

alpha_fmt = strings(nQ, 1);
beta_fmt  = strings(nQ, 1);

for qi = 1:nQ
    alpha_fmt(qi) = format_with_paren(alphas(qi), alpha_std(qi), nd_alpha);
    beta_fmt(qi)  = format_with_paren(betas(qi),  beta_std(qi),  nd_beta);
end

T_out = table(quantiles_q(:), ...
              alphas, alpha_std, alpha_fmt, ...
              betas,  beta_std,  beta_fmt, ...
    'VariableNames', {'q', 'alpha', 'alpha_std', 'alpha_fmt', ...
                          'beta',  'beta_std',  'beta_fmt'});

writetable(T_out, 'alpha_beta_TTS.xlsx');

fprintf('Done. Figure created and alpha_beta_TTS.xlsx written.\n');

%% ---- Helper functions ---------------------------------------------
function prefOut = normalize_prefix_column(prefIn)
    % Normalize a prefix column from readtable into a cell array of char.
    if isstring(prefIn)
        prefOut = cellstr(prefIn);
    elseif iscell(prefIn)
        if ischar(prefIn{1})
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

function s = format_with_paren(value, sigma, ndigits)
    % Format "value +/- sigma" as "v(s)" with ndigits decimal places.
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
