function delta_opt = fit_delta(v_radii, v_mean)
%FIT_DELTA  Estimate analyticity strip width delta from radial spectrum
%
%   delta_opt = FIT_DELTA(v_radii, v_mean)
%
%   Inputs:
%       v_radii : vector of radial wavenumbers k (shell centers)
%       v_mean  : vector of radial energy spectrum E(k)
%
%   Output:
%       delta_opt  : best-fit analyticity strip width delta, assuming
%                    E(k) ~ C * exp(-2 * delta * k) on a suitable k-range.
%
%   Notes:
%   - This function:
%       1) cleans the data (positive energies, finite values),
%       2) discards points below an energy threshold (noise floor),
%       3) selects a mid-range of k indices (~10%–90% of valid radii),
%       4) fits log E(k) with log E ≈ a - 2 * delta * k
%          by minimizing least-squares error over delta using fminsearch.
%       5) For each delta, the optimal a = log C is computed analytically,
%          so only delta is passed to fminsearch (1D optimization).

    vcut = find(v_mean < 1e-15, 1, "first");
    if isempty(vcut)
        vcut = length(v_radii);
    end
    vstop = min(max(10,vcut),ceil(length(v_radii)-5));
    v_radii = v_radii(1:vstop);
    v_mean = v_mean(1:vstop);

    % Ensure column vectors
    k_all = v_radii(:);
    E_all = v_mean(:);

    % Keep only positive, finite values
    mask = (E_all > 0) & isfinite(E_all) & isfinite(k_all);
    k_all = k_all(mask);
    E_all = E_all(mask);

    if numel(k_all) < 5
        error('Not enough valid (k,E) points after initial cleaning.');
    end

    % Sort by k (just in case)
    [k_sorted, idx] = sort(k_all);
    E_sorted = E_all(idx);

    % -------------------------------------------------------------
    % Energy threshold: remove points on/near the numerical noise floor
    % -------------------------------------------------------------
    Emax = max(E_sorted);
    if Emax <= 0
        error('All energies are non-positive after cleaning.');
    end

    % Threshold factor can be tuned
    thresh_factor = 1e-12;
    E_thresh = thresh_factor * Emax;

    mask_thr = (E_sorted > E_thresh);
    if sum(mask_thr) >= 5
        k_use = k_sorted(mask_thr);
        E_use = E_sorted(mask_thr);
    else
        % Fallback: not enough points above threshold, use all data
        k_use = k_sorted;
        E_use = E_sorted;
    end

    if numel(k_use) < 5
        error('Not enough (k,E) points remain after applying energy threshold.');
    end

    % -------------------------------------------------------------
    % Choose a window of available k, filtering out higher wavevectors
    % -------------------------------------------------------------
    N = numel(k_use);
    i1 = 1;
    i2 = max(i1+2, round(1.0 * N));   % ensure at least a few points

    k_fit = k_use(i1:i2);
    E_fit = E_use(i1:i2);

    % Work in log-space
    logE_fit = log(E_fit);

    % 1D objective: only delta is optimized.
    % Model: logE ≈ a - 2 * delta * k
    % For fixed delta, optimal a is:
    %     a(delta) = mean( logE_fit + 2 * delta * k_fit )
    % Residuals: r = logE_fit - (a(delta) - 2 * delta * k_fit)
    % Objective: J(delta) = sum( r.^2 )
    obj = @(delta) objective_delta(delta, k_fit, logE_fit);

    % Initial guess for delta: rough slope estimate from first/last points
    dk    = k_fit(end)    - k_fit(1);
    dlogE = logE_fit(end) - logE_fit(1);
    if dk ~= 0
        % slope ~ dlogE/dk ≈ -2 * delta  => delta ~ -0.5 * slope
        delta0 = max(1e-3, -0.5 * dlogE / dk);
    else
        delta0 = 0.2;  % fallback
    end

    % Run fminsearch
    options = optimset('Display','off', 'TolX',1e-6, 'TolFun',1e-10);
    delta_opt = fminsearch(obj, delta0, options);

end

% -------------------------------------------------------------------------
function J = objective_delta(delta, k_fit, logE_fit)
%OBJECTIVE_DELTA  Least-squares objective for a given delta.
%
%   For each delta, compute a(delta) = mean(logE + 2 * delta * k),
%   then residual r = logE - (a(delta) - 2 * delta * k),
%   and return J(delta) = sum(r.^2).

    % Penalize non-physical negative or zero delta
    if delta <= 0
        J = 1e12 + (abs(delta) + 1)^2;
        return;
    end

    a = mean(logE_fit + 2 * delta * k_fit);   % best logC for this delta
    r = logE_fit - (a - 2 * delta * k_fit);
    J = sum(r.^2);

end
