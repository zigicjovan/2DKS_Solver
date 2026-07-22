function [C_opt,delta_opt] = expfit_delta(v_radii, v_mean, prevdelta)
%   Estimate analyticity strip width delta from radial spectrum
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
%       4) fits log E(k) with log E ≈ log C - 2 * delta * k

    vcut = find(v_mean < 1e-15, 1, "first");
    if isempty(vcut)
        vcut = length(v_radii);
    end
    vstop = min(max(10,vcut),ceil(length(v_radii)-2));
    cutfracrange = 3:-0.5:3;
    cutchoices = NaN(length(cutfracrange),2);

    count = 1;
    for cutfrac = cutfracrange
        vstart = max(1,ceil(vstop/cutfrac));
        if cutfrac == cutfracrange(1)
            v_radii = v_radii(vstart:vstop);
            v_mean = v_mean(vstart:vstop);
        else
            cutstart = length(v_radii)-(vstop-vstart);
            v_radii = v_radii(cutstart:end);
            v_mean = v_mean(cutstart:end);
        end
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
    
        logE_fit = log(E_fit);              % Work in log-space
    
        p = polyfit(k_fit,logE_fit,1);      % exponential fit
        cutchoices(count,1)     = -0.5 * p(1);            % extract delta
        cutchoices(count,2)     = exp(p(2));              % extract coefficient
        count = count + 1;
    end

    idx = find( min(abs(prevdelta-cutchoices(1,:))) ,1);
    delta_opt = cutchoices(idx,1);
    C_opt = cutchoices(idx,2);
end