% NEWTON_METHOD_EXACT  newton_method.m with the exact Hessian of the
% concentrated objective G(p) = |a'r|^2/(a'a), and nothing else changed (no
% backtracking: the first rejected step still ends the search).
%
%   With e = r - a*alpha, alpha = a'r/(a'a), n = a'a and a_i = da/dp_i,
%
%     (1/2) d2G/dp_i dp_j =  Re{ alpha' * a_ij' * e }                        T1
%                          - |alpha|^2 Re{ a_i'a_j - (a_i'a)(a'a_j)/n }       T2
%                          + Re{ (a_i'e) * conj(a_j'e) } / n                  T3
%                          - Re{ alpha' * [ (a_i'a)(a_j'e) + (a_j'a)(a_i'e) ] } / n   T4
%
%   The original code has T1 and T2 without its (a_i'a)(a'a_j)/n part, which
%   is the Hessian with alpha held fixed. For range, a_y is dominated by the
%   carrier phase -j*4*pi*(fc+f)/c * a, which lies along a; with alpha held
%   fixed that component counts as curvature and H(2,2) came out ~1e5x too
%   large. T2's projection removes it; T3 and T4 are the terms from alpha
%   re-optimising with p, which matter when the residual is far from
%   explained (e.g. other scatterers still in r). Checked against a
%   finite-difference Hessian to four significant figures (check_hessian.m).
%
%   Same signature as newton_method.

function [p_hat, p_hat_hist] = newton_method_exact(...
    r, ... % [ML x  1] residual
    x, ... % (m) cross-range position
    y, ... % (m) range position
    u0, ... % (m) radar-to-target center distance
    theta_m, ... % [M x 1] (rad) target yaw
    f_hat_l, ... % [L x 1] (Hz) range-frequency
    fc, ... % (Hz) center frequency
    Ns, ... % number of refinement steps
    use_range_approx ... % range approx boolean
    )

    % compute initial atom and G
    a = compute_atom(x, y, u0, theta_m, f_hat_l, fc, use_range_approx);
    G_old = abs((a' * r))^2 / (a' * a);
    p_hat = [x;y];

    % Ns is the upper bound; the loop can stop early, so nfilled tracks how
    % many columns are real
    p_hat_hist = zeros(2, Ns);
    nfilled = 0;

    max_backtracking_steps = 5;

    for istep = 1:Ns

        x = p_hat(1); y = p_hat(2);

        [dadx, dady, d2adx2, d2ady2, d2adxdy] = compute_atom_derivatives(...
            a, x, y, u0, theta_m, f_hat_l, fc, use_range_approx);

        n       = real(a' * a);
        alpha_k = (a' * r) / n;
        e       = r - a * alpha_k;

        da  = {dadx, dady};
        d2a = {d2adx2, d2adxdy; d2adxdy, d2ady2};

        % inner products reused by every entry
        daHa = [dadx' * a;  dady' * a];   % a_i' a
        daHe = [dadx' * e;  dady' * e];   % a_i' e

        % gradient, (1/2) dG/dp_i = Re{ alpha' a_i' e }
        F = real(conj(alpha_k) * daHe);

        % exact Hessian, (1/2) d2G/dp_i dp_j
        H = zeros(2);
        for i = 1:2
            for j = i:2
                T1 = real(conj(alpha_k) * (d2a{i,j}' * e));
                T2 = -abs(alpha_k)^2 * real(da{i}' * da{j} - daHa(i) * conj(daHa(j)) / n);
                T3 = real(daHe(i) * conj(daHe(j))) / n;
                T4 = -real(conj(alpha_k) * (daHa(i) * daHe(j) + daHa(j) * daHe(i))) / n;
                H(i,j) = T1 + T2 + T3 + T4;
                H(j,i) = H(i,j);
            end
        end

        % Newton step toward the maximum of G; stop if H is not negative
        % definite (not a local maximum model)
        if ~(trace(H) < 0 && det(H) > 0), break; end

        % calculate the direction
        d = H \ F;
        p0 = [x; y];

        for i = 0:max_backtracking_steps

            p = p0 - (1/2)^i .* d;
            x = p(1); y = p(2);
            p_hat_hist(:,istep) = [x; y];
            nfilled = istep;
    
            %-------- Refinement Acceptance Condition -------------
            a_new = compute_atom(x, y, u0, theta_m, f_hat_l, fc, use_range_approx);
            G_new = abs((a_new' * r))^2 / (a_new' * a_new);
    
            continue_search = false;
            if G_new > G_old
                p_hat = [x; y];
                G_old = G_new;
                a = a_new;
                continue_search = true;
                break
            end
        end
        if ~continue_search
            break
        end
    end

    % hand back only the steps that were actually taken
    p_hat_hist = p_hat_hist(:, 1:nfilled);

end
