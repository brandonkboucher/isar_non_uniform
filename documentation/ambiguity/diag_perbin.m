% Per range-frequency-bin ambiguity: each bin has a perfect alias, but at its
% own fold distance Wx(f) = c*prf/(2(fc+f)*w0). Combining bins cancels them.
clear; clc
const = Constants; c = const.c;
fc = 1e9; prf = 200; w0 = pi; M = 16; L = 16; u0 = 1000;
t_m = (0:(1/prf):((M-1)/prf)).'; theta = w0*t_m;
x_true = -6.5395; y_true = -8.1511;

for fs = [300e6 10e6]
    df = fs/L; f_hat_l = (-L/2)*df:df:(L/2-1)*df;
    Wx_l = c*prf./(2*(fc+f_hat_l)*w0);          % fold distance per bin

    fprintf('\n=== fs = %.0f MHz (fractional BW %.3f) ===\n', fs/1e6, fs/fc);
    fprintf('fold distance across the band: %.2f m to %.2f m (spread %.2f m)\n', ...
        min(Wx_l), max(Wx_l), max(Wx_l)-min(Wx_l));

    x_test = linspace(x_true+0.5*min(Wx_l), x_true+1.5*max(Wx_l), 1201);

    % per-bin coherence (M samples each) and the full-band coherence
    Cb = zeros(L, numel(x_test));
    for l = 1:L
        a0 = exp(-1j*4*pi*(fc+f_hat_l(l))/c * rng_of(x_true, y_true, u0, theta));
        for i = 1:numel(x_test)
            ai = exp(-1j*4*pi*(fc+f_hat_l(l))/c * rng_of(x_test(i), y_true, u0, theta));
            Cb(l,i) = abs(ai'*a0)/M;
        end
    end
    Cfull = zeros(1,numel(x_test));
    a0f = compute_atom(x_true,y_true,u0,theta,f_hat_l,fc,false)/sqrt(M*L);
    for i = 1:numel(x_test)
        Cfull(i) = abs((compute_atom(x_test(i),y_true,u0,theta,f_hat_l,fc,false)/sqrt(M*L))'*a0f);
    end

    [~,ib] = max(Cb,[],2);
    fprintf('per-bin ghost peak: min %.3f, max %.3f (every bin still aliases perfectly)\n', ...
        min(max(Cb,[],2)), max(max(Cb,[],2)));
    fprintf('per-bin ghost position spans %.2f m; full-band ghost peak = %.3f\n', ...
        max(x_test(ib))-min(x_test(ib)), max(Cfull));

    f = figure('Visible','off','Position',[0 0 1500 950]);
    subplot(2,1,1)
    imagesc(x_test - x_true, (fc+f_hat_l)/1e9, Cb); axis xy
    hold on; plot(Wx_l, (fc+f_hat_l)/1e9, 'w--', 'LineWidth', 2); hold off
    colormap turbo; cb = colorbar; cb.Label.String = '|coherence| in that bin';
    xlabel('offset from true crossrange [m]'); ylabel('bin frequency f_c + f [GHz]')
    title(sprintf(['each range bin aliases perfectly, at its own fold distance ' ...
        '(dashed = W_x(f)) - f_s = %.0f MHz'], fs/1e6))
    subplot(2,1,2)
    plot(x_test - x_true, Cfull, 'LineWidth', 1.8); hold on
    xline(c*prf/(2*fc*w0), 'r--', 'LineWidth', 1.5); hold off
    ylim([0 1.05]); grid on
    xlabel('offset from true crossrange [m]'); ylabel('|coherence|')
    title(sprintf('all %d bins combined: ghost peak = %.3f (red = W_x at f_c)', L, max(Cfull)))
    saveas(f, sprintf(['/private/tmp/claude-501/-Users-brandonboucher-Documents-MATLAB-' ...
        'research-isar-non-uniform/c28f0005-26c2-492d-b21a-af98ba92c124/scratchpad/' ...
        'perbin_%dMHz.png'], round(fs/1e6)));
end

function r = rng_of(x, y, u0, theta)
    r = zeros(numel(theta),1);
    for m = 1:numel(theta)
        R = [cos(theta(m)) -sin(theta(m)); sin(theta(m)) cos(theta(m))];
        u = R*[x;y];
        r(m) = sqrt((u0+u(2))^2 + u(1)^2);
    end
end
