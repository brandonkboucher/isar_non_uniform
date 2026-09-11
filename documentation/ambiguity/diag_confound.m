% Crossrange ambiguity function |<a(x), a(x_true)>| across three ambiguities,
% for the 2x2 design {wideband, narrowband} x {uniform theta, accelerating}.
clear; clc
const = Constants; c = const.c;

fc = 1e9; prf = 200; w0 = pi; M = 16; L = 16; u0 = 1000;
t_m = (0:(1/prf):((M-1)/prf)).';
T   = M/prf;
Wx  = c*prf/(2*fc*w0);

x_true = -6.5395; y_true = -8.1511;      % the ambiguity -1 scatterer

w1_acc = 5;                              % [rad/s/s] modest yaw acceleration
fprintf('aperture rotation  w0*T           = %.4f rad (%.2f deg)\n', w0*T, rad2deg(w0*T));
fprintf('extra rotation     0.5*w1*T^2     = %.4f rad (%.2f deg), %.1f%% of the aperture\n', ...
    0.5*w1_acc*T^2, rad2deg(0.5*w1_acc*T^2), 100*(0.5*w1_acc*T^2)/(w0*T));
fprintf('Wx = %.4f m\n\n', Wx);

cfg = struct( ...
  'name', {'wideband, uniform', 'wideband, accelerating', ...
           'narrowband, uniform', 'narrowband, accelerating'}, ...
  'fs',   {300e6, 300e6, 10e6, 10e6}, ...
  'w1',   {0, w1_acc, 0, w1_acc});

x_test = linspace(x_true-1.6*Wx, x_true+1.6*Wx, 1601);
f = figure('Visible','off','Position',[0 0 1600 1100]);

for ic = 1:numel(cfg)
    df = cfg(ic).fs/L;
    f_hat_l = (-L/2)*df:df:(L/2-1)*df;
    theta_m = w0*t_m + 0.5*cfg(ic).w1*t_m.^2;

    a_true = compute_atom(x_true, y_true, u0, theta_m, f_hat_l, fc, false)/sqrt(M*L);
    cc = zeros(size(x_test));
    for i = 1:numel(x_test)
        ai = compute_atom(x_test(i), y_true, u0, theta_m, f_hat_l, fc, false)/sqrt(M*L);
        cc(i) = abs(ai'*a_true);
    end

    % coherence at the exact +/- Wx doppelgangers
    gh = zeros(1,2); sh = [-1 1];
    for j = 1:2
        ad = compute_atom(x_true+sh(j)*Wx, y_true, u0, theta_m, f_hat_l, fc, false)/sqrt(M*L);
        gh(j) = abs(ad'*a_true);
    end
    fprintf('%-26s  frac BW %5.3f | ghost(-Wx) %.3f  ghost(+Wx) %.3f\n', ...
        cfg(ic).name, cfg(ic).fs/fc, gh(1), gh(2));

    subplot(2,2,ic)
    plot(x_test, cc, 'LineWidth', 1.6); hold on
    xline(x_true, 'k-',  'LineWidth', 1.2);
    xline(x_true-Wx, 'r--', 'LineWidth', 1.2);
    xline(x_true+Wx, 'r--', 'LineWidth', 1.2);
    hold off
    ylim([0 1.05]); grid on
    xlabel('crossrange of trial atom [m]'); ylabel('|<a(x), a(x_{true})>|')
    title(sprintf('%s  (ghosts: %.2f / %.2f)', cfg(ic).name, gh(1), gh(2)), ...
        'FontSize', 12)
end
sgtitle({'Doppler-ambiguity coherence: solid = true scatterer, dashed = x_{true} \pm W_x', ...
         'the alias only exists when the dashed lines reach ~1'}, 'FontSize', 13);
saveas(f, ['/private/tmp/claude-501/-Users-brandonboucher-Documents-MATLAB-research-isar-non-uniform/' ...
    'c28f0005-26c2-492d-b21a-af98ba92c124/scratchpad/ambiguity_confound.png']);
