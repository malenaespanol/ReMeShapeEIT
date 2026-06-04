% This code produces the figures of the article:
% Qualitative reconstruction methods for imaging 
% interior Robin interfaces in EIT from Robin-to-Dirichlet data
% by Rafael Ceja Ayala, Malena I. Espanol, and Govanni Granados
% June 2026

function run_ReMeShapeEIT_CircularInterface()
    close all; clearvars -except; clc;
    rng(10);

    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');

    x = linspace(-1,1,70);

    examples = {
        'circle_rho0p4_noisefree',      0.4, 1, 1, 32, 0;
        'circle_rho0p4_noise1pct',      0.4, 1, 1, 32, 0.01;
        'circle_rho0p1_noise1pct',      0.1, 1, 1, 32, 0.01;
        'circle_rho0p4_two_sigmas_N32', 0.4, 1, 10, 32, 0.01;
        'circle_rho0p4_two_sigmas_N64', 0.4, 1, 10, 64, 0.01;
    };

    params_nf_M1  = struct('lambda',1e-6,  'alpha_sc',1e-6,  'ttls_k',5, 'norm_p',2);
    params_nf_M2  = struct('alpha_tik',1e-16, 'alpha_sc',1e-16, 'ttls_k',3, 'norm_p',2);
    params_ns_M1  = struct('lambda',1e-6,  'alpha_sc',1e-6,  'ttls_k',5, 'norm_p',2);
    params_ns_M2  = struct('alpha_tik',1e-16, 'alpha_sc',1e-16, 'ttls_k',3, 'norm_p',2);
    params_two_M1 = struct('lambda',1e-7, 'alpha_sc',1e-7, 'ttls_k',8, 'norm_p',2);
    params_two_M2 = struct('alpha_tik',0, 'alpha_sc',0, 'ttls_k',6, 'norm_p',2);

    for e = 1:size(examples,1)
        cname     = examples{e,1};
        rho       = examples{e,2};
        sigma_out = examples{e,3};
        sigma_in  = examples{e,4};
        Nang      = examples{e,5};
        delta     = examples{e,6};

        fprintf('Running example %d/%d: %s  (rho=%.3f, sig_out=%.3g, sig_in=%.3g, N=%d, delta=%.4g)\n',...
            e, size(examples,1), cname, rho, sigma_out, sigma_in, Nang, delta);

        theta = linspace(0,2*pi,Nang+1);
        ang   = theta(1:Nang);

        if contains(cname,'two_sigmas')
            paramsM1 = params_two_M1;
            paramsM2 = params_two_M2;
        elseif delta == 0
            paramsM1 = params_nf_M1;
            paramsM2 = params_nf_M2;
        else
            paramsM1 = params_ns_M1;
            paramsM2 = params_ns_M2;
        end

        A = zeros(Nang);
        for sidx = 1:Nang
            A(sidx,:) = ker(10, ang(sidx), ang, rho, 1, sigma_out, sigma_in);
        end
        A = (1/(2*pi*Nang)) .* A;

        E  = -1 + 2.*rand(Nang,Nang);
        E  = E / norm(E,2);
        Ah = A .* (ones(Nang,Nang) + delta * E);
        normAh = norm(Ah);
        Ah = Ah / normAh;

        [U,S,V] = svd(Ah,'econ');
        sigma_vals = diag(S);

        AT  = Ah.';
        ATA = AT * Ah;
        L   = eye(Nang);
        LTL = L.' * L;

        M_tik = ATA + paramsM1.lambda * LTL;

        methods = {'TIK','SC','TTLS'};
        nX = numel(x);
        M1 = struct('TIK',zeros(nX), 'SC',zeros(nX), 'TTLS',zeros(nX));
        M2 = struct('TIK',zeros(nX), 'SC',zeros(nX), 'TTLS',zeros(nX));

        make_b = @(xi,xj) arrayfun(@(kk) Greens_on_boundary(xi, xj, ang(kk), sigma_out), 1:Nang).';

        for ii = 1:nX
            for jj = 1:nX
                if norm([x(ii), x(jj)]) > 0.9
                    for m = 1:numel(methods)
                        M1.(methods{m})(ii,jj) = 0;
                        M2.(methods{m})(ii,jj) = 0;
                    end
                    continue
                end

                b = make_b(x(ii), x(jj)) / normAh;

                xa_tik = M_tik \ (AT * b);
                M1.TIK(ii,jj) = 1 / norm(xa_tik, paramsM1.norm_p);

                keepM1 = (sigma_vals.^2) >= paramsM1.alpha_sc;
                coeffM1 = zeros(length(sigma_vals),1);
                coeffM1(keepM1) = (U(:,keepM1)' * b) ./ sigma_vals(keepM1);
                xa_sc_M1 = V * coeffM1;
                M1.SC(ii,jj) = 1 / norm(xa_sc_M1, paramsM1.norm_p);

                kM1 = min(paramsM1.ttls_k, length(sigma_vals));
                xa_ttls_M1 = ttls(Ah, b, kM1);
                M1.TTLS(ii,jj) = 1 / norm(xa_ttls_M1, paramsM1.norm_p);

                M2.TIK(ii,jj)  = IndfuncRegFM_vec(U, S, b, paramsM2.alpha_tik, 1);
                M2.SC(ii,jj)   = IndfuncRegFM_vec(U, S, b, paramsM2.alpha_sc,  3);
                M2.TTLS(ii,jj) = IndfuncRegFM_ttls(A, b, U, S, paramsM2.ttls_k);
            end
        end

        xi = rho .* cos(theta);
        yi = rho .* sin(theta);
        titles = {'Tikhonov','Spectral cutoff','TTLS'};

        figure('Name',cname,'Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.8]);

        top_margin     = 0.92;
        bottom_margin  = 0.08;
        left_margin    = 0.05;
        right_margin   = 0.98;
        vertical_gap   = 0.10;
        horizontal_gap = 0.06;
        nCols          = 3;
        title_space    = 0.06;

        usable_top     = top_margin - title_space;
        subplot_height = (usable_top - bottom_margin - vertical_gap) / 2;
        subplot_width  = (right_margin - left_margin - (nCols-1)*horizontal_gap) / nCols;

        for m = 1:3
            left   = left_margin + (m-1)*(subplot_width + horizontal_gap);
            bottom = bottom_margin + subplot_height + vertical_gap;
            ax1    = axes('Position',[left, bottom, subplot_width, subplot_height]);
            plot_panel(M1.(methods{m}), x, xi, yi);
            title(['LSM: ' titles{m}], 'Interpreter','latex', 'FontSize',20);
            ax1.FontSize = 16;
            colorbar;

            bottom = bottom_margin;
            ax2    = axes('Position',[left, bottom, subplot_width, subplot_height]);
            plot_panel(M2.(methods{m}), x, xi, yi);
            title(['RFM: ' titles{m}], 'Interpreter','latex', 'FontSize',20);
            ax2.FontSize = 16;
            colorbar;
        end

        safe_name = strrep(cname,'_','\_');
        sgtitle(sprintf('Figure %d: %s \n\n', e, safe_name),'Interpreter','latex','FontSize',18);
    end

    fprintf('Done: all examples processed.\n');
end

function plot_panel(M, x, xi, yi)
    Mmax = max(M(:));
    if (~isfinite(Mmax)) || (Mmax == 0)
        Mnorm = M;
    else
        Mnorm = M / Mmax;
    end
    if (~any(isfinite(Mnorm(:)))) || (max(Mnorm(:)) - min(Mnorm(:)) == 0)
        imagesc(x, x, Mnorm.'); axis xy;
    else
        contourf(x, x, Mnorm.', 15, 'LineStyle','none');
    end
    colorbar
    hold on; plot(xi, yi, '--r','LineWidth',2); axis square; hold off
end

function y = IndfuncRegFM_vec(Uh, Sh, b, alpha, type)
    if nargin < 5, type = 1; end
    if type == 1
        f = @(t) (t.^2 ./ (alpha + t.^2)).^2;
    elseif type == 2
        f = @(t) (t ./ (alpha + t)).^2;
    elseif type == 3
        f = @(t) double((t.^2) >= alpha);
    elseif type == 4
        a = 1/(2*Sh(1,1)^2); m = floor(1/alpha);
        f = @(t) (1 - (1 - a*t.^2).^m).^2;
    else
        f = @(t) (t.^2 ./ (alpha + t.^2)).^2;
    end
    r    = size(Uh,2);
    Y    = zeros(1,r);
    for iy = 1:r
        t    = Sh(iy,iy);
        Y(iy) = (f(t) .* abs((b')*Uh(:,iy)).^2) ./ abs(t);
    end
    denom = sum(abs(Y));
    if denom == 0, y = 0; else, y = 1/denom; end
end

function x_ttls = ttls(A, b, k)
    [~, nA] = size(A);
    Ab      = [A, b];
    [~, ~, V] = svd(Ab, 'econ');
    k    = max(0, min(floor(k), size(V,2)-1));
    V12  = V(1:nA, k+1:end);
    V22  = V(nA+1, k+1:end-1);
    denom = norm(V22);
    if denom < eps, x_ttls = A \ b; return; end
    x_ttls = -(V12 * V22.') / (denom^2);
end

function y = IndfuncRegFM_ttls(A, b, U, S, k)
    s    = diag(S);
    r    = length(s);
    Utb  = U' * b;
    Utb2 = abs(Utb).^2;
    B    = [A, b];
    [~, Sbar, Vbar] = svd(B, 'econ');
    bar_sigma = diag(Sbar);
    k    = min(max(1, round(k)), length(bar_sigma));
    if k < size(Vbar,2)
        normV22_sq = sum(Vbar(end, k+1:end).^2);
    else
        normV22_sq = eps;
    end
    Bsq    = bar_sigma(1:k).^2;
    v_last = Vbar(end,:);
    w      = (v_last(1:k).^2)' / (normV22_sq + eps);
    f      = @(t) (t.^2) .* (w' * (1 ./ (Bsq - t.^2 + 1e-9)));
    yval   = 0;
    for i = 1:r
        denom = abs(s(i));
        if denom < eps, continue; end
        yval = yval + abs(f(s(i)) * Utb2(i) / denom);
    end
    if yval == 0, y = 0; else, y = 1/yval; end
end

function Gval = Greens_on_boundary(xi, xj, s_ang, sigma)
    zx     = xi; zy = xj;
    rho    = sqrt(zx^2 + zy^2);
    phi_z  = atan2(zy, zx);
    lambda = s_ang - phi_z;
    P      = @(t,lam) (1 - t.^2) ./ (1 - 2.*t.*cos(lam) + t.^2);
    integrand = @(s) s.^(1/sigma - 1) .* P(rho.*s, lambda);
    try
        val = integral(integrand, 0, 1, 'RelTol',1e-6, 'AbsTol',1e-8);
    catch
        ss  = linspace(0,1,300); ss(1) = 1e-9;
        val = trapz(ss, ss.^(1/sigma-1) .* (1-(rho.*ss).^2) ./ ...
              (1 - 2.*rho.*ss.*cos(lambda) + (rho.*ss).^2));
    end
    Gval = (1/(2*pi)) * val;
end

function K = ker(N, t, s, rho, gam, sig1, sig2)
    coef = @(n) (2*rho.^(2*abs(n)).*sig1.*abs(n).*(gam*rho + sig1*abs(n) - sig2*abs(n))) ./ ...
        ((sig1*abs(n)+1) .* (sig1.^2*abs(n).^2 - gam*rho + sig1*abs(n) + sig2*abs(n) ...
        + gam*rho.*rho.^(2*abs(n)) - rho.^(2*abs(n)).*sig1.^2.*abs(n).^2 ...
        + sig1*sig2.*abs(n).^2 + rho.^(2*abs(n)).*sig1.*abs(n) ...
        - rho.^(2*abs(n)).*sig2.*abs(n) + rho.^(2*abs(n)).*sig1.*sig2.*abs(n).^2 ...
        - gam*rho.*sig1.*abs(n) - gam*rho.*rho.^(2*abs(n)).*sig1.*abs(n)));
    k_p = zeros(N,1); k_n = zeros(N,1);
    for n = 1:N
        k_p(n) = coef(n);
        k_n(n) = coef(-n);
    end
    K = (sig1 - gam*rho*log(rho)) / (sig1 - gam*rho*log(rho) + gam*rho*sig1) - 1;
    for n = 1:N
        K = K + (k_p(n).*exp(1i*n*(t-s)) + k_n(n).*exp(1i*n*(s-t)));
    end
    K = real(K);
end