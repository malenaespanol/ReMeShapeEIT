% run_ReMeShapeEIT_Star.m
%
% This code produces the figures of the article:
% Qualitative reconstruction methods for imaging
% interior Robin interfaces in EIT from Robin-to-Dirichlet data
% by Rafael Ceja Ayala, Malena I. Espanol, and Govanni Granados
% June 2026

function run_ReMeShapeEIT_StarInterface()
    close all; clearvars -except; clc;
    rng(10);

    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    x = linspace(-1,1,70);

    Mcoupled = 16;

    shape_params = struct('r0',0.4, 'a',0.25, 'k',5);

    examples = {
        'star_noise1pct', 1, 1, 32, 0.01;
    };

    params_ns_M1 = struct('lambda',1e-7,  'alpha_sc',1e-7,  'alpha_ttls',1/8, 'norm_p',2);
    params_ns_M2 = struct('alpha_tik',0,  'alpha_sc',0,     'alpha_ttls',1/6, 'norm_p',2);

    [rho_fn, drho_fn] = make_interface(shape_params);

    for e = 1:size(examples,1)
        cname     = examples{e,1};
        sigma_out = examples{e,2};
        sigma_in  = examples{e,3};
        Nang      = examples{e,4};
        delta     = examples{e,5};

        fprintf('Running example %d/%d: %s  (sig_out=%.3g, sig_in=%.3g, N=%d, delta=%.4g)\n',...
            e, size(examples,1), cname, sigma_out, sigma_in, Nang, delta);

        theta = linspace(0,2*pi,Nang+1);
        ang   = theta(1:Nang);

        paramsM1 = params_ns_M1;
        paramsM2 = params_ns_M2;

        fprintf('  Building A via coupled-mode solver (M=%d) ...\n', Mcoupled);
        tic;
        A = build_A_coupled(rho_fn, drho_fn, ang, sigma_out, sigma_in, 1, Mcoupled);
        fprintf('  Done (%.1f s)\n', toc);

        E       = -1 + 2.*rand(Nang,Nang);
        E       = E / norm(E,2);
        Ah      = A .* (ones(Nang,Nang) + delta*E);
        A_scale = norm(Ah,2);
        Ahn     = Ah / A_scale;

        [U,S,V]  = svd(Ahn,'econ');
        sigma_vals = diag(S);
        AT       = Ahn.';
        M_tik_M1 = AT*Ahn + paramsM1.lambda*eye(Nang);

        methods = {'TIK','SC','TTLS'};
        nX = numel(x);
        M1 = struct('TIK',zeros(nX),'SC',zeros(nX),'TTLS',zeros(nX));
        M2 = struct('TIK',zeros(nX),'SC',zeros(nX),'TTLS',zeros(nX));

        make_b = @(xi,xj) arrayfun(@(kk) Greens_on_boundary(xi,xj,ang(kk),sigma_out),1:Nang).' / A_scale;

        for ii = 1:nX
            for jj = 1:nX
                if norm([x(ii),x(jj)]) > 0.9
                    for mIdx = 1:numel(methods)
                        M1.(methods{mIdx})(ii,jj) = 0;
                        M2.(methods{mIdx})(ii,jj) = 0;
                    end
                    continue
                end

                b = make_b(x(ii),x(jj)) / A_scale;

                xa_tik        = M_tik_M1 \ (AT*b);
                M1.TIK(ii,jj) = 1 / norm(xa_tik, paramsM1.norm_p);

                keepM1          = (sigma_vals.^2 >= paramsM1.alpha_sc);
                coeffM1         = zeros(length(sigma_vals),1);
                coeffM1(keepM1) = (U(:,keepM1)'*b) ./ sigma_vals(keepM1);
                M1.SC(ii,jj)    = 1 / norm(V*coeffM1, paramsM1.norm_p);

                kM1            = min(round(1/paramsM1.alpha_ttls), length(sigma_vals));
                M1.TTLS(ii,jj) = 1 / norm(ttls(Ahn,b,kM1), paramsM1.norm_p);

                M2.TIK(ii,jj)  = IndfuncRegFM_vec(U, S, b, paramsM2.alpha_tik, 1);
                M2.SC(ii,jj)   = IndfuncRegFM_vec(U, S, b, paramsM2.alpha_sc,  3);

                kM2            = min(round(1/paramsM2.alpha_ttls), length(sigma_vals));
                M2.TTLS(ii,jj) = IndfuncRegFM_ttls(Ahn, b, U, S, kM2);
            end
        end

        xi_plot = arrayfun(@(t) rho_fn(t)*cos(t), theta);
        yi_plot = arrayfun(@(t) rho_fn(t)*sin(t), theta);

        titles = {'Tikhonov','Spectral cutoff','TTLS'};
        figure('Name',cname,'Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.8]);

        top_margin     = 0.92;
        bottom_margin  = 0.08;
        left_margin    = 0.05;
        right_margin   = 0.98;
        vertical_gap   = 0.10;
        horizontal_gap = 0.06;
        title_space    = 0.06;

        usable_top     = top_margin - title_space;
        subplot_height = (usable_top - bottom_margin - vertical_gap) / 2;
        subplot_width  = (right_margin - left_margin - 2*horizontal_gap) / 3;

        for mIdx = 1:3
            left = left_margin + (mIdx-1)*(subplot_width + horizontal_gap);
            ax1  = axes('Position',[left, bottom_margin+subplot_height+vertical_gap, subplot_width, subplot_height]);
            plot_panel(M1.(methods{mIdx}), x, xi_plot, yi_plot);
            title(['LSM: ' titles{mIdx}], 'Interpreter','latex', 'FontSize',20);
            ax1.FontSize = 16; colorbar;

            ax2  = axes('Position',[left, bottom_margin, subplot_width, subplot_height]);
            plot_panel(M2.(methods{mIdx}), x, xi_plot, yi_plot);
            title(['RFM: ' titles{mIdx}], 'Interpreter','latex', 'FontSize',20);
            ax2.FontSize = 16; colorbar;
        end

        safe_name = strrep(cname,'_','\_');
        sgtitle(sprintf('Figure %d: %s\n\n', e, safe_name), 'Interpreter','latex', 'FontSize',18);
    end

    fprintf('Done.\n');
end

function [rho_fn, drho_fn] = make_interface(p)
    rho_fn  = @(t)  p.r0 * (1 + p.a*cos(p.k.*t));
    drho_fn = @(t) -p.r0 * p.a * p.k * sin(p.k.*t);
end

function A = build_A_coupled(rho_fn, drho_fn, ang, sig_out, sig_in, gam, M)
    Nang  = length(ang);
    Ntot  = 2*M+1;
    n_idx = (-M:M)';

    Nq    = max(4*M+4, 64);
    tq    = linspace(0,2*pi,Nq+1); tq = tq(1:Nq);
    rq    = arrayfun(rho_fn,  tq);
    drhot = arrayfun(drho_fn, tq);

    Tx = drhot.*cos(tq) - rq.*sin(tq);
    Ty = drhot.*sin(tq) + rq.*cos(tq);
    ds = sqrt(Tx.^2 + Ty.^2);
    nx =  Ty./ds;
    ny = -Tx./ds;

    r_dot_n =  cos(tq).*nx + sin(tq).*ny;
    t_dot_n = -sin(tq).*nx + cos(tq).*ny;

    absn_row = abs(n_idx)';
    rq_col   = rq(:);
    rq_safe  = max(rq_col, 1e-14);
    rq_pown  = rq_col  .^ absn_row;
    rq_nown  = rq_safe .^ (-absn_row);

    E_phase   = exp(1i*n_idx*tq);
    Phi_plus  = (rq_pown.*E_phase')';
    Phi_minus = (rq_nown.*E_phase')';
    Phi_in    = Phi_plus;

    dPhi_plus_n  = zeros(Ntot,Nq);
    dPhi_minus_n = zeros(Ntot,Nq);
    dPhi_in_n    = zeros(Ntot,Nq);

    for ki = 1:Ntot
        n = n_idx(ki); an = abs(n); phase = exp(1i*n*tq);
        if an == 0
            dPhi_plus_n(ki,:)  = zeros(1,Nq);
            dPhi_minus_n(ki,:) = (1./rq).*r_dot_n;
            dPhi_in_n(ki,:)    = zeros(1,Nq);
        else
            rp_m1 = rq.^(an-1); rm_m1 = rq.^(-(an+1));
            dPhi_plus_n(ki,:)  = rp_m1.*(an.*r_dot_n + 1i*n.*t_dot_n).*phase;
            dPhi_minus_n(ki,:) = rm_m1.*(-an.*r_dot_n + 1i*n.*t_dot_n).*phase;
            dPhi_in_n(ki,:)    = dPhi_plus_n(ki,:);
        end
    end

    wq   = (2*pi/Nq)*ds;
    Psi  = exp(-1i*n_idx*tq);
    PsiW = Psi .* (wq/(2*pi));

    F = zeros(3*Ntot, 3*Ntot);
    for ki = 1:Ntot
        n = n_idx(ki); an = abs(n);
        if an == 0
            F(ki,ki) = 1; F(ki,Ntot+ki) = sig_out;
        else
            F(ki,ki) = sig_out*an+1; F(ki,Ntot+ki) = -sig_out*an+1;
        end
    end

    Gram_plus  = PsiW * Phi_plus.';
    Gram_minus = PsiW * Phi_minus.';
    Gram_in    = PsiW * Phi_in.';

    F(Ntot+1:2*Ntot, 1:Ntot)        =  Gram_plus;
    F(Ntot+1:2*Ntot, Ntot+1:2*Ntot) =  Gram_minus;
    F(Ntot+1:2*Ntot, 2*Ntot+1:end)  = -Gram_in;

    dGram_plus_out  = PsiW * dPhi_plus_n.';
    dGram_minus_out = PsiW * dPhi_minus_n.';
    dGram_in_sig    = PsiW * dPhi_in_n.';

    F(2*Ntot+1:end, 1:Ntot)        = sig_out*dGram_plus_out  - gam*Gram_plus;
    F(2*Ntot+1:end, Ntot+1:2*Ntot) = sig_out*dGram_minus_out - gam*Gram_minus;
    F(2*Ntot+1:end, 2*Ntot+1:end)  = -sig_in*dGram_in_sig;

    A = zeros(Nang,Nang);
    for p = 1:Nang
        f_coeffs = (1/(2*pi))*exp(-1i*n_idx*ang(p));
        rhs      = zeros(3*Ntot,1); rhs(1:Ntot) = f_coeffs;
        coeff    = F \ rhs;
        a_vec    = coeff(1:Ntot);
        b_vec    = coeff(Ntot+1:2*Ntot);
        for ki = 1:Ntot
            n = n_idx(ki);
            if abs(n) == 0
                u_mode_at_bdy = a_vec(ki);
            else
                u_mode_at_bdy = a_vec(ki) + b_vec(ki);
            end
            u0_mode_at_bdy = f_coeffs(ki) / (sig_out*abs(n)+1);
            A(:,p) = A(:,p) + (u_mode_at_bdy - u0_mode_at_bdy)*exp(1i*n*ang');
        end
    end
    A = real(A) * (1/(2*pi*Nang));
end

function y = IndfuncRegFM_vec(Uh, Sh, b, alpha, type)
    if nargin < 5, type = 1; end
    switch type
        case 1,   f = @(t) (t.^2 ./ (alpha + t.^2)).^2;
        case 2,   f = @(t) (t ./ (alpha + t)).^2;
        case 3,   f = @(t) double((t.^2) >= alpha);
        case 4,   a = 1/(2*Sh(1,1)^2); mv = floor(1/alpha); f = @(t) (1-(1-a*t.^2).^mv).^2;
        otherwise, f = @(t) (t.^2 ./ (alpha + t.^2)).^2;
    end
    r = size(Uh,2); Y = zeros(1,r);
    for iy = 1:r
        t    = Sh(iy,iy);
        Y(iy) = f(t) .* abs(b'*Uh(:,iy)).^2 ./ abs(t);
    end
    denom = sum(abs(Y));
    if denom == 0, y = 0; else, y = 1/denom; end
end

function x_ttls = ttls(A, b, k)
    [~, nA] = size(A);
    Ab      = [A, b];
    [~, ~, V] = svd(Ab, 'econ');
    k     = max(0, min(floor(k), size(V,2)-1));
    V12   = V(1:nA, k+1:end);
    V22   = V(nA+1, k+1:end);
    denom = norm(V22);
    if denom < eps, x_ttls = A \ b; return; end
    x_ttls = -(V12 * V22.') / (denom^2);
end

function y = IndfuncRegFM_ttls(A, b, U, S, k)
    s    = diag(S); r = length(s);
    Utb  = U'*b; Utb2 = abs(Utb).^2;
    B    = [A, b];
    [~, Sbar, Vbar] = svd(B, 'econ');
    bar_sigma = diag(Sbar);
    k    = min(max(1,round(k)), length(bar_sigma));
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
    rho_z  = sqrt(xi^2 + xj^2);
    phi_z  = atan2(xj, xi);
    lambda = s_ang - phi_z;
    P      = @(t,lam) (1 - t.^2) ./ (1 - 2.*t.*cos(lam) + t.^2);
    integrand = @(s) s.^(1/sigma-1) .* P(rho_z.*s, lambda);
    try
        val = integral(integrand, 0, 1, 'RelTol',1e-6, 'AbsTol',1e-8);
    catch
        ss  = linspace(0,1,300); ss(1) = 1e-9;
        val = trapz(ss, ss.^(1/sigma-1) .* (1-(rho_z.*ss).^2) ./ ...
              (1 - 2.*rho_z.*ss.*cos(lambda) + (rho_z.*ss).^2));
    end
    Gval = (1/(2*pi)) * val;
end

function plot_panel(M, x, xi, yi)
    Mmax = max(M(:));
    if (~isfinite(Mmax)) || (Mmax == 0), Mnorm = M; else, Mnorm = M/Mmax; end
    if (~any(isfinite(Mnorm(:)))) || (max(Mnorm(:)) - min(Mnorm(:)) == 0)
        imagesc(x, x, Mnorm.'); axis xy;
    else
        contourf(x, x, Mnorm.', 15, 'LineStyle','none');
    end
    colorbar; hold on; plot(xi, yi, '--r', 'LineWidth',2); axis square; hold off
end