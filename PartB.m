%%
clear; close all; clc
exportFigures = false;

%% System Definition
M = 1; L_pend = 1; F = 1; g = 9.81; alpha = 1;

A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(M*L_pend) g/L_pend 0];

B = [0; 1/M; 0; -1/(M*L_pend)];

P = [0 0 0; 1/M 0 0; 0 0 0; -1/(M*L_pend) 0 0];

Ce = [1 0 0 0];
Qe = [0 -1 0];

%% Dynamics Functions
linear_dynamics = @(t, x, u, d1) A * x + B * u + P(:,1) * d1;

non_lin_dynamics = @(t, x, u, d1) [
    x(2);
    (u + d1 - F*x(2)) / M;
    x(4);
    (g/L_pend)*sin(x(3)) - ((u + d1 - F*x(2)) / (M * L_pend)) * cos(x(3))
];

%% Simulation
omegas = [0.1, 1 ,10];
n_periods = 3;
skip_periods = 0; 
period = 50;  % square wave period
initial_state = zeros(4,1);
n_omega = length(omegas);

%% Run simulations and create separate figures for B2, B3, B4

for i = 1:n_omega
    omega_i = omegas(i);
    T_sim = n_periods * (2*pi / omega_i);
    T_start = skip_periods * (2*pi / omega_i);  % Time to skip for steady state

    [K, L_ctrl, Pi, Gamma] = compute_control_gains(omega_i, A, B, P, Ce, Qe);

    [t_lin, ~, u_lin, y_lin, d_all_lin] = run_simulation(T_sim, omega_i, ...
        alpha, period, linear_dynamics, K, L_ctrl, initial_state);

    [t_nl, ~, u_nl, y_nl, d_all_nl] = run_simulation(T_sim, omega_i, ...
        alpha, period, @(t,x,u,d1) non_lin_dynamics(t,x,u,d1), K, L_ctrl, ...
        initial_state);

    % Filter to steady state only
    ss_lin = t_lin >= T_start;
    ss_nl = t_nl >= T_start;

    % Tracking errors
    e_lin = y_lin(:,1) - d_all_lin(:,2);
    e_nl = y_nl(:,1) - d_all_nl(:,2);

    % === B2: Linear System Only ===
    fig_lin = create_single_system_figure(omega_i, 'Linear');
    plot_single_system(t_lin(ss_lin), y_lin(ss_lin,1), d_all_lin(ss_lin,2), ...
        e_lin(ss_lin), u_lin(ss_lin), [0 0.4470 0.7410], 'Linear');
    if exportFigures
        export_figure(fig_lin, ['figures/B2_linear_omega_' num2str(omega_i)]);
    end

    % === B3: Nonlinear System Only ===
    fig_nl = create_single_system_figure(omega_i, 'Nonlinear');
    plot_single_system(t_nl(ss_nl), y_nl(ss_nl,1), d_all_nl(ss_nl,2), ...
        e_nl(ss_nl), u_nl(ss_nl), [0.8500 0.3250 0.0980], 'Nonlinear');
    if exportFigures
        export_figure(fig_nl, ['figures/B3_nonlinear_omega_' num2str(omega_i)]);
    end

    % === B4: Comparison (Both Systems) ===
    fig_both = create_comparison_figure(omega_i);
    plot_comparison(t_lin(ss_lin), y_lin(ss_lin,1), u_lin(ss_lin), e_lin(ss_lin), ...
        t_nl(ss_nl), y_nl(ss_nl,1), u_nl(ss_nl), e_nl(ss_nl), d_all_lin(ss_lin,2));
    if exportFigures
        export_figure(fig_both, ['figures/B4_comparison_omega_' num2str(omega_i)]);
    end
end


%% functions

function [t, x, u, y, d_all] = run_simulation(T, omega, alpha, ...
    period, dynamics, K, L_ctrl, initial_state)

    % Actual disturbance hitting the plant (square wave)
    d1_actual = @(t) 0.5 * square(2*pi*t/period);
    d_exo = @(t) [d1_actual(t); alpha*sin(omega*t); alpha*cos(omega*t)];

    opts = odeset('Events', @angle_event);
    
    % Control uses d_exo, plant receives d1_actual
    [t, x] = ode45(@(t, x) dynamics(t, x, K*x + L_ctrl*d_exo(t), d1_actual(t)), ...
        [0 T], initial_state);

    y = x(:, [1 3]);
    d_all = cell2mat(arrayfun(@(ti) d_exo(ti)', t, 'UniformOutput', false));
    u = (K * x' + L_ctrl * d_all')';
end

function [value, isterminal, direction] = angle_event(~, x)
    value = abs(x(3)) - pi/1.2;
    isterminal = 1;
    direction = 0;
end

function fig = create_single_system_figure(omega, system_name)
    fig = figure('Position', [100 100 800 700], 'Visible', 'on');
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle([system_name ' System: \omega = ' num2str(omega) ' rad/s'], 'FontSize', 16);
end

function plot_single_system(t, s, ref, e, u, color, name)
    % Cart position
    nexttile; hold on;
    plot(t, s, 'LineWidth', 2, 'Color', color, 'DisplayName', name);
    plot(t, ref, 'LineWidth', 1, 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Reference');
    ylabel('s (m)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);

    % Tracking error
    nexttile; hold on;
    plot(t, e, 'LineWidth', 2, 'Color', color, 'DisplayName', name);
    ylabel('e(t) = s(t) - d_2(t)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);

    % Control input
    nexttile; hold on;
    plot(t, u, 'LineWidth', 2, 'Color', color, 'DisplayName', name);
    ylabel('u (N)', 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);
end

function fig = create_comparison_figure(omega)
    fig = figure('Position', [100 100 800 700], 'Visible', 'on');
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(['Linear vs Nonlinear Comparison: \omega = ' num2str(omega) ' rad/s'], 'FontSize', 16);
end

function plot_comparison(t_lin, s_lin, u_lin, e_lin, t_nl, s_nl, u_nl, e_nl, ref)
    COLOR_LIN = [0 0.4470 0.7410];
    COLOR_NL = [0.8500 0.3250 0.0980];
    COLOR_REF = [0.5 0.5 0.5];

    % Cart position
    nexttile; hold on;
    plot(t_lin, s_lin, 'LineWidth', 2, 'Color', COLOR_LIN, 'DisplayName', 'Linear');
    plot(t_nl, s_nl, 'LineWidth', 2, 'Color', COLOR_NL, 'DisplayName', 'Nonlinear');
    plot(t_lin, ref, 'LineWidth', 1, 'LineStyle', '--', 'Color', COLOR_REF, 'DisplayName', 'Reference');
    ylabel('s (m)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);

    % Tracking error
    nexttile; hold on;
    plot(t_lin, e_lin, 'LineWidth', 2, 'Color', COLOR_LIN, 'DisplayName', 'Linear');
    plot(t_nl, e_nl, 'LineWidth', 2, 'Color', COLOR_NL, 'DisplayName', 'Nonlinear');
    ylabel('e(t) = s(t) - d_2(t)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);

    % Control input
    nexttile; hold on;
    plot(t_lin, u_lin, 'LineWidth', 2, 'Color', COLOR_LIN, 'DisplayName', 'Linear');
    plot(t_nl, u_nl, 'LineWidth', 2, 'Color', COLOR_NL, 'DisplayName', 'Nonlinear');
    ylabel('u (N)', 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    grid on; set(gca, 'FontSize', 12);
end

function export_figure(fig, base_name)
    exportgraphics(fig, [base_name '.png'], 'Resolution', 300);
    exportgraphics(fig, [base_name '.pdf'], 'ContentType', 'vector');
end

function [K, L_ctrl, Pi, Gamma] = compute_control_gains(omega, A, B, P, Ce, Qe)

    n = size(A, 1);  % State dimension
    r = 3;           % Exosystem dimension

    % Exosystem matrix for constant disturbance + sinusoidal reference
    S = [0 0 0;
         0 0 omega;
         0 -omega 0];

    % Solve FBI equations: Pi*S = A*Pi + B*Gamma + P
    %                      0 = Ce*Pi + Qe
    A_fbi = [kron(S', eye(n)) - kron(eye(r), A), -kron(eye(r), B);
             kron(eye(r), Ce), zeros(r, r)];
    b_fbi = [P(:); -Qe'];

    X_fbi = A_fbi \ b_fbi;

    Pi = reshape(X_fbi(1:n*r), n, r);
    Gamma = reshape(X_fbi(n*r+1:end), 1, r);
  
  
    desired_poles = [-1,-2,-3,-4];
    K = -place(A, B, desired_poles);

    % Feedforward gain
    L_ctrl = Gamma - K * Pi;
end