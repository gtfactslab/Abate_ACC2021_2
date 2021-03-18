
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Conference: American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/18/2021
% Description:  This script generates Figure 4a.
%               A hexagnol set of intial conditions is partitioned into
%               parallelograms, and this allows for the approximation of
%               the systems reachable set via 3 applications of Theorem 1.

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = .002;   % Timestep for simulation
th = 1;      % Simulation time-horizon
T_size = size(0:dt:th, 2);

global W
% Disturbance Bound
W = [0, 1/2];

colors = eye(3);

figure(1); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 7, -1, 7])
xticks([-1:2:7])
yticks([-1:2:7])
drawnow

global T
X0_Boundary = [];
for i = 1:3
    T = [-cos((i-1)*2*pi/3), cos(i*2*pi/3); ...
         -sin((i-1)*2*pi/3), sin(i*2*pi/3)]
      
    % Initial Set
    Y0 = [-1, 0; ....
           0, 1] + inv(T)*ones(2);

    Y0_Boundary = makeRectangle(Y0);
    X0_Boundary = [X0_Boundary, T*Y0_Boundary];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);
Phi = Phi0;
holder = Phi;
    for t = 1:T_size
        t
        % get next FORWARD TIME reachable set
        holder2 = [];
        for j = 1:size(holder, 2)
                x = holder(:, j);
                for w = W(1):(W(2)-W(1)):W(2)
                    x_next = x + dt*dxdt(x, w);
                    holder2 = [holder2, x_next];
                end
        end
        k = boundary(holder2(1, :)',holder2(2, :)',0.02);
        holder = holder2(:, k);
        Phi = [Phi, holder];
    end
    x_T = holder; % reachable set 
    patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25);
    drawnow

for i = 1:3
    T = [-cos((i-1)*2*pi/3), cos(i*2*pi/3); ...
         -sin((i-1)*2*pi/3), sin(i*2*pi/3)]
      
    % Initial Set
    Y0 = [-1, 0; ....
           0, 1] + inv(T)*ones(2);
    X0 = T*Y0; % rectangle points in original space
    Y0_Boundary = makeRectangle(Y0);
    X0_Boundary = T*makeRectangle(Y0);
    
    % Plot Initial Set X0
    patch(X0_Boundary(1, :), X0_Boundary(2, :), colors(i, :), ...
                                'FaceAlpha', .6, ...
                                'LineWidth', 1.25);
    scatter(X0(1, :), X0(2, :), 'k', 'filled');
    
    T_size = size(0:dt:th, 2);
    yu = Y0(:, 1); % par
    yo = Y0(:, 2);
    for t = 1:T_size
        % Parrallotope Approximation
        yu = yu + dt*dT(yu, W(1), yo, W(2), T);
        yo = yo + dt*dT(yo, W(2), yu, W(1), T);
    end
    xu_T = T*yu;
    xo_T = T*yo;
    RFY = makeRectangle([yu, yo]);
    RFX = T*RFY;
    % Plot R(1, X0_prime) and MM approximation
    patch(RFX(1, :), RFX(2, :), colors(i, :), 'FaceAlpha', .08, 'LineWidth', 1.25);
    scatter([xu_T(1, 1), xo_T(1, 1)], [xu_T(2, 1), xo_T(2, 1)], 'k', 'filled');
    drawnow
end





xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

%xticks([0 1 2 3])
%yticks([0 1 2 3])

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

%matlab2tikz('F4a.tikz', 'width', '6cm', 'height', '4cm')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dxdt(x, w)
    out = [x(2) + sin(x(2)) + w; ...
           x(1) + cos(x(1)) + 1]; 
end

function out = dydt(y, w)
    global T
    out = inv(T)*dxdt(T*y, w);
end

function out = dT(y, w, yhat, what, T)
    if y(1) <= yhat(1) && y(2) <= yhat(2) && w <= what
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * dydt([y(1); x], 0);
        fun_y2 = @(x) [0, 1] * dydt([x; y(2)], 0);
        fun_w1 = @(v) [1, 0] * dydt([0; 0], v);
        fun_w2 = @(v) [0, 1] * dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w, what, options);
        out_w2 = fminbnd(fun_w2, w, what, options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
        out = out - dydt([0; 0], 0);
        
    elseif yhat(1) <= y(1) && yhat(2) <= y(2) && what <= w
        fun_y1 = @(x) -[1, 0]* dydt([y(1); x], 0);
        fun_y2 = @(x) -[0, 1]* dydt([x; y(2)], 0);
        fun_w1 = @(v) -[1, 0]* dydt([0; 0], v);
        fun_w2 = @(v) -[0, 1]* dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what, w, options);
        out_w2 = fminbnd(fun_w2, what, w, options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        out = out - dydt([0; 0], 0);
    else
        error
    end
end
