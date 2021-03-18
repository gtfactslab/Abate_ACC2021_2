
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Conference: American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/18/2021
% Description:  This script generates Figure 2a.
%               Parallelogram approximations of forward time reachable sets
%               are computed for a given system by applying Theorem 1 and
%               Proposition 1.

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global T
T = [1,-2; ...
     1, 1]; 
% Initial Set
Y0 = [ 0, .25; ....
       -.25, 0];

Y0_Boundary = makeRectangle(Y0);
X0_Boundary = T*Y0_Boundary; % Boundary of Inital set in x coordinates
X0 = T*Y0;
   
global W
% Disturbance Bound
W = [0, .25];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Timestep for simulation
th = 1;      % Simulation time-horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute R(1, X0) and MM approximation
Phi0 = Y0_Boundary;
Phi_size = size(Phi0, 2);
Phi = Phi0;

X0_prime = [min(X0_Boundary')', max(X0_Boundary')'];
X0_prime_Boundary = makeRectangle(X0_prime);

T_size = size(0:dt:th, 2);
yu = zeros(2, T_size + 1);  yu(:, 1) = Y0(:, 1); % par
yo = zeros(2, T_size + 1);  yo(:, 1) = Y0(:, 2);
xu = zeros(2, T_size + 1);  xu(:, 1) = X0_prime(:, 1); % rect
xo = zeros(2, T_size + 1);  xo(:, 1) = X0_prime(:, 2);

holder = Phi;
for t = 1:T_size
    t
    % get next FORWARD TIME reachable set
    holder2 = [];
    for i = 1:size(holder, 2)
            y = holder(:, i);
            for w = W(1):.25:W(2)
                y_next = y + dt*dydt(y, w);
                holder2 = [holder2, y_next];
            end
    end
    k = boundary(holder2(1, :)',holder2(2, :)',0.02);
    holder = holder2(:, k);
    Phi = [Phi, holder];
    
    % propegate corners with decomp function
    % Rectangular Approximation
    xu(:, t + 1) = xu(:, t) + dt*d(xu(:, t), W(1), xo(:, t), W(2));
    xo(:, t + 1) = xo(:, t) + dt*d(xo(:, t), W(2), xu(:, t), W(1));
    % Parrallotope Approximation
    yu(:, t + 1) = yu(:, t) + dt*dT(yu(:, t), W(1), yo(:, t), W(2));
    yo(:, t + 1) = yo(:, t) + dt*dT(yo(:, t), W(2), yu(:, t), W(1));
end
y_T = holder; % reachable set

yu_T = yu(:, t + 1);
yo_T = yo(:, t + 1);
RFY = makeRectangle([yu_T, yo_T]);
RFX = T*RFY;
x_T = T*y_T; 
xu_T = T*yu_T;
xo_T = T*yo_T;

xu_T_p = xu(:, t + 1);
xo_T_p = xo(:, t + 1);
RFX_p = makeRectangle([xu_T_p, xo_T_p]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on;

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', 'LineWidth', 1.25);
scatter(X0(1, :), X0(2, :), 'k', 'filled');
patch(X0_prime_Boundary(1, :), X0_prime_Boundary(2, :), 'r', ...
                            'LineWidth', 1.25, ...
                            'FaceAlpha', .2);
scatter(X0_prime(1, :), X0_prime(2, :), 'k', 'filled');



% Plot R(1, X0_prime) and MM approximation
patch(RFX_p(1, :), RFX_p(2, :), 'b', 'FaceAlpha', .08, ...
                                     'LineWidth', 1.25);
scatter([xu_T_p(1, 1), xo_T_p(1, 1)], [xu_T_p(2, 1), xo_T_p(2, 1)], ...
                                     'k', 'filled');

% Plot R(1, X0) and MM approximation
patch(RFX(1, :), RFX(2, :), 'g', 'FaceAlpha', .1, 'LineWidth', 1.25);
patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25);
scatter([xu_T(1, 1), xo_T(1, 1)], [xu_T(2, 1), xo_T(2, 1)], 'k', 'filled');

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-.25, 3.75, -.5, 3]);
xticks([0 1 2 3])
yticks([0 1 2 3])

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

%matlab2tikz('F2a.tikz', 'width', '6cm', 'height', '4cm')


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
    out = [x(1, 1)*x(2, 1) + w; ...
           x(1, 1) + 1]; 
end

function out = dydt(y, w)
    global T
    out = inv(T)*dxdt(T*y, w);
end

function out = d(x, w, xhat, what)
    if x(1) >= 0
        out(1, 1) = x(1)*x(2)+w;
    else
        out(1, 1) = x(1)*xhat(2)+w;
    end
    out(2, 1) = x(1) + 1;
end

function out = dT(y, w, yhat, what)
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
