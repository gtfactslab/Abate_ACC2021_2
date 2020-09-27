
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Submitted to American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/27/2020
% Description:  This script generates Figure 1.
%               Parallelogram approximations of reachable sets are computed
%               for a given system by applying Theorem 1.

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
% Initial Set
Y0 = [ 0, .25; ....
       -.25, 0];
% Disturbance Bound
W = [0, .25];

global T
T = [1,-2; ...
     1, 1]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Timestep for simulation
th = 1;      % Simulation time-horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y0_Boundary = makeRectangle(Y0);
Phi0 = Y0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
nPhi = Phi0;

holder = Phi;
nholder = Phi;

T_size = size(0:dt:th, 2);
yu = zeros(2, T_size + 1);  yu(:, 1) = Y0(:, 1);
yo = zeros(2, T_size + 1);  yo(:, 1) = Y0(:, 2);

nyu = zeros(2, T_size + 1);  nyu(:, 1) = Y0(:, 1);
nyo = zeros(2, T_size + 1);  nyo(:, 1) = Y0(:, 2);

%Compute R(1 ,X0)
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
    % get next BACKWARD TIME reachable set
    nholder2 = [];
    for i = 1:size(nholder, 2)
            ny = nholder(:, i);
            for w = W(1):.25:W(2)
                ny_next = ny - dt*dydt(ny, w);
                nholder2 = [nholder2, ny_next];
            end
    end
    k = boundary(nholder2(1, :)', nholder2(2, :)',0.02);
    nholder = nholder2(:, k);
    nPhi = [nPhi, nholder];
    
    % propegate corners with decomp function
    % Compute Approximation of R(1 ,X0)
    yu(:, t + 1) = yu(:, t) + dt*d(yu(:, t), W(1), yo(:, t), W(2));
    yo(:, t + 1) = yo(:, t) + dt*d(yo(:, t), W(2), yu(:, t), W(1));
    % Compute Approximation of S(1 ,X0)
    nyu(:, t + 1) = nyu(:, t) + dt*nd(nyu(:, t), W(1), nyo(:, t), W(2));
    nyo(:, t + 1) = nyo(:, t) + dt*nd(nyo(:, t), W(2), nyu(:, t), W(1));
end
y_T = holder;
yu_T = yu(:, t + 1); yo_T = yo(:, t + 1);
ny_T = nholder;
nyu_T = nyu(:, t + 1); nyo_T = nyo(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on;

% Plot Initial Set X0
X0_Boundary = T*Y0_Boundary;
X0_corn = T*Y0;
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', 'LineWidth', 1.25);
scatter(X0_corn(1, :), X0_corn(2, :), 'k', 'filled');

% Plot forward time reachable set and MM approximation
RFY = makeRectangle([yu_T, yo_T]);
RFX = T*RFY;
x_T = T*y_T; 
xu_T = T*yu_T;
xo_T = T*yo_T;
patch(RFX(1, :), RFX(2, :), 'g', 'FaceAlpha', .1, 'LineWidth', 1.25);
patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25);
scatter([xu_T(1, 1), xo_T(1, 1)], [xu_T(2, 1), xo_T(2, 1)], 'k', 'filled');

% Plot backward time reachable set and MM approximation
SFY = makeRectangle([nyu_T, nyo_T]);
SFX = T*SFY;
nx_T = T*ny_T; 
nxu_T = T*nyu_T;
nxo_T = T*nyo_T;
patch(SFX(1, :), SFX(2, :), 'b', 'FaceAlpha', .05, 'LineWidth', 1.25);
patch(nx_T(1, :), nx_T(2, :), 'b', 'FaceAlpha', .6, 'LineWidth', 1.25);
scatter(nxu_T(1, 1), nxu_T(2, 1), 'k', 'filled');
scatter(nxo_T(1, 1), nxo_T(2, 1), 'k', 'filled');

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 3, -3, 3]);
xticks([-1 0 1 2 3])
yticks([-3 -1.5 0 1.5 3])

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

% Generate tikz figure
% matlab2tikz('F1.tikz', 'width', '6cm', 'height', '4cm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dynamics
function out = dxdt(x, w)
    out = [x(1, 1)*x(2, 1) + w; ...
           x(1, 1) + 1]; 
end

% Transformed dynamics
function out = dydt(y, w)
    global T
    out = inv(T)*dxdt(T*y, w);
end

% Tight decomposition fucntion for transformed dynamics
function out = d(y, w, yhat, what)
    global T
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

% Tight decomposition fucntion for backward time system
function out = nd(y, w, yhat, what)
    global T
    if y(1) <= yhat(1) && y(2) <= yhat(2) && w <= what
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) -[1, 0] * dydt([y(1); x], 0);
        fun_y2 = @(x) -[0, 1] * dydt([x; y(2)], 0);
        fun_w1 = @(v) -[1, 0] * dydt([0; 0], v);
        fun_w2 = @(v) -[0, 1] * dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w, what, options);
        out_w2 = fminbnd(fun_w2, w, what, options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        out = out + dydt([0; 0], 0);
    elseif yhat(1) <= y(1) && yhat(2) <= y(2) && what <= w
        fun_y1 = @(x) [1, 0]* dydt([y(1); x], 0);
        fun_y2 = @(x) [0, 1]* dydt([x; y(2)], 0);
        fun_w1 = @(v) [1, 0]* dydt([0; 0], v);
        fun_w2 = @(v) [0, 1]* dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what, w, options);
        out_w2 = fminbnd(fun_w2, what, w, options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        out = out + dydt([0; 0], 0);
    else
        error
    end
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end