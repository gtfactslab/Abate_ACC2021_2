
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Conference: American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/18/2021
% Description:  This script generates Figures 2b and 2c.
%               Conservatism in the approximation of reachable sets is 
%               reduced, via the application of Theorem 1 (10 times) with
%               different shape matrixes.  The shape matrixes used are
%               stored in the data set fig2bc_T_set.mat, which is provided.

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

dt = .002;   % Timestep for simulation
th = 1;      % Simulation time-horizon



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('fig2bc_T_set.mat')

% setup figure 1
figure(1); clf;
hold on; grid on;
% plot initial paralellogram
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', 'LineWidth', 1.25);
scatter(X0(1, :), X0(2, :), 'k', 'filled');
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-.25, 3.75, -.5, 3]);

% setup figure 2
figure(2); clf;
hold on; grid on;
% plot initial paralellogram
xlabel('Transform','Interpreter','latex')
ylabel('Area','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([1, 10, 0, 3.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
T_size = size(0:dt:th, 2);


T_now = reshape(T_holder(:, end), [2, 2]);    
X0_prime = [ min((inv(T_now)*X0_Boundary)')', ...
                 max((inv(T_now)*X0_Boundary)')'];
X0_prime_Boundary = T_now*makeRectangle(X0_prime);
x_emb = X0_prime(:);
for t = 1:T_size
    x_emb = x_emb + dt * [dT(x_emb(1:2), W(:, 1), x_emb(3:4), W(:, 2), T_now); ...
                          dT(x_emb(3:4), W(:, 2), x_emb(1:2), W(:, 1), T_now)];
end
[Q1, Q2] =  meshgrid(x_emb(1):.01:x_emb(3), x_emb(2):.01:x_emb(4));

collect = T_now*[Q1(:), Q2(:)]';

k = boundary(collect(1, :)', collect(2, :)');
g = patch(collect(1, k), collect(2, k), 'g', ...
                                        'FaceAlpha', .1, ...
                                        'LineWidth', 1.25);
drawnow

T_set_new = T_holder(:, end);
Area_holder = [];
area_last = 1000;

for i = 1:size(T_holder, 2)
    i/size(T_holder, 2)
    T_now = reshape(T_holder(:, end - i + 1), [2, 2]);
    
    X0_prime = [ min((inv(T_now)*X0_Boundary)')', ...
                 max((inv(T_now)*X0_Boundary)')'];
    X0_prime_Boundary = T_now*makeRectangle(X0_prime);
    
    % simulate embedding system for
    x_emb = X0_prime(:);
    for t = 1:T_size
        x_emb = x_emb + dt * [dT(x_emb(1:2), W(:, 1), x_emb(3:4), W(:, 2), T_now); ...
                              dT(x_emb(3:4), W(:, 2), x_emb(1:2), W(:, 1), T_now)];
    end
    
    % checks to make sure embedding system didn't diverge 
    if norm(x_emb) <= 10000
        Approx = makeRectangle([x_emb(1:2), x_emb(3:4)]);
        
        holder = [];
        for j = 1:size(collect, 2)
            if prod(x_emb(1:2) <= inv(T_now)*collect(:, j)) == 1 && ...
               prod(x_emb(3:4) >= inv(T_now)*collect(:, j)) == 1
                holder = [holder, collect(:, j)];
            end
        end
        collect = holder;
        k = boundary(collect(1, :)', collect(2, :)', 0);
        
        figure(1);
        delete(g);
        g = patch(collect(1, k), collect(2, k), 'g', ...
                                            'FaceAlpha', .1, ...
                                            'LineWidth', 1.25);
        
        figure(2); 
        Area_holder = [Area_holder, polyarea(collect(1,k)',  collect(2, k)')];
        if i ~=1
            delete(q)
        end
        q = plot(1:i, Area_holder, 'Color', 'red', 'LineWidth', 2);
        drawnow;
    end
end


% Compute R(1, X0)
Phi = X0_Boundary;
holder = Phi;
for t = 1:T_size
    t
    % get next FORWARD TIME reachable set
    holder2 = [];
    for i = 1:size(holder, 2)
            x = holder(:, i);
            for w = W(1):.25:W(2)
                x_next = x + dt*dxdt(x, w);
                holder2 = [holder2, x_next];
            end
    end
    k = boundary(holder2(1, :)',holder2(2, :)',0.02);
    holder = holder2(:, k);
    Phi = [Phi, holder];
end
x_T = holder; % reachable set


figure(1)
patch(x_T(1, :), x_T(2, :), 'w', 'FaceAlpha', 1, 'LineWidth', 1.25);
patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25);
axis([-.25, 3.75, -.5, 3]);
xticks([0 1 2 3])
yticks([0 1 2 3])
drawnow


Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

%matlab2tikz('F2b.tikz', 'width', '6cm', 'height', '4cm')

figure(2)
area_true_reach = polyarea(x_T(1, :)', x_T(2, :)')
plot([1, 10], area_true_reach*[1, 1], 'b--', 'LineWidth', 2);

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';


%matlab2tikz('F2c.tikz', 'width', '6cm', 'height', '4cm')


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

function out = dydt(y, w, T)
    out = inv(T)*dxdt(T*y, w);
end

function out = dT(y, w, yhat, what, T)
    if y(1) <= yhat(1) && y(2) <= yhat(2) && w <= what
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * dydt([y(1); x], 0, T);
        fun_y2 = @(x) [0, 1] * dydt([x; y(2)], 0, T);
        fun_w1 = @(v) [1, 0] * dydt([0; 0], v, T);
        fun_w2 = @(v) [0, 1] * dydt([0; 0], v, T);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w, what, options);
        out_w2 = fminbnd(fun_w2, w, what, options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
        out = out - dydt([0; 0], 0, T);
        
    elseif yhat(1) <= y(1) && yhat(2) <= y(2) && what <= w
        fun_y1 = @(x) -[1, 0]* dydt([y(1); x], 0, T);
        fun_y2 = @(x) -[0, 1]* dydt([x; y(2)], 0, T);
        fun_w1 = @(v) -[1, 0]* dydt([0; 0], v, T);
        fun_w2 = @(v) -[0, 1]* dydt([0; 0], v, T);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what, w, options);
        out_w2 = fminbnd(fun_w2, what, w, options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        out = out - dydt([0; 0], 0, T);
    else
        out = [inf; inf];
    end
end
