
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Conference: American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 3/18/2021
% Description:  This script generates Figure 3.
%               A linear transformation is found to transform a given
%               system to a monotone system.  This allows for the
%               computation of the tightest parallelogram containing the
%               reachable set of the system via Theorem 1.  A second
%               approximation is also computed via a different
%               transformation, and this allows for reduced conservatism in
%               the approximation.


clc; clear all;

% simulation parameters
dt = .01;
Sim_Time = 1;

global T W
T1 = [ 1, 1 ;...
      0, 1 ];
  
T = [  1, 4 ;...
      -1, 1 ];
  
W = [-1, 1]%1];

syms x1 x2 w
F([x1; x2], w)
dydt([x1; x2], w)

x0 = [1; 0];
x_ET = [inv(T1)*x0; inv(T1)*x0];
x_ETT = [inv(T)*x0; inv(T)*x0];

holder = x0;
for t = 0:dt:Sim_Time
    holder2 = [];
    for j = 1:1:size(holder, 2)
        xnow = holder(:, j);
        for w = W(1):(W(2)-W(1))/4:W(2)
            xnext = xnow + dt*F(xnow, w);
            holder2 = [holder2, xnext];
        end
    end
    if t > 2*dt
        k = boundary(holder2(1, :)', holder2(2, :)');
        holder = holder2(:, k);
    else
        holder = holder2;
    end
    
    x_ET = x_ET + dt*ET(x_ET(1:2, 1), x_ET(3:4, 1));
    x_ETT = x_ETT + dt*ETT(x_ETT(1:2, 1), x_ETT(3:4, 1));
end
ReachT = holder;
clear k holder xnow xnext

%
figure(1); clf;
hold on; grid on;

ET_ApproxT = makeRectangle([x_ET(1:2), x_ET(3:4)]);
ET_Approx = T1*ET_ApproxT;
patch(ET_Approx(1, :), ET_Approx(2, :), 'r', 'FaceAlpha', .1, 'LineWidth', 1.25)

ETT_ApproxT = makeRectangle([x_ETT(1:2), x_ETT(3:4)]);
ETT_Approx = T*ETT_ApproxT;
patch(ETT_Approx(1, :), ETT_Approx(2, :), 'b', 'FaceAlpha', .1, 'LineWidth', 1.25)

%plot tru reachable set, from exhaustive simulation
disc = 3;
patch(ReachT(1, 1:disc:end), ReachT(2, 1:disc:end), 'w', 'FaceAlpha', 1, 'LineWidth', 1.25)
patch(ReachT(1, 1:disc:end), ReachT(2, 1:disc:end), 'g', 'FaceAlpha', 1, 'LineWidth', 1.25)

scatter(x0(1), x0(2), 80, 'r', 'filled')

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 6, -.25, 2.25]);
xticks([0:2:6])
yticks([0, 1, 2])


Leg = legend();
set(Leg,'visible','off')
grid on;
ax.Layer = 'top';

drawnow

% matlab2tikz('F3.tikz', 'width', '6cm', 'height', '4cm')

function out = F(x, w)
    out = [ x(1) - x(2) + x(2)^3 + w; ...
            x(1) - x(2)];
end

function out = dydt(x, w)
    global T
    out = inv(T)*F(T*x, w);
end

function out = dT(x, w, xhat, what)
    out = [ x(2)^3 + w; ...
            x(1)];
end

function out = ET(x, xhat)
    global W
    out(1:2, 1) = dT(x, W(1), xhat, W(2));
    out(3:4, 1) = dT(xhat, W(2), x, W(1));
end
function out = ETT(x, xhat)
    global W
    out(1:2, 1) = dTT(x, W(1), xhat, W(2));
    out(3:4, 1) = dTT(xhat, W(2), x, W(1));
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dTT(y, w, yhat, what)
    global T
    if y(1) <= yhat(1) && y(2) <= yhat(2) && w <= what
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * dydt([y(1); x], 0);
        fun_y2 = @(x) [0, 1] * dydt([x; y(2)], 0);
        fun_w1 = @(v) [1, 0] * dydt([0; 0], v);
        fun_w2 = @(v) [0, 1] * dydt([0; 0], v);
        
        options = optimset('TolX', 1e-8);
        [~, out_y1, ex(1, 1)] = fminbnd(fun_y1, y(2), yhat(2), options);
        [~, out_y2, ex(2, 1)] = fminbnd(fun_y2, y(1), yhat(1), options);
        [~, out_w1, ex(3, 1)] = fminbnd(fun_w1, w, what, options);
        [~, out_w2, ex(4, 1)] = fminbnd(fun_w2, w, what, options);
        
        if isequal(ex, ones(4, 1))
            out(1, 1) = out_y1 + out_w1;
            out(2, 1) = out_y2 + out_w2;
        else
            out = [inf; inf];
        end
        
    elseif yhat(1) <= y(1) && yhat(2) <= y(2) && what <= w
        fun_y1 = @(x) -[1, 0]* dydt([y(1); x], 0);
        fun_y2 = @(x) -[0, 1]* dydt([x; y(2)], 0);
        fun_w1 = @(v) -[1, 0]* dydt([0; 0], v);
        fun_w2 = @(v) -[0, 1]* dydt([0; 0], v);
        
        options = optimset('TolX', 1e-8);
        [~, out_y1, ex(1, 1)] = fminbnd(fun_y1, yhat(2), y(2), options);
        [~, out_y2, ex(2, 1)] = fminbnd(fun_y2, yhat(1), y(1), options);
        [~, out_w1, ex(3, 1)] = fminbnd(fun_w1, what, w, options);
        [~, out_w2, ex(4, 1)] = fminbnd(fun_w2, what, w, options);
        
        if isequal(ex, ones(4, 1))
            out(1, 1) = - out_y1 - out_w1;
            out(2, 1) = - out_y2 - out_w2;
        else
            out = [-inf; -inf];
        end
    else
        out = [inf; inf];
    end
end





