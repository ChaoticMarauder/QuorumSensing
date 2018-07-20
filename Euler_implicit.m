function [t,x] = Euler_implicit(fun,x0,tf,h)
% Initial value problem ODE solver
% Inputs:
% fun: Function file for dx/dt = fun(t,x)
% x0: initial condition
% tf: final time
% h: step size
% other specs
iter_max = 100; % max iterations for Newton Raphson subroutine
error_tol = 0.001; % convergence tolerance for Newton Raphson subroutine
% compute numner of steps
N = ceil(tf/h);
% initialize outputs
t = zeros(N+1,1);
x = zeros(N+1,length(x0));
x(1,:) = x0';
for i = 1:N
 t(i+1) = t(i) + h;
 % set up nonlinear equation solver
 x_current = x(i,:)';
 x_guess = x_current+h*fun(t(i),x_current);
 disp(['For t = ' num2str(t(i+1))]);
 % Newton Raphson routine
 for ii = 1:iter_max
 N1 = length(x_guess);
 f = x_current + h*fun(t(i+1),x_guess) - x_guess;
 % jacobian computation
 J = h*jacobian(fun,t(i+1),x_guess) - eye(N1);
 dx = GaussElimination_SSJ(J,-f);
 x_new = x_guess +dx;
 % check convergence
 er = norm(x_current + h*fun(t(i+1),x_new)-x_new);
 if er < error_tol
 disp(['Convergence achieved in ' num2str(ii) ' iterations']);
 break; % stops iterations
 else
 % check if max iterations have reached
 if ii == iter_max
 disp('Maximum iterations reached without satisfying convergence criteria');
 else
 % proceed to the next iteration with updated initial
 % condition
 x_guess = x_new;
 end
 end
 end
 x(i+1,:) = x_new';
end