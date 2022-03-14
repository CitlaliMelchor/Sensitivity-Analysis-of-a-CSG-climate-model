function [t,x] = ode_rk4(fh,t,x0)
% ode_rk4: Runge-Kutta fouth order integration for ordinary differential equation(s) (ODE)
%
%         [t,x] = ode_rk4(fh,t,x0)
%
%         Inputs:
%         fh : function handle to function that holds the derivative
%              e.g.: @(t,x)ftank(t,x,u,p)
%              if [dxdt] = ftank(t,x,u,p)
%         t  : time span: [t_start:Dt:t_end]
%         x0 : initial value state(s)
%
%         Outputs:
%         t  : time vector
%         x  : state vector(s)
%
% Programmer: Jan-Eise Vuist
% Date: 2016

%% Input checking (do not change anything here)
x0 = x0.';              % make sure that x0 is a column vector
t = t.';                % make sure that t is a column vector

if ~isnumeric(t)
   error('"t" should be a vector time points for integration.');
end

if ~isnumeric(x0)
   error('"x0" should be a vector of initial conditions.');
end

delta_t = diff(t);      % time steps, should be >0
if any(sign(delta_t(1))*delta_t <= 0)
   error('Entries of time vector "t" are not correct.');
end  

                        % for simplicity: assume time steps should be equidistant
if ~all(diff(delta_t-delta_t(1))<1e-13)
   error('Entries of time vector "t" are not equidistant.');
end
delta_t = delta_t(1);

try                     % test call to function with function handle
   f0 = feval(fh,t(1),x0);
catch err
   msg = ['Unable to evaluate the ODEFUN at t0,x0. ', err];
   error(msg);
end

if ~isequal(size(x0),size(f0))
   error('Inconsistent sizes of x0 and f(t0,x0).');
end

% initialize x (for speed)
neq = length(x0);
N = length(t);
x = zeros(neq,N);

%% Numerical integration (EDIT HERE)
x(:,1) = x0;            % initial value x(t=0)

% compute x for every timestep t(i)
% with time step defined as: delta_t
for i = 1:N-1
                        % compute dx/dt (or k) with function defined in fh
   k1 = feval(fh, t(i),           x(:,i));
   k2 = feval(fh, t(i)+delta_t/2, x(:,i)+k1*delta_t/2);
   k3 = feval(fh, t(i)+delta_t/2, x(:,i)+k2*delta_t/2);
   k4 = feval(fh, t(i)+delta_t,   x(:,i)+k3*delta_t);

   x(:,i+1) = x(:,i) + delta_t*(k1+2*k2+2*k3+k4)/6;
end

%% 
x = x.';                % make sure that x is a column vector

end