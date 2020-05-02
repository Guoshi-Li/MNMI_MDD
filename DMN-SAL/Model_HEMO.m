
% Compute the BOLD signals in all ROIs sequentially!!! 

function [tt, BOLD] = Model_HEMO(N, dt, ST, XX, FLAG_Mean_BOLD) 

rng(66,'twister');

global NR; 

% Hemodynamic model
global kappa;
global gamma;
global tau_hemo;
global alpha;
global rho;
global V0;


NR = N; % Number of brain regions

Tend= ST;  % Sec
DT = dt;   % Integration step (sec)

% Local parameters

kappa = 0.65;  % DEFAULT: 0.65; Oscillation increases with smaller value
gamma = 0.41;  % DEFAULT: 0.41; Oscillation decreases with smaller value

tau_hemo = 0.98;  % Sec
alpha = 0.32;

rho = 0.34;  % 0.34;   0.8
V0 = 0.02;

s0 = 0;
f0 = 1.0;
v0 = 1.0;
q0 = 1.0;

k1 = 7*rho;
k2 = 2;  % 2 / 0.4
k3 = 2*rho-0.2;


NE = 4;  % Number of equations

%=================================== 
%     Hemodynamic model
%=================================== 
  
for i = 1:NR
  E = XX(:,i);
  I = XX(:,i+NR);
  EI = (2/3)*E + (1/3)*I;
  
  if (FLAG_Mean_BOLD==1)
    XM = mean(EI);
    EI = EI - XM;
  end
  
  X0 = [s0 f0 v0 q0];
  X = X0;
  
  tt = [];
  yy = [];

  t = 0;
  count = 1;

  while (t<=Tend) 
    
    tt = [tt t];
    yy = [yy; X]; 
    Z = EI(count);
    
    X = RK2(NE, DT, t,X, Z); 
    t = t+DT;
    count = count + 1;
  end

  S  = yy(:,1);
  F  = yy(:,2);
  V  = yy(:,3);
  Q  = yy(:,4);

  BOLD(:,i) = V0*( k1*(1-Q) + k2*(1-Q./V) + k3*(1-V) );

end




end




%=======================================================================
%       Fourth-order Runge-Kutta [RK(4)] method (Hemodynamic Model)
%=======================================================================

function y = RK2(n, h, t, y, z)
  h_half = h/2; 
  s = zeros(1,n);
  yk = zeros(1,n);
    
  %========================
  f = HEMO(t, y, z);
  s = f;
  
  %========================
    
  tk = t + h_half;
  yk = y + h_half*f;
  f = HEMO(tk, yk, z);
  s = s + 2*f;
  
  %========================
   
  yk = y + h_half*f;
  f = HEMO(tk, yk, z);
  s = s + 2*f;
  
 %========================
  
  tk = t + h;
  yk = y + h*f;  
  f = HEMO(tk, yk, z);
  
%========================

  y = y + h/6*(s+f);
  
  
  
end




%==============================================
%              Hemo-dynamics
%==============================================

function dS=HEMO(t,XS, Z)

global kappa;
global gamma;
global tau_hemo;
global alpha;
global rho;

s=XS(1);
f=XS(2);
v=XS(3);
q=XS(4);

dsdt = Z - kappa*(s) - gamma*(f-1);
dfdt = s;
dvdt = 1/tau_hemo*( f - v^(1/alpha) );
dqdt = 1/tau_hemo*( f*OE(f, rho)/rho - v^(1/alpha)*q/v );

dS=[dsdt dfdt dvdt dqdt];
   
end


function y=OE(flow, rho)

y = 1-(1-rho)^(1/flow);

end










