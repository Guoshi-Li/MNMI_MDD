
function [tt, xx] = Model_NEURAL(N, dt, ST, SC, gC, WEE, WEI, WIE, PI, FLAG_Noise) 

rng(66,'twister');

% seed = 6;
% rng(seed);

global NR; 
global A;
global C;
global P;

global tau_e;
global tau_i;
global u;
global sigma; 
global Wee;
global Wei;
global Wie;

global u_noise;
global sigma_noise; 

NR = N; % Number of brain regions

Tend= ST;  % Sec
DT = dt;   % Integration step (sec)

% Local parameters
Wee = diag(WEE);    % 
Wei = diag(WEI);    % 
Wie = diag(WIE);    % 

u = 1;      % 1
sigma = 0.25;  % 0.25


tau_e = 20*1e-3; % sec
tau_i = 20*1e-3; % sec

u_noise = 0;
sigma_noise = 0.3;   % 

A = SC;
C = gC;  % global scaling factor 

% Spontanoues input 
P = PI*ones(N,1);

NE = 2*NR;  % Number of equations

X0 = zeros(1, NE);
X = X0;

tt = [];
xx = [];

t = 0;
count = 1;

while (t<=Tend) 
    tt = [tt t];
    xx = [xx; X]; 
    
    if (FLAG_Noise == 0)
        Noise_E = 0 + 0*randn(NR, 1);
    else
        Noise_E = u_noise + sigma_noise*randn(NR, 1);
        Noise_I = u_noise + sigma_noise*randn(NR, 1);
    end

    X = RK(NE, DT, t,X, Noise_E); 
    t = t+DT;
    count = count + 1;
end

end

  
%=======================================================================
%       Fourth-order Runge-Kutta [RK(4)] method (Neural Model)
%=======================================================================

function y = RK(n, h, t, y, noise)
  h_half = h/2; 
  s = zeros(1,n);
  yk = zeros(1,n);
    
  %========================
  f = WC(t, y, noise);
  s = f;
  
  %========================
    
  tk = t + h_half;
  yk = y + h_half*f;
  f = WC(tk, yk, noise);
  s = s + 2*f;
  
  %========================
   
  yk = y + h_half*f;
  f = WC(tk, yk, noise);
  s = s + 2*f;
  
 %========================
  
  tk = t + h;
  yk = y + h*f;  
  f = WC(tk, yk, noise);
  
%========================

  y = y + h/6*(s+f);
  
  
  
end




%==============================================
%      Neural model dynamics
%==============================================

function dS=WC(t, XS, noise)

global NR; 
global C;
global P;
global A;
global tau_e;
global tau_i;
global Wee;
global Wei;
global Wie;
% global u_noise;
% global sigma_noise; 


E=XS(1:NR)';
I=XS(NR+1:end)';


dEdt = 1/tau_e*( -E + FS( C.*A*E +  Wee*E - Wie*I + P  + noise )   );

dIdt = 1/tau_i*( -I + FS(    Wei*E  + noise )  );   

% dI1dt = 1/tau_i*( -I1 + 0  );   

dS=[dEdt' dIdt'];
   
end





function sf=FS(x)
 
global u;
global sigma; 

sf = 1./( 1 + exp(-(x-u)/sigma) );

end







