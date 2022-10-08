% -------------------------------------------------------------------------
% 
% 04.11.2021
% Yiğit Günsür ELMACIOĞLU
% 
% https://www.youtube.com/watch?v=m5sEln5bWuM&list=PLxdnSsBqCrrEx3A6W94s...
% QGClk6Q4YCg-h&index=16&ab_channel=ChristopherLum
% 
% adresindeki uçak modeli kullanılarak oluşturuldu.
% -------------------------------------------------------------------------
function x_dot = RCAM_MODEL(x,U,rho)
u = x(1) ;    % x velocity in body frame
v = x(2) ;    % y velocity in body frame
w = x(3) ;    % z velocity in body frame
p = x(4) ;    % x angular velocity
q = x(5) ;    % y angular velocity
r = x(6) ;    % z angular velocity
phi  = x(7) ; % roll angle
teta = x(8) ; % pitch angle
ksi  = x(9) ; % yaw angle

u1 = U(1) ;   %aileron
u2 = U(2) ;   %elevator
u3 = U(3) ;   %rudder
u4 = U(4) ;   %throttle 1
u5 = U(5) ;   %throttle 2

m = 120000 ;  % mass of the plane
cbar = 6.6 ;
lt = 24.8 ;   %distance btw aerodynamic center of the tail and CG
S = 260 ;     %wing planform area
St = 64 ;     %tail unit planform area

% Center of gravity position
Xcg = 0.23*cbar ;
Ycg = 0 ;
Zcg = 0.1*cbar ;

% Aerodynamic center position
Xac = 0.12*cbar ;
Yac = 0 ;
Zac = 0 ;

% Motor 1 position
Xapt1 = 0 ;
Yapt1 = -7.94 ;
Zapt1 = -1.9 ;

%  Motor 2 position
Xapt2 = 0 ;
Yapt2 = 7.94 ;
Zapt2 = -1.9 ;

g = 9.81 ;
depsda = 0.25 ;       %derivative of epsilon wrt alpha
alpha_L0 = -11.5*pi/180 ;   %zero lift angle
n = 5.5 ;
a0 = 15.212 ;
a1 = -155.2 ;
a2 = 609.2  ;
a3 = -768.5 ;
alpha_switch = 14.5*pi/180 ;

Va = sqrt( u^2 + v^2 + w^2 ) ;  % Velocity
alpha = atan2( w,u ) ;          % Angle of attack
beta = asin( v/Va ) ;           % Side-slip angle
Q = 0.5*rho*Va^2 ;              % Dynamic pressure
wbe_b = [p; q; r] ;  % angular velocity matrix in body frame
v_b = [u; v; w] ;    % velocity matrix in body frame

% Lift coefficient of wings
if ( alpha<=alpha_switch )
    CL_wb = n*(alpha - alpha_L0) ;
else
    CL_wb = a3*alpha^3 + a2*alpha^2 + a1*alpha + a0 ;
end

% Lift coefficient
epsilon = depsda*(alpha - alpha_L0) ;
alpha_t = alpha - epsilon + u2 + 1.3*q*lt/Va ;
CL_t = 3.1*(St/S)*alpha_t ;

CL = CL_wb + CL_t ;                       % resultant lift coefficient
CD = 0.13 + 0.07*(5.5*alpha + 0.654)^2 ;  % Drag coefficient
CY = -1.6*beta + 0.24*u3 ;                % Side force coefficient

FA_s = [ -CD*Q*S; CY*Q*S; -CL*Q*S ] ;     % Aerodynamic forces in wind axis

% Direction cosine matrix for wind-body frame
C_bs = [ cos(alpha) 0 -sin(alpha);
         0          1       0    ; 
         sin(alpha) 0  cos(alpha)] ;

FA_b = C_bs*FA_s ;  % Aerodynamic forces in body frame

% equation 2.33
eta11 = -1.4*beta ;
eta21 = -0.59 - 3.1*(St*lt)/(S*cbar)*(alpha-epsilon) ;
eta31 = (1 - alpha*180/15/pi)*beta ;

eta = [ eta11 ; eta21 ; eta31 ] ;

%moment change with respect to states
dCMdx = (cbar/Va)*[-11             0              5;
                    0 -4.03*(St*lt^2)/(S*cbar^2) 0;
                    1.7            0            -11.5 ] ;

%control effectiveness
dCMdu = [-0.6            0         0.22;
          0  (-3.1*(St*lt))/(S*cbar) 0;
          0              0        -0.63] ;     
          
CMac_b = eta + dCMdx*wbe_b + dCMdu*[ u1; u2; u3 ] ;
MAac_b = CMac_b*S*Q*cbar ;

rcg_b = [ Xcg; Ycg; Zcg ] ; % Center of gravity position
rac_b = [Xac; Yac; Zac ] ;  % Aerodynamic center 
MAcg_b = MAac_b + cross(FA_b, rcg_b - rac_b) ; % Resultant moment about CG

% Forces of motors
F1 = 10*u4*m*g ;
F2 = 10*u5*m*g ;
FE1_b = [F1;0;0] ;
FE2_b = [F2;0;0] ;
FE_b = FE1_b + FE2_b ;

% Moment arms of the motors
new1 = [Xcg - Xapt1 ;
        Yapt1 - Ycg ;
        Zcg - Zapt1 ] ;    
new2 = [Xcg - Xapt2 ;
        Yapt2 - Ycg ;
        Zcg - Zapt2 ] ;    

% Moments of the motors about CG
MEcg1_b = cross(new1,FE1_b) ;
MEcg2_b = cross(new2,FE2_b) ;

% Resultant moment of the motors
MEcg_b = MEcg1_b + MEcg2_b ;

% Gravity in body frame
g_b = [-g*sin(teta);
       g*cos(teta)*sin(phi);
       g*cos(teta)*cos(phi)] ;
Fg_b = m*g_b ;

% Inertia matrix of the plane in body frame , video 17.58
Ib = m*[40.07  0  -2.0923;
          0    64   0   ;
      -2.0923  0   99.92  ] ;

% Inverse of the inertia matrix,  video 17.58
invIb = 1/m*[0.0249836  0   0.000523151 ;
             0       0.015625    0      ;
             0.000523151 0   0.010019   ] ;

% Resultant force in body frame
F_b = Fg_b + FE_b + FA_b ;

% Speed changes in body frame
x13dot = 1/m*F_b - cross(wbe_b,v_b) ;

% Resultant moment about CG
Mcg_b = MAcg_b + MEcg_b ;

% Angular speed changes in body frame
x46dot = invIb*(Mcg_b - cross(wbe_b,Ib*wbe_b)) ;

H_phi = [1  sin(phi)*tan(teta)  cos(phi)*tan(teta) ;
         0       cos(phi)          -sin(phi)       ;
         0  sin(phi)/cos(teta)  cos(phi)/cos(teta) ] ;

% Euler angle changes 
x79dot = H_phi*wbe_b ;

% Direction cosine matrices
C1v = [cos(ksi) sin(ksi) 0;
       -sin(ksi) cos(ksi) 0;
       0        0        1] ;      % Yaw turn

C21 = [cos(teta) 0 -sin(teta);
       0         1      0    ;
       sin(teta) 0  cos(teta)] ;   % Pitch turn
   
Cb2 = [1         0     0     ;
       0     cos(phi) sin(phi);
       0     -sin(phi) cos(phi)] ; % Roll turn

% Yaw-pitch-roll DCM   from global to body frame
Cbv = Cb2*C21*C1v ;

% from body to global frame
Cvb = Cbv' ;

% Velocity matrix in global frame
x1012dot = Cvb*v_b ;

% Total derivative matrix
x_dot = [x13dot; x46dot; x79dot; x1012dot(1); x1012dot(2); -x1012dot(3) ] ;  
end
  


