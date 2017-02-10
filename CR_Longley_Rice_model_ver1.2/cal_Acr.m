function [ A_cr,func,loc1,loc2 ] = cal_Acr( terrain_data, data_GAP, f_MHz, Tx_height, Rx_height, TGC_type, Ns, H_V )
% This program is used to calculate the Transmission Loss over Irregular
% Terrain; the data of the terrain must be known exactly.
%   terrain_data  must be array
%   data_Gap      stands for the distance between terrain data point
%                   (meter)
%   f_MHz         frequency in MHz (can be an array!)
%   Tx_height     Transmiting antenna height in meter
%   Rx_height     Reciving antenna height in meter
%   TGC_type:
%        1 stands for Poor ground
%        2 stands for Average ground
%        3 stands for Good ground (default)
%        4 stands for Sea water
%        5 stands for Fresh water
%   Ns            corresponding surface refractivity (default Ns = 301)
%   H_V           horizontal/vertical polarzition(default = 1(H))
%
%   A_cr          Reference attenuation in dB
%   func          stands for which Function is used to predict attenuation
%        1 stands for Two-Ray Optics
%        2 stands for Diffraction Attenuation
%        3 stands for Scatter Attenuation


% Precoding
if nargin < 8, H_V = 1; end
if nargin < 7, Ns = 301; end
if nargin < 6, TGC_type = 3; end

% Calculating useful parameter
global a Sigma Epsilon fun
a = 6370*(1-0.04665*exp(0.005577*Ns))^(-1); %effective earth's radius
switch TGC_type            % ¦Ò & ¦Å
    case 1
        Sigma = 0.001; Epsilon = 4;
    case 2
        Sigma = 0.005; Epsilon = 15;
    case 3
        Sigma = 0.02; Epsilon = 25;
    case 4
        Sigma = 5; Epsilon = 81;
    case 5
        Sigma = 0.01; Epsilon = 81;
end

L_of_T = size(terrain_data,2);
d_km = (L_of_T-1)*data_GAP/1000; % kilometer
h_s1 = Tx_height + terrain_data(1); % the antenna height above sea level in meters
h_s2 = Rx_height + terrain_data(L_of_T);

[ ele_E1,d_L1,h_L1,loc1 ] = cal_ele_angle( terrain_data, data_GAP, h_s1, 1);
[ ele_E2,d_L2,h_L2,loc2 ] = cal_ele_angle( terrain_data, data_GAP, h_s2, 2);
delta_h = cal_deltah ( terrain_data );
[ h_e1,h_e2 ] = cal_h_effective( Tx_height, Rx_height, h_L1, h_L2, h_s1, h_s2 );
d_L = d_L1 + d_L2; % the total distance between the antennas and their horizon.
ele_E = max(ele_E1+ele_E2, -d_L/a); % elevation angle in radians(all)


% ¦Èe>0.2 prediction unavailable
if (ele_E1>0.2 || ele_E2>0.2)
    A_cr = NaN;
else
    for i = 1:size(f_MHz,2)
        f_single = f_MHz(i);
        [ A_ed,m_d ] = cal_ademd( f_single, h_e1, h_e2, Tx_height, Rx_height,...
            ele_E, delta_h, d_L1, d_L2, H_V );
% Two-Ray Optics    
        if d_L>d_km
            A_cr(i) = com_two_ray ( f_single, d_km, h_e1, h_e2, delta_h,...
                d_L1, d_L2, A_ed, m_d, H_V );
            fun = 1;
% Diffraction Attenuation/Scatter Attenuation
        else
            A_cr(i) = com_scatter ( f_single, d_km, h_e1, h_e2, Tx_height, ...
                Rx_height, d_L1, d_L2, A_ed, m_d, ele_E, Ns, H_V );
        end
    end
end
func = fun;

end


function [ ele_Et,d_Lt,h_Lt,loc ] = cal_ele_angle( h_l, data_GAP, h_s ,direction)
% calculate elevation angle in radians
%           actual horizon distance in kilometer
%           obstacle height in meters
%       and location of the point
% all these parameters
global a
    N = size(h_l,2);
    if direction == 1
        d_l = linspace(1,N-1,N-1)*data_GAP/1000;
        ele_angle = 0.001.*(h_l(2:N) - h_s)./d_l - d_l./2/a;
        [ele_Et,loc] = max(ele_angle);
        d_Lt = d_l(loc);
        loc = loc + 1;
        h_Lt = h_l(loc);
        
    elseif direction == 2
        d_l = linspace(N-1,1,N-1)*data_GAP/1000;
        ele_angle = 0.001.*(h_l(1:N-1) - h_s)./d_l - d_l./2/a;
        [ele_Et,loc] = max(ele_angle);
        d_Lt = d_l(loc);
        h_Lt = h_l(loc);
    end
        
end

function [ delta_h ] = cal_deltah( terrain_data )
% calculate delta_h
    datas = sort( terrain_data );
    N = size( datas,2 );
    N1 = int32(N*0.1);
    N9 = int32(N*0.9);
    delta_h = datas(N9) - datas(N1);
end

% there are something wrong in effective antenna height calculation
function [ h_e1,h_e2 ] = cal_h_effective( Tx_height, Rx_height, h_L1, h_L2, h_s1, h_s2 )
% calculate effective antenna height in meters
h_e1 = max(Tx_height,h_s1-h_L1);
h_e2 = max(Rx_height,h_s2-h_L2);

end

function [ Aed,m_d ] = cal_ademd( f_MHz, h_e1, h_e2, Tx_height, Rx_height, ele_E, delta_h, d_L1, d_L2, H_V )
% calculate A_ed & m_d
global a
h_g1 = Tx_height; %m
h_g2 = Rx_height; %m
d_L = d_L1 + d_L2; %km
d_Ls1 = sqrt(0.002*a*h_e1); %m
d_Ls2 = sqrt(0.002*a*h_e2); %m
d_Ls = d_Ls1 + d_Ls2; %m
d_3 = d_L + 0.5*(a^2/f_MHz)^(1/3); %km 
if d_3 < d_Ls, d_3 = d_Ls; end
d_4 = d_3 + (a^2/f_MHz)^(1/3); %km
Lambda = 300/f_MHz; %m


if (h_g1<5)||(h_g2<5)
    C=10;
else
    C=0;
end
w_3 = (1 +0.1*(delta_h_d(delta_h,d_3)/Lambda *(sqrt((h_e1*h_e2+C)/...
    (h_g1*h_g2+C))+(a*ele_E+d_L)/d_3))^(1/2))^(-1);
w_4 = (1 +0.1*(delta_h_d(delta_h,d_4)/Lambda *(sqrt((h_e1*h_e2+C)/...
    (h_g1*h_g2+C))+(a*ele_E+d_L)/d_4))^(1/2))^(-1);
ele_E3 = ele_E + d_3/a;
ele_E4 = ele_E + d_4/a;
v_13 = 1.2915 * ele_E3 * sqrt(f_MHz*d_L1*(d_3-d_L)/(d_3-d_L2));
v_23 = 1.2915 * ele_E3 * sqrt(f_MHz*d_L2*(d_3-d_L)/(d_3-d_L1));
v_14 = 1.2915 * ele_E4 * sqrt(f_MHz*d_L1*(d_4-d_L)/(d_3-d_L2));
v_24 = 1.2915 * ele_E4 * sqrt(f_MHz*d_L2*(d_4-d_L)/(d_3-d_L1));
A_k3 = A_v(v_13)+A_v(v_23); %dB
A_k4 = A_v(v_14)+A_v(v_24); %dB

a_1 = d_L1^2/(0.002*h_e1);
a_2 = d_L2^2/(0.002*h_e2);
a_3 = (d_3-d_L)/ele_E3;
a_4 = (d_4-d_L)/ele_E4;

B_1 = 416.4*f_MHz^(1/3)*(1.607-K_hv(a_1, H_V, f_MHz));
B_2 = 416.4*f_MHz^(1/3)*(1.607-K_hv(a_2, H_V, f_MHz));
B_3 = 416.4*f_MHz^(1/3)*(1.607-K_hv(a_3, H_V, f_MHz));
B_4 = 416.4*f_MHz^(1/3)*(1.607-K_hv(a_4, H_V, f_MHz));

x_1 = B_1*a_1^(-2/3)*d_L1;
x_2 = B_2*a_2^(-2/3)*d_L2;
x_3 = B_3*a_3^(-2/3)*(d_3-d_L)+x_1;
x_4 = B_4*a_4^(-2/3)*(d_4-d_L)+x_2;

A_r3 = G_x(x_3) - F_x(x_1, K_hv(a_1, H_V, f_MHz)) - F_x(x_2, K_hv(a_2, H_V, f_MHz)) - 20; %dB
A_r4 = G_x(x_4) - F_x(x_1, K_hv(a_1, H_V, f_MHz)) - F_x(x_2, K_hv(a_2, H_V, f_MHz)) - 20; %dB

A_3 = (1-w_3)*A_k3 + w_3*A_r3;
A_4 = (1-w_4)*A_k4 + w_4*A_r4;

A_fo = min(5*log10(1+h_g1*h_g2*f_MHz*Sigma_h(delta_h_d(delta_h,d_Ls))*1e-5),15); %dB
m_d = (A_4-A_3) / (d_4-d_3);
Aed = A_fo + A_4 - m_d*d_4;
end

function [ Acr ] = com_two_ray ( f_MHz, d_km, h_e1, h_e2, delta_h, d_L1, d_L2, A_ed, m_d, H_V )
% define Gp = 10log10(g_o1*g_o2)
global a Epsilon Sigma
% confirm parameter from Diffraction Attenuaion
d_L = d_L1 + d_L2; %km
d_Ls1 = sqrt(0.002*a*h_e1); %m
d_Ls2 = sqrt(0.002*a*h_e2); %m
d_Ls = d_Ls1 + d_Ls2; %m
Lambda = 300/f_MHz; %m

% d_0 d_1
if A_ed>=0
    d_0 = min(4e-5*h_e1*h_e2*f_MHz, 0.5*d_L); %km
else
    d_01 = min(-A_ed/m_d, d_L-2);
    d_0 = max(d_01, 0.5*d_L);
end
d_1 = d_0 + 0.5*(d_L-d_0);

% 3.9 defined parameter
 x = 18000*Sigma/f_MHz;
 Phi_0 = atan((h_e1+h_e2)/1000/d_0);
 Phi_1 = atan((h_e1+h_e2)/1000/d_1);
 p_0 = sqrt((((Epsilon-cos(Phi_0)^2)^2+x^2)^(1/2)+(Epsilon-cos(Phi_0)^2))/2);
 p_1 = sqrt((((Epsilon-cos(Phi_1)^2)^2+x^2)^(1/2)+(Epsilon-cos(Phi_1)^2))/2);
 q_0 = x/2/p_0;
 q_1 = x/2/p_1;
 dSigma_h0 = Sigma_h( delta_h_d ( delta_h, d_0 ) );
 dSigma_h1 = Sigma_h( delta_h_d ( delta_h, d_1 ) );
 
 %b m c R
 if H_V == 1 %h
     b_0 = 1/(p_0^2 + q_0^2); %radians
     m_0 = 2*p_0/(p_0^2 + q_0^2);
     c_0 = atan(q_0/(p_0+sin(Phi_0))) - atan(q_0/(p_0-sin(Phi_0))); %radians
     b_1 = 1/(p_1^2 + q_1^2); %radians
     m_1 = 2*p_1/(p_1^2 + q_1^2);
     c_1 = atan(q_1/(p_1+sin(Phi_1))) - atan(q_1/(p_1-sin(Phi_1))); %radians
 else %v
     b_0 = (Epsilon^2 + x^2)/(p_0^2 + q_0^2);
     m_0 = 2*(p_1*Epsilon + q_0*x)/(p_0^2 + q_0^2);
     y10 = (x*sin(Phi_0) + q_0)/(Epsilon*sin(Phi_0) + p_0);
     y20 = (x*sin(Phi_0) - q_0)/(Epsilon*sin(Phi_0) - p_0);
     if sin(Phi_0)>=p_0
         c_0 = atan(y10) - atan(y20) + pi;
     else if p_0*sin(Phi_0)>0.5
             c_0 = atan(y10) + atan(y20);
         else 
             c_0 = atan(y10) - atan(y20);
         end
     end
     b_1 = (Epsilon^2 + x^2)/(p_1^2 + q_1^2);
     m_1 = 2*(p_1*Epsilon + q_1*x)/(p_1^2 + q_1^2);
     y11 = (x*sin(Phi_1) + q_1)/(Epsilon*sin(Phi_1) + p_1);
     y21 = (x*sin(Phi_1) - q_1)/(Epsilon*sin(Phi_1) - p_1);
     if sin(Phi_1)>=p_1
         c_1 = atan(y11) - atan(y21) + pi;
     else if p_1*sin(Phi_1)>0.5
             c_1 = atan(y11) + atan(y21);
         else 
             c_1 = atan(y11) - atan(y21);
         end
     end
 end
 R_0 = (1+b_0*sin(Phi_0)^2-m_0*sin(Phi_0))/(1+b_0*sin(Phi_0)^2+m_0*sin(Phi_0));
 R_1 = (1+b_1*sin(Phi_1)^2-m_1*sin(Phi_1))/(1+b_1*sin(Phi_1)^2+m_1*sin(Phi_1));


%Re PLD = 2*pi*deltaR/Lambda
Re__0 = R_0 * exp(-2*pi*dSigma_h0*sin(Phi_0)/Lambda); %g_r1*g_r2 =g_o1*g_o2 as default
if Re__0>0.5 && Re__0>sqrt(sin(Phi_0))
    Re_0 = Re__0;
else
    Re_0 = sqrt(sin(Phi_0));
end
Re__1 = R_1 * exp(-2*pi*dSigma_h1*sin(Phi_1)/Lambda); %g_r1*g_r2 =g_o1*g_o2 as default
if Re__1>0.5 && Re__1>sqrt(sin(Phi_1))
    Re_1 = Re__1;
else
    Re_1 = sqrt(sin(Phi_1));
end
PLD_0 = 4.1917e-5*f_MHz*h_e1*h_e2/d_0; % path length difference in radians
PLD_1 = 4.1917e-5*f_MHz*h_e1*h_e2/d_1;

%A_t A_d A_Ls
A_0t = -10*log10(1+Re_0^2 - 2*Re_0*cos(PLD_0 - c_0)); %attenuation relative to free space in dB
A_1t = -10*log10(1+Re_1^2 - 2*Re_1*cos(PLD_1 - c_1)); %Gp=10*log10(g_o1*g_o2) as default
A_0d = A_ed + m_d * d_0;
A_1d = A_ed + m_d * d_1;
A_Ls = A_ed + m_d * d_Ls;

%A_0 A_1 w_0
w_0 = (1+f_MHz*delta_h*1e-4)^-1;
A_0 = min(w_0*A_0t + (1-w_0)*A_0d, A_0d); %dB
A_1 = min(w_0*A_1t + (1-w_0)*A_1d, A_1d);

%k_1 k_2
k_2 = max(((A_Ls-A_0)*(d_1-d_0)-(A_1-A_0)*(d_Ls-d_0))/...
    ((d_1-d_0)*log10(d_Ls/d_0)-(d_Ls-d_0)*log10(d_1/d_0)), 0); %dB
k_1 = ((A_Ls-A_0)-k_2*log10(d_Ls/d_0))/(d_Ls-d_0); %dB/km
if k_1<0
    k_1 = 0;
    k_2 = (A_Ls-A_0)/log10(d_Ls/d_0);
end

% Reference Attenuation Acr
Acr = max( A_0 + k_1*(d_km-d_0) + k_2*log10(d_km/d_0), 0); %dB
end

function [ Ads ] = com_scatter ( f_MHz, d_km, h_e1, h_e2, Tx_height, Rx_height, d_L1, d_L2, A_ed, m_d, ele_E, Ns, H_V )
global a fun
% confirm parameters
d_L = d_L1 + d_L2; %km
d_5 = d_L + 200;
d_6 = d_L + 400;
ele_E5 = ele_E + d_5/a;
ele_E6 = ele_E + d_6/a;

%H_5,6 S_5,6
H_5 = min((1/h_e1 + 1/h_e2)/(ele_E5*f_MHz*abs(0.007 - 0.058*ele_E5)), 15); %dB
H_6 = min((1/h_e1 + 1/h_e2)/(ele_E6*f_MHz*abs(0.007 - 0.058*ele_E6)), 15); %dB
S_5 = H_5 + 10*log10(f_MHz*ele_E5^4)-0.1*(Ns - 301)*exp(-ele_E5*d_5/40); %dB
S_6 = H_6 + 10*log10(f_MHz*ele_E6^4)-0.1*(Ns - 301)*exp(-ele_E6*d_6/40); %dB

As_5 = as_56( S_5, ele_E5, d_5 );
As_6 = as_56( S_6, ele_E6, d_6 );
m_s = (As_6-As_5)/(d_6-d_5);
A_es = As_5 - m_s * d_5;
d_x = max((A_es - A_ed)/(m_d - m_s), d_L+0.25*(a^2/f_MHz)^(1/3)*log10(f_MHz));
%A_es m_s
if H_5 < 10
    A_5 = As_5;
    A_6 = As_6;
    m_s = (A_6-A_5)/(d_6-d_5);
    A_es = A_5 - m_s * d_5;
else
	[ A_ed_,m_d_ ] = cal_ademd( f_MHz, h_e1, h_e2, Tx_height, Rx_height,...
            ele_E, 0, d_L1, d_L2, H_V );

    d_x1 = (As_5 - m_s * d_5 - A_ed_)/(m_d_ - m_s);
    d_x2 = d_L+0.25*(a^2/f_MHz)^(1/3)*log10(f_MHz);
    d_xo = d_x1*(3-0.2*H_5) + d_x2 *(0.2*H_5-2); %km
    A_xo = A_ed_ + m_d_ * d_xo;
    A_so = as_56( S_5, ele_E5, d_xo );
    A_sx = A_xo + (As_5 - A_so);
    A_es = A_sx - m_s * d_xo;
end
%Final
if d_km<d_x
    Ads = A_ed + m_d * d_km; %dB
    fun = 2;
else
    Ads = A_es + m_s * d_km; %dB
    fun = 3;
end
if Ads<0
    Ads = 0;
end
end

function [ delta_hd ] = delta_h_d ( delta_h, d )
    delta_hd = delta_h * (1-0.8*exp(-0.02*d));
end

function [ A ] = A_v( v )
if v>=0 && v<=2.4    
    A = 6.02+9.11*v-1.27*v^2;
elseif v>2.4
    A = 12.953+20*log10(v);
end
end

function [ K ] = K_hv( aa, H_V, f_MHz )
global Epsilon Sigma
    x = 18000*Sigma/f_MHz;
    K = 0.36278*(aa*f_MHz)^(-1/3)*((Epsilon-1)^2 + x^2)^(-1/4);
if H_V == 2
    K = K*(Epsilon^2 + x^2);
end
end

function [ F ] = F_x( x, K )
if x>0 && x<=200
    if K>=0 && K<=1e-5
        F = min(40*log10(x)-117,-117);
    elseif K>=1e-5 && K<1 && x>= -450/(log10(K))^3
        F = 40*log10(x)-117;
    else
        F = 20*log10(K)+2.5e-5*x^2/K-15;
    end
elseif x>=200 && x<=2000
    w = 0.0134*x*exp(-0.005*x);
    F = w*(40*log10(x)-117)+(1-w)*(0.05751*x-10*log10(x));
else
    F = 0.05751*x-10*log10(x);
end
end

function [ G ] = G_x( x )
    G = 0.05751*x-10*log10(x);
end

function [ Sigma_h ] = Sigma_h( delta_hd )
    if delta_hd > 4
        Sigma_h = 0.78*delta_hd*exp(-0.5*(delta_hd)^(1/4)); %m
    else
        Sigma_h = 0.39*delta_hd; %m
    end
end

function [ AS_56 ] = as_56( S, ele, d )
    if ele*d<=10
        AS_56 = S + 103.4 + 0.332*ele*d - 10*log10(ele*d); %dB
    elseif ele*d<=70
        AS_56 = S + 97.1 + 0.212*ele*d - 2.5*log10(ele*d); %dB
    else
        AS_56 = S + 86.8 + 0.157*ele*d + 5*log10(ele*d); %dB
    end
end
