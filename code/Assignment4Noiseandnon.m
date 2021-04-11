%% Part 3: Circuit With Noise
% In this section, a current source and capacitor was added in parallel
% with R3 to simulate noise. Therefore, a new G and C matrix was
% formulated.

G = zeros(6, 6); 

%Conductances(1/R):
G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 0.001;

%Additional Parameters:
alpha = 100;
Cval = 0.25;
L = 0.2;

Cn = 0.00001; 

G(1, 1) = 1;                                    % 1
G(2, 1) = -G1; G(2, 2) = G1 + G2;               % 2
G(3 ,2) = -1; G(3, 4) = 1;                      % iL
G(4, 3) = -1; G(4, 4) = G3; G(4, 7) = 1;        % 3
G(5, 5) = 1; G(5, 4) = -alpha*G3;               % 4
G(6, 6) = G4 + G0; G(6, 5) = -G4;               % 5
G(7, 7) = 1;                                    % In

C = zeros(7);

C(2, 1) = -Cval; C(2, 2) = Cval;
C(3, 3) = L;
C(4, 4) = Cn;

%% 
% The updated C matrix is:
C
%%

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/1000; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.001:1
    Vold = V;
    vin4(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin4(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo4(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(12);
plot(t, vin4);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo4);
title('(A) V0 vs t for guassian pulse input with noise');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

% Guassian with noise
fo = fft(vo4);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(13);
plot(fs,p);
hold on
fo = fft(vin4);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Guassian pulse Input with noise');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

%%
% Next, the value of the newly added capacitor was varied to analyze its
% effect on the bandwidth. 

C(4, 4) = 0.0001; 

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/1000; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.001:1
    Vold = V;
    vin5(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin5(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo5(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(14);
plot(t, vin5);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo5);
title('(A) V0 vs t for guassian pulse input with noise (Cn = 0.0001)');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

% Guassian with noise
fo = fft(vo5);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(15);
plot(fs,p);
hold on
fo = fft(vin5);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Guassian pulse Input with noise (Cn = 0.0001)');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

C(4, 4) = 0.001; 

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/1000; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.001:1
    Vold = V;
    vin6(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin6(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo6(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(16);
plot(t, vin4);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo6);
title('(A) V0 vs t for guassian pulse input with noise (Cn = 0.001)');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

% Guassian with noise
fo = fft(vo6);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(17);
plot(fs,p);
hold on
fo = fft(vin6);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Guassian pulse Input with noise (Cn = 0.001)');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

C(4, 4) = 0.01; 

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/1000; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.001:1
    Vold = V;
    vin7(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin7(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo7(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(18);
plot(t, vin7);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo7);
title('(A) V0 vs t for guassian pulse input with noise (Cn = 0.01)');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

% Guassian with noise
fo = fft(vo7);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(19);
plot(fs,p);
hold on
fo = fft(vin7);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Guassian pulse Input with noise (Cn = 0.01)');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

%%
% The capacitance begins to more clearly affect the bandwidth when the
% value approached 0.01 F. For this value, the bandwidth decreased, evident
% in comparing the frequency spectrums for this case (see figure 19) and
% the previous cases (see figures 17 and 15). The simulation results
% suggest that the bandwidth decreases as the capacitance increases. 

%% 
% Next, the number of timesteps were varied. 

C(4, 4) = 0.00001; 

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/100; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.01:1
    Vold = V;
    vin8(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin8(ii);
    A = (C / 0.01) + G;
    V = A \ ((C * Vold / 0.01) + F');
    
    vo8(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.01 : 1;
figure(20);
plot(t, vin8);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo8);
title('(A) V0 vs t for guassian pulse input with noise (100 timesteps)');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

A = zeros(7); 
F = zeros(1, 7);
ii = 1; 
dt = 1.0/10000; 
Vold = zeros(7,1);
V = zeros(7,1);

for t = 0.001:0.0001:1
    Vold = V;
    vin9(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    In = randn * 0.001;
    F(7) = In; 
    F(1) = vin9(ii);
    A = (C / 0.0001) + G;
    V = A \ ((C * Vold / 0.0001) + F');
    
    vo9(ii) = V(6);
    ii = ii + 1; 
end

t = 0.001 : 0.0001 : 1;
figure(21);
plot(t, vin9);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo9);
title('(A) V0 vs t for guassian pulse input with noise (10000 timesteps)');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
hold off

%% Non-linearity
% In order to solve for a non-linear transconductance equation for the
% voltage source on the output stage, a new B matrix would need to be
% implemented. In addition, a jacobian matrix may be implemented in the
% program for the derivatives, which may present nested loops. Also, a H 
% matrix would be required in order to solve for the V matrix. 
