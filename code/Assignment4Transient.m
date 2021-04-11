%% Part 2: Transient Circuit Simulation
% In this part, the circuit was analyzed in the time domain for various
% inputs for a time period of 1 second. The first input was a simple unit 
% step input that goes high after 0.03 seconds.

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

vin = zeros(1, 1000);
vo = zeros(1, 1000);
v3 = zeros(1, 1000);

G(1, 1) = 1;                                    % 1
G(2, 1) = -G1; G(2, 2) = G1 + G2;               % 2
G(3 ,2) = -1; G(3, 4) = 1;                      % iL
G(4, 3) = -1; G(4, 4) = G3;                     % 3
G(5, 5) = 1; G(5, 4) = -alpha*G3;               % 4
G(6, 6) = G4 + G0; G(6, 5) = -G4;               % 5

C = zeros(6);

C(2, 1) = -Cval; C(2, 2) = Cval;
C(3, 3) = L;
 
F = zeros(1, 6);

ii = 1; 


V = zeros(6,1);
for t = 0.001:0.001:1
    Vold = V;
    if t < 0.03
        vin(ii) = 0; 
        F(1) = 1; 
    else
        vin(ii) = 1;
        F(1) = 0; 
    end
    
    F(1) = vin(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo(ii) = V(6);
    %v3(ii) = V(4);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(6);
plot(t, vin);
hold on
plot(t, vo);
title('(A) V0 vs t for unit step input');
ylabel('V0 (red) Vin (blue)');
xlabel('t');
ylim([-2 10])


%%
% The second input was a sinusoidal input at a frequency of 1/0.03 Hz. 

V = zeros(6,1);
F = zeros(1, 6);
ii = 1; 

for t = 0.001:0.001:1
    Vold = V;
    vin2(ii) = sin(2*pi*(1/0.03) * t);
    
    F(1) = vin2(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo2(ii) = V(6);
    %v3(ii) = V(4);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(7);
plot(t, vin2);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo2);
title('(A) V0 vs t for sinusoidal input');
ylabel('V0 (red) Vin (blue)');
xlabel('t');

%%
% The third input was a guassian pulse with a magnitude of 1, std dev. of 
% 0.03 seconds and a delay of 0.06 seconds.

A = zeros(6); 
F = zeros(1, 6);
v = -10;
ii = 1; 
dt = 1.0/1000; 
Vold = zeros(6,1);
V = zeros(6,1);

for t = 0.001:0.001:1
    Vold = V;
    vin3(ii) = exp(-0.5 * power(((t - 0.06) / (0.03)), 2)); % guassian pulse
    
    F(1) = vin3(ii);
    A = (C / 0.001) + G;
    V = A \ ((C * Vold / 0.001) + F');
    
    vo3(ii) = V(6);
    %v3(ii) = V(4);
    ii = ii + 1; 
end

t = 0.001 : 0.001 : 1;
figure(8);
plot(t, vin3);
%title('(A) V1 vs t');
%ylabel('V1');
xlabel('t');
hold on
plot(t, vo3);
title('(A) V0 vs t for guassian pulse input');
ylabel('V0 (red) Vin (blue)');
xlabel('t');

%%
% Next, the frequency content of each output and input that was analyzed in
% the previous part was plotted. 
 
% Frequency content for each input type:
% Unit step:
fo = fft(vo);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(9);
plot(fs,p);
title('Frequency Content of Vo for Unit Step Input');
hold on
fo = fft(vin);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Unit Step Input');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

% Sine:
fo = fft(vo2);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(10);
plot(fs,p);
hold on
fo = fft(vin2);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Sinusoidal Input');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off

% Guassian
fo = fft(vo3);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
figure(11);
plot(fs,p);
hold on
fo = fft(vin3);
fs = (-1000/2:1000/2-1); 
p = abs(fftshift(fo)) .^ 2/1000;     
plot(fs,p);
title('Frequency Content of Vo & Vin for Guassian pulse Input');
ylabel('Magnitude: Vin (red) Vo (blue)')
hold off
