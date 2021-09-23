% Filter implementation 
% Filter specifiations
% index no : 180427n
clc;
clear all;
A = 4;
B = 2;
C = 7;
Ap  = 0.03 + (0.01*A);  % Maximum passband ripple, dB
Aa  = 45 + B;           % Maximum stopband ripple, dB
wp1 = C*100 + 300;      % Lower passband edge, rad/s
wp2 = C*100 + 700;      % Upper passband edge, rad/s
wa1 = C*100 + 150;      % Lower stopband edge, rad/s
wa2 = C*100 + 800;      % Upper stopband edge, rad/s
ws  = 2*(C*100 + 1200); %sampling frequency rad/s
T   = 2*pi/ws;          %sampling period  

bt1 = wp1 - wa1;        %Lower transition width
bt2 = wa2 - wp2;        %Upper transition width
bt  = min(bt1, bt2);    %Critical transition width 
wc1 = wp1 - bt/2;       %Lower cutoff frequency
wc2 = wp2 + bt/2;       %Upper cutoff frequency

%deriving delta
delta_p = (10^(0.05*Ap) - 1)/(10^(0.05*Ap) + 1);
delta_a = 10^(-0.05*Aa);
delta   = min(delta_p, delta_a);

%Actual stopband attenuation 
act_Aa = 20*log10(1/delta);

%calculating alpha
if act_Aa <= 21
    alpha = 0;
elseif act_Aa > 50
    alpha = 0.1102*(act_Aa - 8.7);
else
    alpha = 0.5842*(act_Aa-21)^0.4 + 0.07886*(act_Aa-21);
end

%calculating D
if act_Aa <= 21
    D = 0.9222;
else
    D = (act_Aa - 7.95)/14.36; 
end 

% calculating lowest odd value for N(lenth of the filter)
N = ceil(ws*D/bt + 1);
if mod(N,2) == 0
    N = N + 1;
end 
disp(N);
% determining kaiser window parameters 
nmax = (N - 1)/2; %maximum time steps (because of symmetry only..
               ... considered positive part)

kwn = zeros(1, nmax+ 1);   % row vector to store kaiser window value

% obtaining Io(alpha)
Io_alpha = 1;
k = 1;
while true
   bassel = (1/factorial(k) * (alpha/2)^k)^2; % current bassel value
   Io_alpha = Io_alpha + bassel;
   k = k + 1;
   if bassel < 10e-6      % bassel limit taught in lecture 
       break
   end 
end

% obtaining Io(beta)
for k = 1:nmax+1 
    beta = alpha * sqrt(1 - (2*(k-1)/(N-1))^2); % calculating beta
    Io_beta = 1;    % initializing 
    j = 1;
    while true
        bassel = (1/factorial(j) * (beta/2)^j)^2;
        Io_beta = Io_beta + bassel;
        j = j + 1;
        if bassel < 10e-6      % bassel limit taught in lecture 
           break
        end 
    end
    kwn(k) = Io_beta/Io_alpha; % uptataing window rowvector
end

figure 
%[fkw, f] = freqz(kwn);
stem(0:nmax, kwn);

% idealized impulse response of filter 
hi_nT = zeros(1,nmax+ 1); 
for n = 1:nmax+1
    if n == 1
        hi_nT(n) = (2/ws)*(wc2 - wc1);
    else
        hi_nT(n) = ((sin(wc2*(n-1)*T) - sin(wc1*(n-1)*T))/((n-1)*pi)); 
    end
end
figure;
stem(0:nmax, hi_nT);
    
% determining causal impulse response 
htemp = hi_nT .* kwn;
hn = [fliplr(htemp(2:end)), htemp]; %fliping positive value to negative 

figure;
stem(0:2*nmax, hn);
title("Causal impulse response of the filter"); 
xlim([0 N-1]); grid on;
ylabel("Amplitude"), xlabel("\it{n} (samples)");

% frequency response of the filter 
[freq_res, f] = freqz(hn, 1, 2048);  % f - rad/sample H(w) = freq_res
w = (f*ws)/(2*pi);                   % w - rad/sec

% magnitute response in the range of (0, ws/2)
mag_freq_res = 20*log10(abs(freq_res));
figure;
plot(w, mag_freq_res);
title('Magnitude Response in the range (0,ws/2)'); 
xlabel('Angular Frequency (rad/s)');
ylabel('Magnitude (dB)'); 
axis([0, ws/2, -90, 20]); 
grid on;

% magnitude response of the passband 
figure;
plot(w, mag_freq_res);
title('Magnitude response of the passband '); 
xlabel('Angular Frequency (rad/s)');
ylabel('Magnitude (dB)'); 
axis([wc1, wc2, -Ap/2, Ap/2]); 
grid on;

% Input signal generation 
w1 = wa1/2;
w2 = wp1 + (wp2 - wp1)/2;
w3 = wa2 + (ws/2 - wa2)/2;

samples = 300;
n = 0:samples;
xnT = sin(w1*n*T) + sin(w2*n*T) + sin(w3*n*T);

% Expected output signal 
yexp = sin(w2*n*T);

% Filtered signal  
N_point = length(xnT) + length(hn) -1;
Xw = fft(xnT,N_point);
Hw = fft(hn, N_point);
filtered_out = ifft(Xw.*Hw);

figure;
subplot(3,1,1); stem(n, xnT); title("Input Signal");
xlabel('n (samples)'); ylabel('x(nT)'); 
axis([0, (samples+1), -4, 4]);

subplot(3,1,2); stem(n, yexp); title("Expected output Signal");
xlabel('n (samples)'); ylabel('yexp(nT)'); 
axis([0, (samples+1), -1.5, 1.5]);

subplot(3,1,3); stem(n, filtered_out(1:samples+1)); 
title("Obtained output from designed filter");
xlabel('n (samples)'); ylabel('yo(nT)'); 
axis([0, (samples+1), -1.5, 1.5]);

% freqency domain analisis 
Xw = fft(xnT); 
n1 = length(Xw);
Yexpw = fft(yexp);
f_n = (0 : n1-1)/(ws/n1);%freqency vector for input and expected output

Xshift = fftshift(Xw); 
Yshift = fftshift(Yexpw);
fshift = (-n1/2 : n1/2 -1)*(ws/n1);

YoW = fft(filtered_out, N_point);
Ye_shift = fftshift(YoW);
fe_shift = (-N_point/2 : N_point/2 -1)*(ws/N_point); % freq vector 

figure;
subplot(3,1,1); plot(fshift, abs(Xshift));
title("fft of Input Signal");
xlabel('frequency (rad/s)'); ylabel('X(w)'); 

subplot(3,1,2); plot(fshift, abs(Yshift)); 
title("fft of expected output Signal");
xlabel('frequency (rad/s)'); ylabel('Yexp(w)'); 

subplot(3,1,3); plot(fe_shift, abs(Ye_shift)); 
title("fft of Obtained output from designed filter");
xlabel('frequency (rad/s)'); ylabel('Yo(w)');  



