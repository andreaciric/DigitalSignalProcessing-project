clear all
close all
clc

%% 3

fs = 600;
f1 = 60;
f2 = 120;
f3 = 180;
f4 = 240;

F1 = f1/fs;
F2 = f2/fs;
F3 = f3/fs;
F4 = f4/fs;

N = 1024; 
n = 0:N-1;

x = sin(2*pi*n*F1) + sin(2*pi*n*F2) + sin(2*pi*n*F3)+ sin(2*pi*n*F4);
x = x/(max([abs(min(x)) abs(max(x))])); % normalizacija

figure;
plot(n, x);
title('Input signal');
xlabel('n');
axis([0 100 min(x)-0.5 max(x)+0.5]);

fc = 60;
Aa = 40; Ap = 0.5;
[b, a] = power_line_noise_filter(fs, fc, Aa, Ap);

[sos] = tf2sos(b, a);
bsos = sos(:, 1:end/2);
asos = sos(:, (end/2+1):end);

y = IIR_direct_II_cascade(bsos, asos, x);

Nfreqz = 4096; 
[H, w_in] = freqz(x, Nfreqz); 
Ha = abs(H);
[H1, w_out] = freqz(y, Nfreqz); 
Ha1 = abs(H1);

figure;
plot(w_in*fs/2/pi, Ha);
title('Input signal - amplitude');
xlabel('w');
ylabel('|H(e^(jW))|');

figure;
plot(w_out*fs/2/pi, Ha1);
title('Output signal - amplitude');
xlabel('w');
ylabel('|H(e^(jW))|');

%% 4.a

bit_length = 2;
frac_length = 1; 
word_length = frac_length + bit_length;

while 1

    FixedPointAttributes = fimath ( 'RoundingMethod', 'Floor', 'OverflowAction', 'Saturate', ...
        'ProductMode', 'SpecifyPrecision', 'ProductWordLength', 2*word_length, 'ProductFractionLength', 2*frac_length, ...
        'SumMode', 'SpecifyPrecision', 'SumWordLength', 2*word_length, 'SumFractionLength', 2*frac_length) ;
               
    FI_a = fi(asos, true, word_length, frac_length, FixedPointAttributes);

    FixedPointAttributes.OverflowAction = 'Wrap';
    FI_a = fi(FI_a, true, word_length, frac_length, FixedPointAttributes);
    
    FI_zp = zeros(2, length(FI_a));
    for i = 1:length(FI_a)
        FI_zp(:,i) = roots(double(FI_a(i,:)));
    end
    
    max_pole = max(abs(double(FI_zp(:,:))));
    if max_pole >= 1
        frac_length = frac_length+1;
        word_length = frac_length+bit_length;
    else        
        disp('Potreban broj bita za razlomljeni deo: ')
        disp(frac_length);
        break
    end
    
end

FI_a = fi (asos, true, word_length, frac_length, FixedPointAttributes);
FI_b = fi (bsos, true, word_length, frac_length, FixedPointAttributes);
FI_zp = zeros(2, length(FI_a));
FI_zo = zeros(2, length(FI_a));

FI_a1 = fi (asos, true, word_length-1, frac_length, FixedPointAttributes);
FI_b1 = fi (bsos, true, word_length-1, frac_length, FixedPointAttributes);
FI_zp1 = zeros(2, length(FI_a));
FI_zo1 = zeros(2, length(FI_a));

for i = 1:length(FI_a)
        FI_zo(:,i) = roots(double(FI_b(i,:)));
        FI_zp(:,i) = roots(double(FI_a(i,:)));
        
        FI_zo1(:,i) = roots(double(FI_b1(i,:)));
        FI_zp1(:,i) = roots(double(FI_a1(i,:)));
end
    
figure;
subplot(121);
[hz, hp, ht] = zplane(FI_zo, FI_zp); 
title('zeros and poles');
xlabel('Re(z)');
ylabel('Im(z)');
subplot(122);
[hz1, hp1, ht1] = zplane(FI_zo1, FI_zp1); 
title('zeros and poles - not stable');
xlabel('Re(z)');
ylabel('Im(z)');

%% 4.b

frac_length1 = 10; 
bit_length1 = frac_length1+bit_length;

while 1
    
    FixedPointAttributes = fimath ( 'RoundingMethod', 'Floor', 'OverflowAction', 'Saturate', ...
            'ProductMode', 'SpecifyPrecision', 'ProductWordLength', 2*bit_length1, 'ProductFractionLength', 2*frac_length1, ...
            'SumMode', 'SpecifyPrecision', 'SumWordLength', 2*bit_length1, 'SumFractionLength', 2*frac_length1) ;
    
    FI_a = fi(asos, true, bit_length1, frac_length1, FixedPointAttributes);
    FI_b = fi(bsos, true, bit_length1, frac_length1, FixedPointAttributes);

    FixedPointAttributes.OverflowAction = 'Wrap';
    FI_a = fi(FI_a, true, bit_length1, frac_length1, FixedPointAttributes);
    FI_b = fi(FI_b, true, bit_length1, frac_length1, FixedPointAttributes);
    
    Nfreqz = 4096; 
    [H, w1] = freqz(b, a, Nfreqz); 
    Ha = abs(H);
    
    sos = [double(FI_b) double(FI_a)];
    [H1, w2] = freqz(sos, Nfreqz); 
    Ha1 = abs(H1);
    
    dHa = max(abs(Ha - Ha1));
    
    if dHa >= 0.01
        frac_length1 = frac_length1+1;
        bit_length1 = frac_length1+bit_length;
    else
        disp('Minimalan broj bita za razlomljeni deo: ')
        disp(frac_length1);
        break
    end
    
end

figure;
plot(w1, Ha);
title('Amplitude filter characteristics'), grid on, hold on
plot(w2, Ha1, 'r');
xlabel('w');
ylabel('|H(e^(jw))|');
ylim([0.85 1.05]);
xlim([0.5 0.75]);
