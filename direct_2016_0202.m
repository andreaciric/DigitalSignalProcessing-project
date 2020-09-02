clear all
close all
clc

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

y = IIR_direct_II(b, a, x);

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

bit_length = 5;

frac_length = 10; 
word_length = frac_length + bit_length;

while 1

    FixedPointAttributes = fimath ( 'RoundingMethod', 'Floor', 'OverflowAction', 'Saturate', ...
        'ProductMode', 'SpecifyPrecision', 'ProductWordLength', 2*word_length, 'ProductFractionLength', 2*frac_length, ...
        'SumMode', 'SpecifyPrecision', 'SumWordLength', 2*word_length, 'SumFractionLength', 2*frac_length) ;
               
    FI_a = fi(a, true, word_length, frac_length, FixedPointAttributes);

    FixedPointAttributes.OverflowAction = 'Wrap';
    FI_a = fi(FI_a, true, word_length, frac_length, FixedPointAttributes);
    
    FI_zp = roots(double(FI_a));
    
    max_pole = max(abs(double(FI_zp)));
    if max_pole >= 1
        frac_length = frac_length+1;
        word_length = frac_length+bit_length;
    else        
        disp('Potreban broj bita za razlomljeni deo: ')
        disp(frac_length);
        break
    end
    
end

FI_a = fi (a, true, word_length, frac_length, FixedPointAttributes);
FI_b = fi (b, true, word_length, frac_length, FixedPointAttributes);
FI_zp = roots(double(FI_a));
FI_zo = roots(double(FI_b));

FI_a1 = fi (a, true, word_length-1, frac_length, FixedPointAttributes);
FI_b1 = fi (b, true, word_length-1, frac_length, FixedPointAttributes);
FI_zp1 = roots(double(FI_a1));
FI_zo1 = roots(double(FI_b1));

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
    
    FI_a = fi(a, true, bit_length1, frac_length1, FixedPointAttributes);
    FI_b = fi(b, true, bit_length1, frac_length1, FixedPointAttributes);

    FixedPointAttributes.OverflowAction = 'Wrap';
    FI_a = fi(FI_a, true, bit_length1, frac_length1, FixedPointAttributes);
    FI_b = fi(FI_b, true, bit_length1, frac_length1, FixedPointAttributes);
    
    Nfreqz = 4096; 
    [H, w1] = freqz(b, a, Nfreqz); 
    Ha = abs(H);
    
    [H1, w2] = freqz(double(FI_b), double(FI_a), Nfreqz); 
    Ha1 = abs(H1);
    
    dHa = max(abs(Ha - Ha1));
    
    if dHa >= 0.02
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

%% 5

PWL = 63;
PFL = 46;
SWL = 63;
SFL = 46;

FCB = 28;
FCF = 23;
SB = 24;
SF = 23;

FixedPointAttributes = fimath ( 'RoundingMethod', 'Floor', 'OverflowAction', 'Saturate', ...
        'ProductMode', 'SpecifyPrecision', 'ProductWordLength', PWL, 'ProductFractionLength', PFL, ...
        'SumMode', 'SpecifyPrecision', 'SumWordLength', SWL, 'SumFractionLength', SFL) ;
fi_params = struct('FILTER_COEFITIENTS_BITLENGTH', FCB, 'FILTER_COEFITIENTS_FRAC', FCF, ...
                   'SIGNAL_BITLENGTH', SB, 'SIGNAL_FRAC', SF);

FI_a = fi(a, true, fi_params.FILTER_COEFITIENTS_BITLENGTH, fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_b = fi(b, true, fi_params.FILTER_COEFITIENTS_BITLENGTH, fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_x = fi(x, true, fi_params.SIGNAL_BITLENGTH, fi_params.SIGNAL_FRAC, FixedPointAttributes);

FixedPointAttributes.OverflowAction = 'Wrap';
FI_a = fi(FI_a, true, fi_params.FILTER_COEFITIENTS_BITLENGTH, fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_b = fi(FI_b, true, fi_params.FILTER_COEFITIENTS_BITLENGTH, fi_params.FILTER_COEFITIENTS_FRAC, FixedPointAttributes);
FI_x = fi(FI_x, true, fi_params.SIGNAL_BITLENGTH, fi_params.SIGNAL_FRAC, FixedPointAttributes);

y = IIR_direct_II(b, a, x);

y_fixed_point = FI_IIR_direct_II(FI_b, FI_a, FI_x);

figure
subplot(411);
plot(n, x);
title('Ulazni signal');
xlabel('n');
xlim([0 1050]);

subplot(412);
plot(n, y);
title('Izlazni signal sa double preciznoscu');
xlabel('n');
xlim([0 1050]);

subplot(413);
plot(n, y_fixed_point);
title('Izlazni signal sa fixed-point preciznoscu');
xlabel('n');
xlim([0 1050]);

err = y - double(y_fixed_point);
subplot(414);
plot(n, err);
title('Razlika izlaznih signala');
xlabel('n');
xlim([0 1050]);


