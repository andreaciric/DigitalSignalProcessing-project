close all
clear all
clc

%% input signal

load('ecg_corrupted.mat');  %loading ecg file
fs = 360;                   %sampling frequency (Hz)
t = (1:length(val))/fs;     %time

figure
plot(t,val)
title('ECG input signal - time domain')
xlim([0 (length(t)-1)/fs])
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%d s';
fig.XMinorTick = 'on';
% savefig('ECG input signal - time domain')
% print('ECG input signal - time domain','-dsvg','-r0')

figure
plot(t,val)
title('ECG input signal - time domain')
xlim([0 (length(t)-1)/fs/12])
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
% savefig('ECG input signal - time domain zoom')
% print('ECG input signal - time domain zoom','-dsvg','-r0')


val_ft = fft(val);
 
figure
stem((0:length(val_ft)/2-1)/length(val_ft)*fs,val_ft(1:length(val_ft)/2))
title('ECG input signal - frequency domain')
xlim([0 95])
xlabel('frequency')
fig = gca;
fig.XAxis.TickLabelFormat = '%d Hz';
fig.XMinorTick = 'on';
% savefig('ECG input signal - frequency domain')
% print('ECG input signal - frequency domain','-dsvg','-r0')

figure
stem((0:length(val_ft)/2-1)/length(val_ft)*fs,val_ft(1:length(val_ft)/2))
title('ECG input signal - frequency domain')
xlim([54 66])
ylim([-5 5])
xlabel('frequency')
fig = gca;
fig.XAxis.TickLabelFormat = '%d Hz';
fig.XMinorTick = 'on';
hold on
line([60-2 60-2],[-5 5],'Color','red','LineStyle','--')
line([60+2 60+2],[-5 5],'Color','red','LineStyle','--')
line([60-0.5 60-0.5],[-2.5 2.5],'Color','red','LineStyle','--')
line([60+0.5 60+0.5],[-2.5 2.5],'Color','red','LineStyle','--')
hold off
% savefig('ECG input signal - frequency domain zoom 2')
% print('ECG input signal - frequency domain3 zoom 2','-dsvg','-r0')

figure
stem((0:length(val_ft)/2-1)/length(val_ft)*fs,val_ft(1:length(val_ft)/2))
title('ECG input signal - frequency domain')
xlim([0 6])
xlabel('frequency')
fig = gca;
fig.XAxis.TickLabelFormat = '%d Hz';
fig.XMinorTick = 'on';
hold on
line([0.5 0.5],[-100 150],'Color','red','LineStyle','--')
hold off
% savefig('ECG input signal - frequency domain zoom 1')
% print('ECG input signal - frequency domain zoom 1','-dsvg','-r0')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% baseline drift filter - highpass IIR filter

fa = 0.4; 
fp = 1; 
wa = 2*pi*fa;
wp = 2*pi*fp;
Aa = 30;
Ap = 0.5;

[b, a] = baseline_drift_filter(fs, fa, fp, Aa, Ap);

y_baseline_drift_filtered = filter(b, a, val);

figure
plot(t,y_baseline_drift_filtered)
title('ECG signal after baseline drift filter - time domain')
xlim([0 (length(t)-1)/fs])
ylim([-1.5 1.5]);
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%d s';
fig.XMinorTick = 'on';
% savefig('ECG signal after baseline drift filter - time domain')
% print('ECG signal after baseline drift filter - time domain','-dsvg','-r0')

figure
plot(t,y_baseline_drift_filtered)
title('ECG signal after baseline drift filter - time domain')
xlim([0 (length(t)-1)/fs/12])
ylim([-1.5 1.5]);
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%d s';
fig.XMinorTick = 'on';
% savefig('ECG signal after baseline drift filter - time domain zoom')
% print('ECG signal after baseline drift filter - time domain zoom','-dsvg','-r0')

y_baseline_drift_filtered_ft = fft(y_baseline_drift_filtered);
 
figure
stem((0:length(y_baseline_drift_filtered_ft)/2-1)/length(y_baseline_drift_filtered_ft)*fs,y_baseline_drift_filtered_ft(1:length(y_baseline_drift_filtered_ft)/2))
title('ECG signal after baseline drift filter - frequency domain')
xlim([0 6])
xlabel('frequency')
fig = gca;
fig.XAxis.TickLabelFormat = '%d Hz';
fig.XMinorTick = 'on';
hold on
line([0.5 0.5],[-100 100],'Color','red','LineStyle','--')
hold off
% savefig('ECG signal after baseline drift filter - frequency domain')
% print('ECG signal after baseline drift filter - frequency domain','-dsvg','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% power line noise filter - highpass IIR filter

fc = 60;
Aa1 = 40; 
Ap1 = 0.5;
 
[b1, a1] = power_line_noise_filter(fs, fc, Aa1, Ap1);

y_power_line_noise_filtered = filter(b1, a1, y_baseline_drift_filtered);

figure
plot(t,y_power_line_noise_filtered)
title('ECG signal after power line noise filter - time domain')
xlim([0 (length(t)-1)/fs])
ylim([-1.5 1.5]);
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%d s';
fig.XMinorTick = 'on';
% savefig('ECG signal after power line noise filter - time domain')
% print('ECG signal after power line noise filter - time domain','-dsvg','-r0')

figure
plot(t,y_power_line_noise_filtered)
title('ECG signal after power line noise filter - time domain')
xlim([0 (length(t)-1)/fs/12])
ylim([-1.5 1.5]);
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
% savefig('ECG signal after power line noise filter - time domain zoom')
% print('ECG signal after power line noise filter - time domain zoom','-dsvg','-r0')

y_power_line_noise_filtered_ft = fft(y_power_line_noise_filtered);
 
figure
stem((0:length(y_power_line_noise_filtered_ft)/2-1)/length(y_power_line_noise_filtered_ft)*fs,y_power_line_noise_filtered_ft(1:length(y_power_line_noise_filtered_ft)/2))
title('ECG signal after power line noise filter - frequency domain')
xlim([54 66])
xlabel('frequency')
fig = gca;
fig.XAxis.TickLabelFormat = '%d Hz';
fig.XMinorTick = 'on';
hold on
line([60-2 60-2],[-5 5],'Color','red','LineStyle','--')
line([60+2 60+2],[-5 5],'Color','red','LineStyle','--')
line([60-0.5 60-0.5],[-2.5 2.5],'Color','red','LineStyle','--')
line([60+0.5 60+0.5],[-2.5 2.5],'Color','red','LineStyle','--')
hold off
% savefig('ECG signal after power line noise filter - frequency domain')
% print('ECG signal after power line noise filter - frequency domain','-dsvg','-r0')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% baseline drift filter characteristics

zo = roots(b);
zp = roots(a);
Wp = wp/fs; 
Wa = wa/fs;

figure;
[hz, hp, ht] = zplane(zo, zp);
title('baseline drift filter - zeros and poles')
xlabel('Re(z)');
ylabel('Im(z)');
% savefig('baseline drift filter - zeros and poles')
% print('baseline drift filter - zeros and poles','-dsvg','-r0')

figure;
[hz, hp, ht] = zplane(zo, zp);
title('baseline drift filter - zeros and poles')
xlabel('Re(z)');
ylabel('Im(z)');
xlim([0.97 1.03]);
ylim([-0.02 0.02]);
% savefig('baseline drift filter - zeros and poles zoom')
% print('baseline drift filter - zeros and poles zoom','-dsvg','-r0')

Nfreqz = 15000; 
[h, w] = freqz(b, a, Nfreqz);
H = abs(h);

figure;
plot(w, 20*log10(H));
title('baseline drift filter - magnitude');
xlabel('W');
ylabel('Magnitude [dB]');
xlim([0 0.03]);
ylim([-70 10]);
grid on
fig = gca;
fig.YAxis.TickLabelFormat = '%d dB';
fig.YMinorTick = 'on';
hold on
line([Wp Wp], [-Ap 0], 'Color', 'r');
line([Wa Wa], [-70 -Aa], 'Color', 'r');
line([Wp 0.03], [-Ap -Ap], 'Color', 'r');
line([0 Wa], [-Aa -Aa], 'Color', 'r');
hold off
% savefig('baseline drift filter - magnitude')
% print('baseline drift filter - magnitude','-dsvg','-r0')

figure()
semilogx(w,H), grid
xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('baseline drift filter - magnitude');
% savefig('baseline drift filter - magnitude2')
% print('baseline drift filter - magnitude2','-dsvg','-r0')

H_phase = unwrap(angle(h));

figure;
plot(w, H_phase); 
title('baseline drift filter - phase');
ylabel('Phase [rad]');
xlabel('W');
xlim([0 3.14]);
% savefig('baseline drift filter - phase')
% print('baseline drift filter - phase','-dsvg','-r0')

figure;
plot(w*fs/2/pi, H_phase); 
title('baseline drift filter - phase');
ylabel('Phase [rad]');
xlabel('Frequency [Hz]');
xlim([0 5]);
% savefig('baseline drift filter - phase2')
% print('baseline drift filter - phase2','-dsvg','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% power line noise drift filter characteristics

zo1 = roots(b1);
zp1 = roots(a1);

fp1 = fc - 2; 
fp2 = fc + 2;
fa1 = fc - 0.5; 
fa2 = fc + 0.5;
Wp1 = 2*pi*fp1/fs; 
Wp2 = 2*pi*fp2/fs;
Wa1 = 2*pi*fa1/fs; 
Wa2 = 2*pi*fa2/fs;

figure;
[hz1, hp1, ht1] = zplane(zo1, zp1); 
title('power line noise filter - zeros and poles')
xlabel('Re(z)');
ylabel('Im(z)');
% savefig('power line noise filter - zeros and poles')
% print('power line noise filter - zeros and poles','-dsvg','-r0')

figure;
[hz1, hp1, ht1] = zplane(zo1, zp1); 
title('power line noise filter - zeros and poles')
xlabel('Re(z)');
ylabel('Im(z)');
xlim([0.4 0.6]);
ylim([0.76 0.96]);
% savefig('power line noise filter - zeros and poles zoom')
% print('power line noise filter - zeros and poles zoom','-dsvg','-r0')

[h1, w1] = freqz(b1, a1, Nfreqz);
H1 = abs(h1);

figure;
plot(w1, 20*log10(H1));
title('power line noise filter - magnitude');
xlabel('W');
ylabel('Magnitude [dB]');
xlim([0.9 1.2]);
ylim([-80 10]);
grid on
hold on
line([Wa1 Wa1], [-Aa1 -80], 'Color', 'r');
line([Wa2 Wa2], [-Aa1 -80], 'Color', 'r');
line([Wp1 Wp1], [-Ap1 0], 'Color', 'r');
line([Wp2 Wp2], [-Ap1 0], 'Color', 'r');
line([Wa1 Wa2], [-Aa1 -Aa1], 'Color', 'r');
line([0 Wp1], [-Ap1 -Ap1], 'Color', 'r');
line([Wp2 1.2], [-Ap1 -Ap1], 'Color', 'r');
hold off
% savefig('power line noise filter - magnitude')
% print('power line noise filter - magnitude','-dsvg','-r0')

H1_phase = unwrap(angle(h1));

figure;
plot(w, H1_phase); 
title('power line noise filter - phase');
ylabel('Phase [rad]');
xlabel('W');
xlim([0 3.14]);
% savefig('power line noise filter - phase')
% print('power line noise filter - phase','-dsvg','-r0')

figure;
plot(w*fs/2/pi, H1_phase); 
title('power line noise filter - phase');
ylabel('Phase [rad]');
xlabel('Frequency [Hz]');
% xlim([0 3.14]);
% savefig('power line noise filter - phase2')
% print('power line noise filter - phase2','-dsvg','-r0')

smtlb = sgolayfilt(y_power_line_noise_filtered,5,23);
figure;
plot(t, smtlb); 
title('ECG signal after Savitzky-Golay filtering - time domain')
xlim([0 (length(t)-1)/fs/12])
ylim([-1.5 1.5]);
xlabel('time')
fig = gca;
fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
savefig('ECG signal after Savitzky-Golay filtering - time domain')
print('ECG signal after Savitzky-Golay filtering - time domain','-dsvg','-r0')