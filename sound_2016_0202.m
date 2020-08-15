clear all
close all
clc

[audio, fs] = audioread('sound_corrupted.wav');

nfft = 2^12;
overlap = 3/4*nfft;
window = kaiser(nfft, 7);

[s,f,t] = spectrogram(audio, window, overlap, nfft, fs);
s_dB = 20*log10(abs(s));

figure
imagesc(t, f, s_dB), title('Spektrogram signala');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig = gca; fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
fig.YAxis.TickLabelFormat = '%d Hz';
fig.YMinorTick = 'on';
ylim([0 16000]);

% savefig('Spektrogram signala')
% print('Spektrogram signala','-dsvg','-r0')

%NO:
fa1 = 1430;
fa2 = 1570;
fp1 = 1380;
fp2 = 1620;

fa3 = 2600;
fa4 = 2790;
fp3 = 2550;
fp4 = 2840;

%NF:
fp5 = 6400;
fa5 = 6450;

Ap = 0.5;
Aa = 80;

Wa1 = 2*pi*fa1/fs;
Wp1 = 2*pi*fp1/fs; 

Wa2 = 2*pi*fa2/fs; 
Wp2 = 2*pi*fp2/fs;

Wa3 = 2*pi*fa3/fs; 
Wp3 = 2*pi*fp3/fs;

Wa4 = 2*pi*fa4/fs;
Wp4 = 2*pi*fp4/fs;

Wa5 = 2*pi*fa5/fs; 
Wp5 = 2*pi*fp5/fs;

Wa = [Wa1 Wa2 Wa3 Wa4 Wa5];
Wp = [Wp1 Wp2 Wp3 Wp4 Wp5];

dp = (10^(0.05*Ap)-1 )/(10^(0.05*Ap)+1);
da = 10^(-0.05*Aa);

x_g1 = [0 Wp1]; y_g1 = [1-dp 1-dp];
x_g2 = [0 Wp1]; y_g2 = [1+dp 1+dp];
x_g3 = [Wp1 Wp1]; y_g3 = [da/10 1+dp];

x_g4 = [Wa1 Wa2]; y_g4 = [da da];
x_g5 = [Wa1 Wa1]; y_g5 = [da 10*da];
x_g6 = [Wa2 Wa2]; y_g6 = [da 10*da];

x_g7 = [Wp2 Wp3]; y_g7 = [1-dp 1-dp];
x_g8 = [Wp2 Wp3]; y_g8 = [1+dp 1+dp];
x_g9 = [Wp2 Wp2]; y_g9 = [da/10 1+dp];
x_g10 = [Wp3 Wp3]; y_g10 = [da/10 1+dp];

x_g11 = [Wa3 Wa4]; y_g11 = [da da];
x_g12 = [Wa3 Wa3]; y_g12 = [da 10*da];
x_g13 = [Wa4 Wa4]; y_g13 = [da 10*da];

x_g14 = [Wp4 Wp5]; y_g14 = [1-dp 1-dp];
x_g15 = [Wp4 Wp5]; y_g15 = [1+dp 1+dp];
x_g16 = [Wp4 Wp4]; y_g16 = [da/10 1+dp];
x_g17 = [Wp5 Wp5]; y_g17 = [da/10 1+dp];

x_g18 = [Wa5 pi]; y_g18 = [da da];
x_g19 = [Wa5 Wa5]; y_g19 = [da/10 1+dp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h, N] = FIR_optim(Wp, Wa, Ap, Aa);

Nfreqz = 15000; 
[H, w] = freqz(h, 1, Nfreqz);
Ha = abs(H);
Hp = unwrap(angle(H));

figure;
n = 0:N;
stem(n, h);
title('FIR filter: optimization method -> impulse response'), grid
xlabel('n');
xlim([0 N]);

savefig('FIR filter optimization method impulse response')
print('FIR filter optimization method impulse response','-dsvg','-r0')

figure;
stem(n, h);
title('FIR filter: optimization method -> impulse response'), grid
xlabel('n');
xlim([15*(N+1)/32 17*(N+1)/32]);

savefig('FIR filter optimization method impulse response zoom')
print('FIR filter optimization method impulse response zoom','-dsvg','-r0')

figure;
plot(w, 20*log10(Ha));
title('FIR filter: optimization method -> amplitude filter characteristics');
xlabel('W');
ylabel('Magnitude [dB]');
xlim([0 pi]);
grid on
hold on
line(x_g1, 20*log10(y_g1), 'Color', 'r');
line(x_g2, 20*log10(y_g2), 'Color', 'r');
line(x_g3, 20*log10(y_g3), 'Color', 'r');
line(x_g4, 20*log10(y_g4), 'Color', 'r');
line(x_g5, 20*log10(y_g5), 'Color', 'r');
line(x_g6, 20*log10(y_g6), 'Color', 'r');
line(x_g7, 20*log10(y_g7), 'Color', 'r');
line(x_g8, 20*log10(y_g8), 'Color', 'r');
line(x_g9, 20*log10(y_g9), 'Color', 'r');
line(x_g10, 20*log10(y_g10), 'Color', 'r');
line(x_g11, 20*log10(y_g11), 'Color', 'r');
line(x_g12, 20*log10(y_g12), 'Color', 'r');
line(x_g13, 20*log10(y_g13), 'Color', 'r');
line(x_g14, 20*log10(y_g14), 'Color', 'r');
line(x_g15, 20*log10(y_g15), 'Color', 'r');
line(x_g16, 20*log10(y_g16), 'Color', 'r');
line(x_g17, 20*log10(y_g17), 'Color', 'r');
line(x_g18, 20*log10(y_g18), 'Color', 'r');
line(x_g19, 20*log10(y_g19), 'Color', 'r');
hold off

savefig('FIR filter optimization method amplitude')
print('FIR filter optimization method amplitude','-dsvg','-r0')

figure;
plot(w, Hp); 
title('FIR filter: optimization method -> phase filter characteristics');
axis([0 pi min(Hp)-0.1 max(Hp)+0.1])
ylabel('Phase [rad]');
xlabel('W');

savefig('FIR filter optimization method phase')
print('FIR filter optimization method phase','-dsvg','-r0')

figure;
plot(w*fs/2/pi, Hp); 
title('FIR filter: optimization method -> phase filter characteristics');
ylabel('Phase [rad]');
xlabel('Frequency [Hz]');


savefig('FIR filter optimization method phase2')
print('FIR filter optimization method phase2','-dsvg','-r0')

audio_filtered1 = filter(h, 1, audio);

[s,f,t] = spectrogram(audio_filtered1, window, overlap, nfft, fs);
s1_dB = 20*log10(abs(s));

figure
imagesc(t, f, s1_dB);
title('FIR filter: optimization method -> spectrogram');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig = gca; fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
fig.YAxis.TickLabelFormat = '%d Hz';
fig.YMinorTick = 'on';
ylim([0 16000]);

savefig('FIR filter optimization method spectrogram')
print('FIR filter optimization method spectrogram','-dsvg','-r0')


audiowrite('out_signal1_2015_0202.wav', audio_filtered1, fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h, N] = FIR_limit(Wp, Wa, Ap, Aa);

Nfreqz = 15000;
[H, w] = freqz(h, 1, Nfreqz);
Ha = abs(H);
Hp = unwrap(angle(H));

figure;
n = 0:N;
stem(n, h);
title('FIR filter: limited impulse response method -> impulse response'), grid
xlabel('n');
xlim([0 N]);

savefig('FIR filter limited impulse response method impulse response')
print('FIR filter limited impulse response method impulse response','-dsvg','-r0')

figure;
stem(n, h);
title('FIR filter: limited impulse response method -> impulse response'), grid
xlabel('n');
xlim([15*(N+1)/32 17*(N+1)/32]);

savefig('FIR filter limited impulse response method impulse response zoom')
print('FIR filter limited impulse response method impulse response zoom','-dsvg','-r0')

figure;
plot(w, 20*log10(Ha));
title('FIR filter: limited impulse response method -> amplitude filter characteristics');
xlabel('W');
ylabel('Magnitude [dB]');
xlim([0 pi]);
grid on
hold on
line(x_g1, 20*log10(y_g1), 'Color', 'r');
line(x_g2, 20*log10(y_g2), 'Color', 'r');
line(x_g3, 20*log10(y_g3), 'Color', 'r');
line(x_g4, 20*log10(y_g4), 'Color', 'r');
line(x_g5, 20*log10(y_g5), 'Color', 'r');
line(x_g6, 20*log10(y_g6), 'Color', 'r');
line(x_g7, 20*log10(y_g7), 'Color', 'r');
line(x_g8, 20*log10(y_g8), 'Color', 'r');
line(x_g9, 20*log10(y_g9), 'Color', 'r');
line(x_g10, 20*log10(y_g10), 'Color', 'r');
line(x_g11, 20*log10(y_g11), 'Color', 'r');
line(x_g12, 20*log10(y_g12), 'Color', 'r');
line(x_g13, 20*log10(y_g13), 'Color', 'r');
line(x_g14, 20*log10(y_g14), 'Color', 'r');
line(x_g15, 20*log10(y_g15), 'Color', 'r');
line(x_g16, 20*log10(y_g16), 'Color', 'r');
line(x_g17, 20*log10(y_g17), 'Color', 'r');
line(x_g18, 20*log10(y_g18), 'Color', 'r');
line(x_g19, 20*log10(y_g19), 'Color', 'r');
hold off

savefig('FIR filter limited impulse response method amplitude')
print('FIR filter limited impulse response method amplitude','-dsvg','-r0')

figure;
plot(w, Hp); 
title('FIR filter: limited impulse response method -> phase filter characteristics');
axis([0 pi min(Hp)-0.1 max(Hp)+0.1])
ylabel('Phase [rad]');
xlabel('W');

savefig('FIR filter limited impulse response method phase')
print('FIR filter limited impulse response method phase','-dsvg','-r0')

figure;
plot(w*fs/2/pi, Hp); 
title('FIR filter: limited impulse response method -> phase filter characteristics');
ylabel('Phase [rad]');
xlabel('Frequency [Hz]');

savefig('FIR filter limited impulse response method phase2')
print('FIR filter limited impulse response method phase2','-dsvg','-r0')

audio_filtered2 = filter(h, 1, audio);

[s,f,t] = spectrogram(audio_filtered2, window, overlap, nfft, fs);
s1_dB = 20*log10(abs(s));

figure
imagesc(t, f, s1_dB);
title('FIR filter: limited impulse response method -> spectrogram');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig = gca; fig.XAxis.TickLabelFormat = '%.1f s';
fig.XMinorTick = 'on';
fig.YAxis.TickLabelFormat = '%d Hz';
fig.YMinorTick = 'on';
ylim([0 16000]);

savefig('FIR filter limited impulse response method spectrogram')
print('FIR filter limited impulse response method spectrogram','-dsvg','-r0')

audiowrite('out_signal2_2015_0202.wav', audio_filtered2, fs);
