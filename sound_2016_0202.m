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
imagesc(t, f, s_dB), title('Spektogram signala');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig2 = gca; fig2.XAxis.TickLabelFormat = '%.1f s';
fig2.XMinorTick = 'on';
fig2.YAxis.TickLabelFormat = '%d Hz';
fig2.YMinorTick = 'on';
ylim([0 16000]);

savefig('Spektogram signala')
print('Spektogram signala','-dsvg','-r0')

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
line([Wp1 Wp1], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa1 Wa1], 20*log10([da 10*da]), 'Color', 'r');
line([Wp2 Wp2], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa2 Wa2], 20*log10([da 10*da]), 'Color', 'r');
line([Wp3 Wp3], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa3 Wa3], 20*log10([da 10*da]), 'Color', 'r');
line([Wp4 Wp4], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa4 Wa4], 20*log10([da 10*da]), 'Color', 'r');
line([Wp5 Wp5], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa5 Wa5], 20*log10([da/10 1+dp]), 'Color', 'r');
line([0 Wp1], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([0 Wp1], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa1 Wa2], 20*log10([da da]), 'Color', 'r');
line([Wp2 Wp3], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([Wp2 Wp3], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa3 Wa4], 20*log10([da da]), 'Color', 'r');
line([Wp4 Wp5], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([Wp4 Wp5], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa5 pi], 20*log10([da da]), 'Color', 'r');
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

audio_filtered1 = filter(h, 1, audio);

[s,f,t] = spectrogram(audio_filtered1, window, overlap, nfft, fs);
s1_dB = 20*log10(abs(s));

figure
imagesc(t, f, s1_dB);
title('FIR filter: optimization method -> spectrogram');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig2 = gca; fig2.XAxis.TickLabelFormat = '%.1f s';
fig2.XMinorTick = 'on';
fig2.YAxis.TickLabelFormat = '%d Hz';
fig2.YMinorTick = 'on';
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
line([Wp1 Wp1], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa1 Wa1], 20*log10([da 10*da]), 'Color', 'r');
line([Wp2 Wp2], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa2 Wa2], 20*log10([da 10*da]), 'Color', 'r');
line([Wp3 Wp3], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa3 Wa3], 20*log10([da 10*da]), 'Color', 'r');
line([Wp4 Wp4], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa4 Wa4], 20*log10([da 10*da]), 'Color', 'r');
line([Wp5 Wp5], 20*log10([da/10 1+dp]), 'Color', 'r');
line([Wa5 Wa5], 20*log10([da/10 1+dp]), 'Color', 'r');
line([0 Wp1], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([0 Wp1], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa1 Wa2], 20*log10([da da]), 'Color', 'r');
line([Wp2 Wp3], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([Wp2 Wp3], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa3 Wa4], 20*log10([da da]), 'Color', 'r');
line([Wp4 Wp5], 20*log10([1-dp 1-dp]), 'Color', 'r');
line([Wp4 Wp5], 20*log10([1+dp 1+dp]), 'Color', 'r');
line([Wa5 pi], 20*log10([da da]), 'Color', 'r');
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

audio_filtered2 = filter(h, 1, audio);

[s,f,t] = spectrogram(audio_filtered2, window, overlap, nfft, fs);
s1_dB = 20*log10(abs(s));

figure
imagesc(t, f, s1_dB);
title('FIR filter: limited impulse response method -> spectrogram');
axis('xy');
xlabel('Time');
ylabel('Frequency');
fig2 = gca; fig2.XAxis.TickLabelFormat = '%.1f s';
fig2.XMinorTick = 'on';
fig2.YAxis.TickLabelFormat = '%d Hz';
fig2.YMinorTick = 'on';
ylim([0 16000]);

savefig('FIR filter limited impulse response method spectrogram')
print('FIR filter limited impulse response method spectrogram','-dsvg','-r0')


audiowrite('out_signal2_2015_0202.wav', audio_filtered2, fs);
