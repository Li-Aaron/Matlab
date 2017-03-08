clear;clc;close all;

%% INIT
Nfft = 2048;
Nsig = 1200;
% % RANDOM PSK
% PrimSync_temp = randi(63,62,1)-1;
% PrimSync = exp(-1i*pi*PrimSync_temp/63);
% QPSK
Signal_Temp = randi(4,Nsig,1);
QPSK_Signal = exp(-1i*pi*(Signal_Temp-0.5)/2);

NSignal = length(QPSK_Signal); %Length of PSS (62 default)

Signal_DC_Last_Zero_First = [zeros(1,1);QPSK_Signal(end/2+1:end);zeros(Nfft-NSignal-1,1);QPSK_Signal(1:end/2)];
Signal_Zero_First_DC_Last_Zero = [zeros((Nfft-NSignal)/2,1);QPSK_Signal(1:end/2);zeros(1,1);QPSK_Signal(end/2+1:end);zeros((Nfft-NSignal)/2-1,1)];
Signal_First_DC_Last_Zero = [QPSK_Signal(1:end/2);zeros(1,1);QPSK_Signal(end/2+1:end);zeros(Nfft-NSignal-1,1)];


IFFT_Signal_DLZF = ifft(Signal_DC_Last_Zero_First,Nfft);
IFFT_Signal_ZFDLZ = ifft(Signal_Zero_First_DC_Last_Zero,Nfft);
IFFT_Signal_FDLZ = ifft(Signal_First_DC_Last_Zero,Nfft);

FFT_Signal_DLZF = fft(IFFT_Signal_DLZF,Nfft);
FFT_Signal_Shifted_DLZF = fftshift(FFT_Signal_DLZF);
FFT_Signal_ZFDLZ = fft(IFFT_Signal_ZFDLZ,Nfft);
FFT_Signal_Shifted_ZFDLZ = fftshift(FFT_Signal_ZFDLZ);
FFT_Signal_FDLZ = fft(IFFT_Signal_FDLZ,Nfft);
FFT_Signal_Shifted_FDLZ = fftshift(FFT_Signal_FDLZ);

%% Original Signal Transfered
figure;
subplot(3,2,1);
stem(real(Signal_DC_Last_Zero_First),'x','MarkerSize',3);
title('DC Last Zero First(SORA/TD-LTE)');
subplot(3,2,3);
stem(real(Signal_Zero_First_DC_Last_Zero),'x','MarkerSize',3);
title('Zero First DC Last Zero(reversed TD-LTE)');
subplot(3,2,5);
stem(real(Signal_First_DC_Last_Zero),'x','MarkerSize',3);
title('First DC Last Zero(default IFFT)');

%% IFFTED
%figure;
subplot(3,2,2);
plot(real(IFFT_Signal_DLZF));
title('IFFTED DLZF(SORA/TD-LTE)');
subplot(3,2,4);
plot(real(IFFT_Signal_ZFDLZ));
title('IFFTED ZFDLZ(reversed TD-LTE)');
subplot(3,2,6);
plot(real(IFFT_Signal_FDLZ));
title('IFFTED FDLZ(default IFFT)');

%% FFTED
figure;
%DLZF
subplot(3,2,1);
plot(abs(FFT_Signal_DLZF),'b');
title('DC Last Zero First');
subplot(3,2,2);
plot(abs(FFT_Signal_Shifted_DLZF),'b');
title('DC Last Zero First (Shifted)');
%ZFDLZ
subplot(3,2,3);
plot(abs(FFT_Signal_ZFDLZ),'b');
title('Zero First DC Last Zero');
subplot(3,2,4);
plot(abs(FFT_Signal_Shifted_ZFDLZ),'b');
title('Zero First DC Last Zero (Shifted)');
%FDLZ
subplot(3,2,5);
plot(abs(FFT_Signal_FDLZ),'b');
title('First DC Last Zero');
subplot(3,2,6);
plot(abs(FFT_Signal_Shifted_FDLZ),'b');
title('First DC Last Zero (Shifted)');

%% COSINE (Single Frequency)
% all points = Nfft
Period = 128; %points
% if consider all points are within one sec then Num_of_Period can also
% be considered as Frequency 
% -> thus FFT_Offset is the max_specturm_point - zero which can also be
% considered as Frequency
% if these two value is similar, the fft is right.
Num_of_Period = Nfft/Period; 
x = 0:Nfft-1; %2000 points
Cos_Signal = cos(x*2*pi/Period); 
FFT_Cos_Signal = fft(Cos_Signal,Nfft);
FFT_Cos_Signal_Shifted = fftshift(FFT_Cos_Signal);
%PLOT
figure;
subplot(3,1,1);
plot(Cos_Signal);
title('COSINE Signal Time Domin (Single Frequency)');
subplot(3,1,2);
plot(abs(FFT_Cos_Signal));
title('COSINE Spectrum Frequency Domin(Single Frequency)');
subplot(3,1,3);
plot(abs(FFT_Cos_Signal_Shifted));
title('COSINE Spectrum Frequency Domin Shifted(Single Frequency)');

[PLOTVECTOR,idx_of_max_spectrum] = max(abs(fftshift(FFT_Cos_Signal)));
FFT_Offset = abs(idx_of_max_spectrum - Nfft/2);
Offset_direction = sign(idx_of_max_spectrum - Nfft/2);
%PLOT
subplot(3,1,3);
hold on    
line([idx_of_max_spectrum idx_of_max_spectrum],[0 PLOTVECTOR],'Color',[1 0 0],'LineWidth',3);
line([Nfft/2+Offset_direction*Num_of_Period Nfft/2+Offset_direction*Num_of_Period],[0 PLOTVECTOR],'Color',[1 0.5 0.5],'LineWidth',3);
line([Nfft/2 Nfft/2],[0 PLOTVECTOR],'Color',[1 0.8 0.8],'LineWidth',3);
hold off

% then we convert the specturm signal back to time domin
IFFT_Cos_Signal = ifft(FFT_Cos_Signal);
IFFT_Cos_Signal_Shifted = ifft(FFT_Cos_Signal_Shifted,Nfft);

figure;
subplot(3,1,1);
plot(Cos_Signal);
title('COSINE Signal Time Domin (Single Frequency)');
subplot(3,1,2);
plot(IFFT_Cos_Signal);
title('COSINE Signal Time Domin(Recovered)');
subplot(3,1,3);
plot(IFFT_Cos_Signal_Shifted);
title('COSINE Signal Time Domin (Shifted then Recover)');


