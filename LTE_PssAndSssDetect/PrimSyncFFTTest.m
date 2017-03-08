clear;clc;close all;

%% INIT
load PrimSync
Nfft = 2048;
Nsub = 1200;
% % RANDOM PSK
% PrimSync_temp = randi(63,62,1)-1;
% PrimSync = exp(-1i*pi*PrimSync_temp/63);
% QPSK
PrimSync_temp = randi(4,Nsub,1);
PrimSync = exp(-1i*pi*(PrimSync_temp-0.5)/2);

Nprim = length(PrimSync); %Length of PSS (62 default)
PrimSync_DC_Last_Zero_First = [zeros(1,1);PrimSync(end/2+1:end);zeros(Nfft-Nprim-1,1);PrimSync(1:end/2)];
%PrimSync_Last_Zero_First = [PrimSync(end/2+1:end);zeros(Nfft-Nprim,1);PrimSync(1:end/2)];
PrimSync_Zero_First_DC_Last_Zero = [zeros((Nfft-Nprim)/2,1);PrimSync(1:end/2);zeros(1,1);PrimSync(end/2+1:end);zeros((Nfft-Nprim)/2-1,1)];
PrimSync_First_DC_Last_Zero = [PrimSync(1:end/2);zeros(1,1);PrimSync(end/2+1:end);zeros(Nfft-Nprim-1,1)];

IFFT_Prim_DLZF = Nfft*ifft(PrimSync_DC_Last_Zero_First,Nfft);
IFFT_Prim_ZFDLZ = Nfft*ifft(PrimSync_Zero_First_DC_Last_Zero,Nfft);
IFFT_Prim_FDLZ = Nfft*ifft(PrimSync_First_DC_Last_Zero,Nfft);
IFFT_Prim_shift_DLZF = Nfft*ifft(fftshift(PrimSync_DC_Last_Zero_First),Nfft);

figure;
%ABS
subplot(3,1,1);
plot(abs(IFFT_Prim_DLZF),'b');
hold on;
plot(abs(IFFT_Prim_ZFDLZ),'r');
plot(abs(IFFT_Prim_FDLZ),'k');
hold off;
title('Waveform ABS');
legend('DLZF','ZFDLZ','FDLZ');

%I
subplot(3,1,2);
plot(real(IFFT_Prim_DLZF),'b');
hold on;
plot(real(IFFT_Prim_ZFDLZ),'r');
plot(real(IFFT_Prim_FDLZ),'k');
hold off;
title('Waveform I');

%Q
subplot(3,1,3);
plot(imag(IFFT_Prim_DLZF),'b');
hold on;
plot(imag(IFFT_Prim_ZFDLZ),'r');
plot(imag(IFFT_Prim_FDLZ),'k');
hold off;
title('Waveform Q');

%% Original Signal Transfered
figure;
subplot(4,2,1);
stem(real(PrimSync_DC_Last_Zero_First),'x','MarkerSize',3);
title('DC Last Zero First(SORA/TD-LTE)');
subplot(4,2,3);
stem(real(PrimSync_Zero_First_DC_Last_Zero),'x','MarkerSize',3);
title('Zero First DC Last Zero(reversed TD-LTE)');
subplot(4,2,5);
stem(real(PrimSync_First_DC_Last_Zero),'x','MarkerSize',3);
title('First DC Last Zero(default IFFT)');
subplot(4,2,7);
stem(real(fftshift(PrimSync_DC_Last_Zero_First)),'x','MarkerSize',3);
title('fftshift with First DC Last Zero(Shifted TD-LTE)');



%% IFFTED
%figure;
subplot(4,2,2);
plot(real(IFFT_Prim_DLZF));
title('IFFTED DLZF(SORA/TD-LTE)');
subplot(4,2,4);
plot(real(IFFT_Prim_ZFDLZ));
title('IFFTED ZFDLZ(reversed TD-LTE)');
subplot(4,2,6);
plot(real(IFFT_Prim_FDLZ));
title('IFFTED FDLZ(default IFFT)');
subplot(4,2,8);
plot(real(IFFT_Prim_shift_DLZF));
title('IFFTED fftshifted DLZF(Shifted TD-LTE)');






