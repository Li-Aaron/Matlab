%Nofdm = 128  256 512 1024 2048 2048
%Nsub  = 72   180 300 600  900  1200
%BW    : 1.4  3   5   10   15   20   MHz
%RealBW: 1.08 2.7 4.5 9    13.5 18   MHz
% TD-LTE Signal --> Random Data
% 1 Frame

clear;clc;close all;

%% INIT
load PrimSync
fig = 1;
BW = 20;
Time_GI = [5.208333333333333e-3 4.687500000000000e-3];%ms
Subcarrier_Spacing = 15; %KHz
switch BW
    case 1.4
        Nfft = 128;
        Nsub = 72;
    case 3        
        Nfft = 256;
        Nsub = 180;
    case 5        
        Nfft = 512;
        Nsub = 300;
    case 10        
        Nfft = 1024;
        Nsub = 600;
    case 15      
        Nfft = 2048;
        Nsub = 900;
    case 20      
        Nfft = 2048;
        Nsub = 1200;
    otherwise
        error('No Such BandWidth, Please Try 1.4/3/5/10/15/20 MHz')
end
Ns = 14; %Number of Subcarriers ( MUST BE EVAL!! )
Nprim = length(PrimSync); %Length of PSS (62 default)
Location = 7; %Location of PSS
SNR = 5;
% CP Length and Index
CPlength = round(Nfft*Subcarrier_Spacing*Time_GI);
CPIdx{1} = [Nfft-CPlength(1)+1:Nfft 1:Nfft];
CPIdx{2} = [Nfft-CPlength(2)+1:Nfft 1:Nfft];
% PSS with GI or All OFDM Symbol = PSS
GI = 1;
if GI
    PrimGILength = 5;%PSS GI for Each Side
else
    PrimGILength = (Nsub-Nprim)/2; %Fill All Blanks In Symbol
end
    
    
QpskSignalTemp = randi(4,Nsub*Ns-Nprim-PrimGILength*2,1); %Generated QPSK Signal (Easy Mode)
QpskSignal = exp(-1i*pi*(QpskSignalTemp-0.5)/2);


%% PRIMARY SYNC SIGNAL (925K in Real Channel)
figure(fig);fig=fig+1;
subplot(2,1,1)
plot(real(PrimSync));
title('PRIMARY SYNC SIGNAL (I)')
subplot(2,1,2)
plot(imag(PrimSync));
title('PRIMARY SYNC SIGNAL (Q)')

%% Direct Insert
% AddedPrimSyncSignal = [QpskSignal(1:end/2);PrimSync;QpskSignal(end/2+1:end)];
% CorredSignal = xcorr(PrimSync,AddedPrimSyncSignal);
% [~,location_max] = max(abs(CorredSignal));
% [~,location_min] = min(abs(CorredSignal));
% 
% 
% figure(fig);fig=fig+1;
% subplot(3,1,1)
% stem(abs(CorredSignal),'Marker','None');
% axis([0 Nsub*(Ns-1)+Nprim 0 max(abs(CorredSignal))*1.1])
% title('CORRED SIGNAL')
% subplot(3,1,2)
% plot(real(AddedPrimSyncSignal));
% axis([0 Nsub*(Ns-1)+Nprim -1 1])
% title('WAVE ADDED PRIMARY SYNC SIGNAL (I)')
% subplot(3,1,3)
% axis([0 Nprim -1 1])
% plot(real(PrimSync));
% title('PRIMARY SYNC SIGNAL (I)')

%% TX
IFFT_PrimSync = Nfft*ifft([zeros(1,1);PrimSync(end/2+1:end);zeros(Nfft-Nprim-1,1);PrimSync(1:end/2)],Nfft);
IFFT_PrimSync_Save(:,1) = real(IFFT_PrimSync);
IFFT_PrimSync_Save(:,2) = imag(IFFT_PrimSync);
IFFT_PrimSync_Save = round(IFFT_PrimSync_Save*6);
%stem(real(IFFT_PrimSync));

figure(fig);fig=fig+1;
subplot(3,1,1)
plot(1/Nfft*abs(fftshift(fft(IFFT_PrimSync,Nfft))));
title('IFFTED PRIMARY SYNC SPECTRUM (Frequency Domain)')
subplot(3,1,2)
plot(real(IFFT_PrimSync));
title('IFFTED PRIMARY SYNC SIGNAL (I-Time Domain)')
subplot(3,1,3)
plot(imag(IFFT_PrimSync));
title('IFFTED PRIMARY SYNC SIGNAL (Q-Time Domain)')


PreIFFT_Signal_Temp = zeros(Nsub,Ns);
PrimMapping = false(Nsub,Ns); %PSS MAP
PrimMapping((Nsub-Nprim)/2+1:Nsub-(Nsub-Nprim)/2,Location) = 1;
PrimGIMapping = false(Nsub,Ns); %PSS GI MAP
PrimGIMapping([(Nsub-Nprim)/2-PrimGILength+1:(Nsub-Nprim)/2,...
    Nsub-(Nsub-Nprim)/2+1:Nsub-(Nsub-Nprim)/2+PrimGILength],Location) = 1;
DataIdx = false(Nsub,Ns); %DATA MAP
%DataIdx(:,[1:Location-1,Location+1:Ns])=1; %Not Used Anymore
DataIdx(~(PrimMapping+PrimGIMapping))=1;

PreIFFT_Signal_Temp(PrimMapping) = PrimSync; %ADD PSS
PreIFFT_Signal_Temp(DataIdx) = QpskSignal; %ADD DATA(SSS/REF SIGNAL be instead by DATA)
%dc
% ©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´
% ©¦FRAME(DATA)        ©¦FRAME(PSS+[DATA])   ©¦
% ©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È
% ©¦DC                 ©¦DC                  ©¦
% ©¦SIGNAL(LAST HALF)  ©¦PSS(LAST HALF)      ©¦
% ©¦EMPTY              ©¦[SIGNAL(LAST HALF)] ©¦
% ©¦                   ©¦EMPTY               ©¦
% ©¦SIGNAL(FIRST HALF) ©¦[SIGNAL(FIRST HALF)]©¦
% ©¦                   ©¦PSS(FIRST HALF)     ©¦
% ©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼
PreIFFT_Signal = [zeros(1,Ns);PreIFFT_Signal_Temp(end/2+1:end,:);zeros(Nfft-Nsub-1,Ns);PreIFFT_Signal_Temp(1:end/2,:)];
IFFT_Signal = Nfft*ifft(PreIFFT_Signal,Nfft);
%Tx_Signal = reshape(IFFT_Signal,[],1);

% ©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´
% ©¦OFDMSYMBOL(1)->CP(1)        ©¦OFDMSYMBOL(2:Ns/2)->CP(2)        ©¦
% ©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È
% ©¦OFDMSYMBOL(Ns/2+1)->CP(1)   ©¦OFDMSYMBOL(Ns/2+2:end)->CP(2)    ©¦
% ©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼
Tx_Signal = [        IFFT_Signal(CPIdx{1},1);...
            reshape(IFFT_Signal(CPIdx{2},2:Ns/2),[],1);...
            IFFT_Signal(CPIdx{1},Ns/2+1);...
            reshape(IFFT_Signal(CPIdx{2},Ns/2+2:end),[],1) ];

%% BB SELF TEST
% xcorr
[Corred_TX_Signal,idx] = xcorr(Tx_Signal,IFFT_PrimSync);
shifted_idx=find(idx==0);


% CP increase the length of all symbols (From Bottom to Top)
if Location <= Ns/2
    AddCPlen = CPlength(1)*1 + CPlength(2) * (Location-1);
else
    AddCPlen = CPlength(1)*2 + CPlength(2) * (Location-2);
end


% PLOT
[PLOTVECTOR,location_max_ifft] = max(abs(Corred_TX_Signal));
[~,location_min_ifft] = min(abs(Corred_TX_Signal));
figure(fig);fig=fig+1;
subplot(3,1,1)
stem(abs(Corred_TX_Signal),'Marker','None');
axis([shifted_idx length(Tx_Signal)+shifted_idx 0 max(abs(Corred_TX_Signal))*1.1])
title('CORRED SIGNAL')
subplot(3,1,2)
plot(real(Tx_Signal));
axis([0 length(Tx_Signal) -max(real(Tx_Signal)) max(real(Tx_Signal))])
title([num2str(BW),'M TD-LTE WAVE (I)'])
subplot(3,1,3)
axis([0 Nfft -1 1])
plot(real(IFFT_PrimSync));
title('LOCAL IFFTED PRIMARY SYNC SIGNAL (I)')


if abs(location_max_ifft-shifted_idx - (Nfft*(Location-1) + AddCPlen) )<3 % 3*35ns=105ns(Offset)
    fprintf('SUCCESS Find PSS\n');
    
    subplot(3,1,1)
    hold on
    line([location_max_ifft location_max_ifft],[PLOTVECTOR PLOTVECTOR],'Color',[1 0 0],'LineWidth',3,'Marker','o');
    hold off
    
    subplot(3,1,2)
    hold on
    line([Nfft*(Location-1)+AddCPlen Nfft*(Location-1)+AddCPlen],...
        [-max(real(Tx_Signal)) max(real(Tx_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    line([Nfft*Location+AddCPlen Nfft*Location+AddCPlen],...
        [-max(real(Tx_Signal)) max(real(Tx_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off
else 
    fprintf('FAIL Find PSS\n');
    fprintf('Frame Detect in Point %d, which is not Point %d\n',location_max_ifft-shifted_idx,Nfft*(Location-1)+AddCPlen);
end


%% AWGN Channel
Channel_Signal = awgn(Tx_Signal,SNR,'measured');
[Corred_Channel_Signal,idx] = xcorr(Channel_Signal,IFFT_PrimSync);
shifted_idx=find(idx==0);

% PLOT
[PLOTVECTOR,location_max_ifft] = max(abs(Corred_Channel_Signal));
[~,location_min_ifft] = min(abs(Corred_Channel_Signal));
figure(fig);fig=fig+1;
subplot(3,1,1)
stem(abs(Corred_Channel_Signal),'Marker','None');
axis([shifted_idx length(Channel_Signal)+shifted_idx 0 max(abs(Corred_Channel_Signal))*1.1])
title('CORRED SIGNAL')
subplot(3,1,2)
plot(real(Channel_Signal));
axis([0 length(Channel_Signal) -max(real(Channel_Signal)) max(real(Channel_Signal))])
title([num2str(BW),'M TD-LTE WAVE PASSED AWGN CHANNEL (I)'])
subplot(3,1,3)
axis([0 Nfft -1 1])
plot(real(IFFT_PrimSync));
title('LOCAL IFFTED PRIMARY SYNC SIGNAL (I)')

if abs(location_max_ifft-shifted_idx - (Nfft*(Location-1) + AddCPlen) )<5 % 5*35ns=175ns(Offset)
    fprintf('SUCCESS Find PSS\n');
    
    subplot(3,1,1)
    hold on
    line([location_max_ifft location_max_ifft],[PLOTVECTOR PLOTVECTOR],'Color',[1 0 0],'LineWidth',3,'Marker','o');
    hold off
    
    subplot(3,1,2)
    hold on
    line([Nfft*(Location-1)+AddCPlen Nfft*(Location-1)+AddCPlen],...
        [-max(real(Channel_Signal)) max(real(Channel_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    line([Nfft*Location+AddCPlen Nfft*Location+AddCPlen],...
        [-max(real(Channel_Signal)) max(real(Channel_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off
else 
    fprintf('FAIL Find PSS\n');
    fprintf('Frame Detect in Point %d, which is not Point %d\n',location_max_ifft-shifted_idx,Nfft*(Location-1)+AddCPlen);
end

