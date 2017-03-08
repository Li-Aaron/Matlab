clear;clc;close all;

%% INIT
load PrimSync
fig = 1;
BW = 1.4;
switch BW
    case 1.4
        Nofdm = 128;
        Nsub = 72;
    case 3        
        Nofdm = 256;
        Nsub = 180;
    case 5        
        Nofdm = 512;
        Nsub = 300;
    case 10        
        Nofdm = 1024;
        Nsub = 600;
    case 15      
        Nofdm = 2048;
        Nsub = 900;
    case 20      
        Nofdm = 2048;
        Nsub = 1200;
    otherwise
        error('No Such BandWidth, Please Try 1.4/3/5/10/15/20 MHz')
end
Ns = 14; %Number of Subcarriers
Nprim = length(PrimSync); %Length of PSS (62 default)
Location = 7; %Location of PSS
% PSS with GI or All OFDM Symbol = PSS
GI = 0;
if GI
    PrimGILength = 5;%PSS GI for Each Side
else
    PrimGILength = (Nsub-Nprim)/2; %Fill All Blanks In Symbol
end
    
    

%% AWGN SNR -- FRAME DETECT SUCCESS RATE
SNR_TEST = -20:5;
SNRmin = min(SNR_TEST);
MC_max = 2000;
Success = zeros(1,length(SNR_TEST));
h = waitbar(0,['CURRENT SNR = ',num2str(SNRmin)]);
for i = 1:length(SNR_TEST)
    waitbar((i-1)/length(SNR_TEST),h,['CURRENT SNR = ',num2str(i+SNRmin-1)]);
    for MC = 1:MC_max
        QpskSignalTemp = randi(4,Nsub*Ns-Nprim-PrimGILength*2,1); %Generated QPSK Signal (Easy Mode)
        QpskSignal = exp(-1i*pi*(QpskSignalTemp-0.5)/2);


        IFFT_PrimSync = Nofdm*ifft([zeros(1,1);PrimSync(end/2+1:end);zeros(Nofdm-Nprim-1,1);PrimSync(1:end/2)],Nofdm);
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

        PreIFFT_Signal = [zeros(1,Ns);PreIFFT_Signal_Temp(end/2+1:end,:);zeros(Nofdm-Nsub-1,Ns);PreIFFT_Signal_Temp(1:end/2,:)];
        IFFT_Signal = Nofdm*ifft(PreIFFT_Signal,Nofdm);
        Tx_Signal = reshape(IFFT_Signal,[],1);

        Channel_Signal = awgn(Tx_Signal,SNR_TEST(i),'measured');

        [Corred_Channel_Signal,idx] = xcorr(Channel_Signal,IFFT_PrimSync);
        shifted_idx=find(idx==0);
        [PLOTVECTOR,location_max_ifft] = max(abs(Corred_Channel_Signal));
        if abs(location_max_ifft-shifted_idx - Nofdm*(Location-1))<5 %offset<175ns
            Success(i) = Success(i) + 1;
        end
    end
end
delete(h);
Success_Rate = Success/MC_max;
plot(SNR_TEST,Success_Rate);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('PSS Successful Detect Rate');