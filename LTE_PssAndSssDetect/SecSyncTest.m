%Nofdm = 128  256 512 1024 2048 2048
%Nsub  = 72   180 300 600  900  1200
%BW    : 1.4  3   5   10   15   20   MHz
%RealBW: 1.08 2.7 4.5 9    13.5 18   MHz
% TD-LTE Signal --> Random Data
% 10 Frames

clear;clc;close all;

%% INIT
load PrimSync
load SecSync
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
Time_Sample = 1/Nfft/Subcarrier_Spacing;%ms
Nframe = 10; %Number of Frames ( MUST BE 10!! )
Ns = 14; %Number of Subcarriers ( MUST BE EVAL!! )
Nprim = length(PrimSync); %Length of PSS (62 default)
Nsec = length(SecSync); %Length of SSS (62 default)
PrimLoc = 3; %Location of PSS
PrimSubLoc = [1,6]; %Subframe Location of PSS (ONLY 2 VALUE)
SecLoc = 14; %Location of SSS
SecSubLoc = [0,5]; %Subframe Location of SSS (ONLY 2 VALUE)
SNR = 5;
% CP Length and Index
CPlength = round(Nfft*Subcarrier_Spacing*Time_GI);
CPIdx{1} = [Nfft-CPlength(1)+1:Nfft 1:Nfft];
CPIdx{2} = [Nfft-CPlength(2)+1:Nfft 1:Nfft];
% PSS with GI or All OFDM Symbol = PSS
GI = 1;
if GI
    PrimGILength = 5;%PSS GI for Each Side
    SecGILength = 5;%SSS GI for Each Side
else
    PrimGILength = (Nsub-Nprim)/2; %Fill All Blanks In Symbol
    SecGILength = (Nsub-Nsec)/2;
end
%% BPSK Signal With PSS/SSS fill banks
BpskSignal = cell(Nframe,1);
for SubF = 0:Nframe-1
    if any(PrimSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        BpskSignalTemp = randi(2,Nsub*Ns-Nprim-PrimGILength*2,1); 
    elseif any(SecSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        BpskSignalTemp = randi(2,Nsub*Ns-Nsec-SecGILength*2,1); 
    elseif any(intersect(PrimSubLoc,SecSubLoc) == SubF) % included in both Locations
        BpskSignalTemp = randi(2,Nsub*Ns-Nprim-Nsec-PrimGILength*2-SecGILength*2,1); 
    else
        BpskSignalTemp = randi(2,Nsub*Ns,1); 
    end
    BpskSignal{SubF+1} = (BpskSignalTemp-1.5)*2; %Generated BPSK Signal (Easy Mode)
end
clear BpskSignalTemp;
%% QPSK Signal With PSS/SSS fill banks
QpskSignal = cell(Nframe,1);
for SubF = 0:Nframe-1
    if any(PrimSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        QpskSignalTemp = randi(4,Nsub*Ns-Nprim-PrimGILength*2,1); %Generated QPSK Signal (Easy Mode)
    elseif any(SecSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        QpskSignalTemp = randi(4,Nsub*Ns-Nsec-SecGILength*2,1); %Generated QPSK Signal (Easy Mode)
    elseif any(intersect(PrimSubLoc,SecSubLoc) == SubF) % included in both Locations
        QpskSignalTemp = randi(4,Nsub*Ns-Nprim-Nsec-PrimGILength*2-SecGILength*2,1); %Generated QPSK Signal (Easy Mode)
    else
        QpskSignalTemp = randi(4,Nsub*Ns,1); %Generated QPSK Signal (Easy Mode)
    end
    QpskSignal{SubF+1} = exp(-1i*pi*(QpskSignalTemp-0.5)/2);
end
clear QpskSignalTemp;

%% DATA SIGNAL
% DataSignal = BpskSignal;
DataSignal = QpskSignal;
%% SECOND SYNC SIGNAL (925K in Real Channel)
figure(fig);fig=fig+1;
subplot(2,1,1)
stem(SecSync);
axis([0 Nsec+1 -1.5 1.5]);
title('SECOND SYNC SIGNAL (1)')
subplot(2,1,2)
stem(SecSync2);
axis([0 Nsec+1 -1.5 1.5]);
title('SECOND SYNC SIGNAL (2)')

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

%% IFFTED SYNC SIGNALS
IFFT_PrimSync = Nfft*ifft([zeros(1,1);PrimSync(end/2+1:end);zeros(Nfft-Nprim-1,1);PrimSync(1:end/2)],Nfft);
IFFT_SecSync = Nfft*ifft([zeros(1,1);SecSync(end/2+1:end);zeros(Nfft-Nsec-1,1);SecSync(1:end/2)],Nfft);
IFFT_SecSync2 = Nfft*ifft([zeros(1,1);SecSync2(end/2+1:end);zeros(Nfft-Nsec-1,1);SecSync2(1:end/2)],Nfft);
IFFT_PrimSync_Save(:,1) = real(IFFT_PrimSync);
IFFT_PrimSync_Save(:,2) = imag(IFFT_PrimSync);
IFFT_PrimSync_Save = round(IFFT_PrimSync_Save*6);
IFFT_SecSync_Save(:,1) = real(IFFT_SecSync(1:16:end));
IFFT_SecSync_Save(:,2) = imag(IFFT_SecSync(1:16:end));
IFFT_SecSync_Save = round(IFFT_SecSync_Save*6);
%stem(real(IFFT_PrimSync));

figure(fig);fig=fig+1;
subplot(3,1,1)
plot(1/Nfft*abs(fftshift(fft(IFFT_SecSync,Nfft))));
title('IFFTED SECOND SYNC SPECTRUM (Frequency Domain(1))')
subplot(3,1,2)
plot(real(IFFT_SecSync));
title('IFFTED SECOND SYNC SIGNAL (1)')
subplot(3,1,3)
plot(imag(IFFT_SecSync2));
title('IFFTED SECOND SYNC SIGNAL (2)')

%% GENMAP
PrimMapping = false(Nsub,Ns); %PSS MAP
PrimMapping((Nsub-Nprim)/2+1:Nsub-(Nsub-Nprim)/2,PrimLoc) = 1;
PrimGIMapping = false(Nsub,Ns); %PSS GI MAP
PrimGIMapping([(Nsub-Nprim)/2-PrimGILength+1:(Nsub-Nprim)/2,...
    Nsub-(Nsub-Nprim)/2+1:Nsub-(Nsub-Nprim)/2+PrimGILength],PrimLoc) = 1;
SecMapping = false(Nsub,Ns); %SSS MAP
SecMapping((Nsub-Nsec)/2+1:Nsub-(Nsub-Nsec)/2,SecLoc) = 1;
SecGIMapping = false(Nsub,Ns); %SSS GI MAP
SecGIMapping([(Nsub-Nsec)/2-SecGILength+1:(Nsub-Nsec)/2,...
    Nsub-(Nsub-Nsec)/2+1:Nsub-(Nsub-Nsec)/2+PrimGILength],SecLoc) = 1;

%% TX
PreIFFT_Signal_Temp = zeros(Nsub,Ns);
IFFT_Signal = cell(Nframe,1);
Tx_Signal = [];
for SubF = 0:Nframe-1
    DataIdx = false(Nsub,Ns); %DATA MAP
    %DataIdx(:,[1:Location-1,Location+1:Ns])=1; %Not Used Anymore
    if any(PrimSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        DataIdx(~(PrimMapping+PrimGIMapping))=1; % + PSS
        PreIFFT_Signal_Temp(PrimMapping) = PrimSync; %ADD PSS
        PreIFFT_Signal_Temp(DataIdx) = DataSignal{SubF+1}; %ADD DATA(REF SIGNAL be instead by DATA)
    elseif any(SecSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
        DataIdx(~(SecMapping+SecGIMapping))=1; % + SSS
        if SubF <= Nframe/2-1
            PreIFFT_Signal_Temp(SecMapping) = SecSync; %ADD PSS
        else
            PreIFFT_Signal_Temp(SecMapping) = SecSync2; %ADD PSS
        end
        PreIFFT_Signal_Temp(DataIdx) = DataSignal{SubF+1}; %ADD DATA(REF SIGNAL be instead by DATA)        
    elseif any(intersect(PrimSubLoc,SecSubLoc) == SubF) % included in both Locations
        DataIdx(~(PrimMapping+PrimGIMapping+SecMapping+SecGIMapping))=1; % + PSS / SSS
        PreIFFT_Signal_Temp(PrimMapping) = PrimSync; %ADD PSS
        PreIFFT_Signal_Temp(SecMapping) = SecSync; %ADD PSS
        PreIFFT_Signal_Temp(DataIdx) = DataSignal{SubF+1}; %ADD DATA(REF SIGNAL be instead by DATA)       
    else
        DataIdx(:)=1; % ALL Data
        PreIFFT_Signal_Temp(DataIdx) = DataSignal{SubF+1}; %ADD DATA(REF SIGNAL be instead by DATA)     
    end



%dc
% ©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´
% ©¦SYMBOL(DATA)       ©¦SYMBOL(PSS+[DATA])  ©¦SYMBOL(SSS+[DATA])  ©¦
% ©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È
% ©¦DC                 ©¦DC                  ©¦DC                  ©¦
% ©¦SIGNAL(LAST HALF)  ©¦PSS(LAST HALF)      ©¦SSS(LAST HALF)      ©¦
% ©¦EMPTY              ©¦[SIGNAL(LAST HALF)] ©¦[SIGNAL(LAST HALF)] ©¦
% ©¦                   ©¦EMPTY               ©¦EMPTY               ©¦
% ©¦SIGNAL(FIRST HALF) ©¦[SIGNAL(FIRST HALF)]©¦[SIGNAL(FIRST HALF)]©¦
% ©¦                   ©¦PSS(FIRST HALF)     ©¦SSS(FIRST HALF)     ©¦
% ©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼
PreIFFT_Signal = [zeros(1,Ns);PreIFFT_Signal_Temp(end/2+1:end,:);zeros(Nfft-Nsub-1,Ns);PreIFFT_Signal_Temp(1:end/2,:)];
IFFT_Signal{SubF+1} = Nfft*ifft(PreIFFT_Signal,Nfft);
%Tx_Signal = reshape(IFFT_Signal,[],1);

% ©°©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ð©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©´
% ©¦OFDMSYMBOL(1)->CP(1)        ©¦OFDMSYMBOL(2:Ns/2)->CP(2)        ©¦
% ©À©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©à©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©È
% ©¦OFDMSYMBOL(Ns/2+1)->CP(1)   ©¦OFDMSYMBOL(Ns/2+2:end)->CP(2)    ©¦
% ©¸©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©Ø©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¼
Tx_Signal = [Tx_Signal;
            IFFT_Signal{SubF+1}(CPIdx{1},1);...
            reshape(IFFT_Signal{SubF+1}(CPIdx{2},2:Ns/2),[],1);...
            IFFT_Signal{SubF+1}(CPIdx{1},Ns/2+1);...
            reshape(IFFT_Signal{SubF+1}(CPIdx{2},Ns/2+2:end),[],1) ];
end
x_timeline = linspace(0,Time_Sample*(length(Tx_Signal)-1),length(Tx_Signal));

%PLOT
figure(fig);fig=fig+1;
subplot(2,1,1);
plot(x_timeline,real(Tx_Signal));
title('TD-LTE WAVEFORM (I)');
subplot(2,1,2);
plot(x_timeline,imag(Tx_Signal));
title('TD-LTE WAVEFORM (Q)');
%% BB SELF TEST PSS
% xcorr
[Corred_TX_Signal_Prim,idx] = xcorr(Tx_Signal,IFFT_PrimSync);
shifted_idx=find(idx==0);


% CP increase the length of all symbols (From Bottom to Top)
if PrimLoc <= Ns/2
    AddCPlenPrim = CPlength(1)*1 + CPlength(2) * (PrimLoc-1);
else
    AddCPlenPrim = CPlength(1)*2 + CPlength(2) * (PrimLoc-2);
end
AllCPlen = CPlength(1)*2 + CPlength(2) * (Ns-2);

% PLOT
[PLOTVECTOR,location_max_ifft] = max(abs(Corred_TX_Signal_Prim));
%[~,location_min_ifft] = min(abs(Corred_TX_Signal));
figure(fig);fig=fig+1;
subplot(3,1,1)
stem(abs(Corred_TX_Signal_Prim),'Marker','None');
axis([shifted_idx length(Tx_Signal)+shifted_idx 0 max(abs(Corred_TX_Signal_Prim))*1.1])
title('CORRED SIGNAL')
subplot(3,1,2)
plot(x_timeline,real(Tx_Signal));
axis([0 ceil(max(x_timeline)) -max(real(Tx_Signal)) max(real(Tx_Signal))])
title([num2str(BW),'M TD-LTE WAVE (I)'])
subplot(3,1,3)
axis([0 Nfft -1 1])
plot(real(IFFT_PrimSync));
title('LOCAL IFFTED PRIMARY SYNC SIGNAL (I)')

Loc_Prim = (Nfft*(PrimLoc-1) + AddCPlenPrim) + (Nfft*Ns+AllCPlen).*PrimSubLoc;
Loc_maxcorr = location_max_ifft-shifted_idx;
if any(abs(Loc_Prim - Loc_maxcorr)<3) % 3*35ns=105ns(Offset)
    fprintf('SUCCESS Find PSS\n');
    
    % PLOT
    subplot(3,1,1)
    hold on
    line([location_max_ifft location_max_ifft],[PLOTVECTOR PLOTVECTOR],'Color',[1 0 0],'LineWidth',3,'Marker','o');
    hold off
    
    subplot(3,1,2)
    hold on
    PrimFound = find(abs(Loc_Prim - Loc_maxcorr)<3);
    line([Time_Sample*Loc_Prim(PrimFound) Time_Sample*Loc_Prim(PrimFound)],...
        [-max(real(Tx_Signal)) max(real(Tx_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off
else 
    PrimFound = 0;
    fprintf('FAIL Find PSS\n');
    fprintf('Frame Detect in Point %d, which is not Point %d or %d\n',Loc_maxcorr,Loc_Prim(1),Loc_Prim(2));
end

%% BB SELF TEST SSS
if PrimFound
    % xcorr
    [Corred_TX_Signal_Sec1,idx1] = xcorr(Tx_Signal,IFFT_SecSync);
    [Corred_TX_Signal_Sec2,idx2] = xcorr(Tx_Signal,IFFT_SecSync2);
    shifted_idx1=find(idx1==0);
    shifted_idx2=find(idx2==0);
    % CP increase the length of all symbols (From Bottom to Top)
    if SecLoc <= Ns/2
        AddCPlenSec = CPlength(1)*1 + CPlength(2) * (SecLoc-1);
    else
        AddCPlenSec = CPlength(1)*2 + CPlength(2) * (SecLoc-2);
    end

    % PLOT
    [PLOTVECTOR(1),location_max_sec(1)] = max(abs(Corred_TX_Signal_Sec1));
    [PLOTVECTOR(2),location_max_sec(2)] = max(abs(Corred_TX_Signal_Sec2));
    figure(fig);fig=fig+1;
    subplot(3,1,1)
    stem(abs(Corred_TX_Signal_Sec1),'Marker','None');
    axis([shifted_idx1 length(Tx_Signal)+shifted_idx1 0 max(abs(Corred_TX_Signal_Sec1))*1.1])
    title('CORRED SIGNAL WITH SSS (1)')
    subplot(3,1,2)
    stem(abs(Corred_TX_Signal_Sec2),'Marker','None');
    axis([shifted_idx2 length(Tx_Signal)+shifted_idx2 0 max(abs(Corred_TX_Signal_Sec2))*1.1])
    title('CORRED SIGNAL WITH SSS (2)')
    subplot(3,1,3)
    plot(x_timeline,real(Tx_Signal));
    axis([0 ceil(max(x_timeline)) -max(real(Tx_Signal)) max(real(Tx_Signal))])
    title([num2str(BW),'M TD-LTE WAVE (I)'])
    hold on
    line([Time_Sample*Loc_Prim(PrimFound) Time_Sample*Loc_Prim(PrimFound)],...
        [-max(real(Tx_Signal)) max(real(Tx_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off

    Loc_Sec = (Nfft*(SecLoc-1) + AddCPlenSec) + (Nfft*Ns+AllCPlen).*SecSubLoc;
    Loc_maxcorr_SSS(1) = location_max_sec(1)-shifted_idx1;
    Loc_maxcorr_SSS(2) = location_max_sec(2)-shifted_idx2;
    % SSS
    for Secidx = [1 2]
        if abs(Loc_Sec(Secidx) - Loc_maxcorr_SSS(Secidx))<3 % 3*35ns=105ns(Offset)
            fprintf('SUCCESS Find SSS_1\n');

            % PLOT
            subplot(3,1,Secidx)
            hold on
            line([location_max_sec(Secidx) location_max_sec(Secidx)],[PLOTVECTOR(Secidx) PLOTVECTOR(Secidx)],...
                'Color',[0.5 0.5 0],'LineWidth',3,'Marker','o');
            hold off
            if PrimFound == Secidx
                subplot(3,1,3)
                hold on
                SecFound = Secidx;
                line([Time_Sample*Loc_Sec(SecFound) Time_Sample*Loc_Sec(SecFound)],...
                    [-max(real(Tx_Signal)) max(real(Tx_Signal))],...
                    'Color',[0.5 0.5 0],'LineWidth',3);
                hold off
            end
        else
            fprintf('FAIL Find SSS_%d\n',Secidx);
            fprintf('Frame Detect in Point %d, which is not Point %d\n',Loc_maxcorr_SSS(1),Loc_Sec(1));
        end
    end
end
%% AWGN Channel PSS
Channel_Signal = awgn(Tx_Signal,SNR,'measured');

% xcorr
[Corred_Channel_Signal,idx] = xcorr(Channel_Signal,IFFT_PrimSync);
shifted_idx=find(idx==0);

% PLOT
[PLOTVECTOR,location_max_ifft] = max(abs(Corred_Channel_Signal));
%[~,location_min_ifft] = min(abs(Corred_TX_Signal));
figure(fig);fig=fig+1;
subplot(3,1,1)
stem(abs(Corred_Channel_Signal),'Marker','None');
axis([shifted_idx length(Channel_Signal)+shifted_idx 0 max(abs(Corred_Channel_Signal))*1.1])
title('CORRED SIGNAL')
subplot(3,1,2)
plot(x_timeline,real(Channel_Signal));
axis([0 ceil(max(x_timeline)) -max(real(Channel_Signal)) max(real(Channel_Signal))])
title([num2str(BW),'M TD-LTE WAVE OVER AWGN CHANNEL(I) SNR =',num2str(SNR),'(dB)'])
subplot(3,1,3)
axis([0 Nfft -1 1])
plot(real(IFFT_PrimSync));
title('LOCAL IFFTED PRIMARY SYNC SIGNAL (I)')

Loc_Prim = (Nfft*(PrimLoc-1) + AddCPlenPrim) + (Nfft*Ns+AllCPlen).*PrimSubLoc;
Loc_maxcorr = location_max_ifft-shifted_idx;
if any(abs(Loc_Prim - Loc_maxcorr)<3) % 3*35ns=105ns(Offset)
    fprintf('SUCCESS Find PSS\n');
    
    % PLOT
    subplot(3,1,1)
    hold on
    line([location_max_ifft location_max_ifft],[PLOTVECTOR PLOTVECTOR],'Color',[1 0 0],'LineWidth',3,'Marker','o');
    hold off
    
    subplot(3,1,2)
    hold on
    PrimFound = find(abs(Loc_Prim - Loc_maxcorr)<3);
    line([Time_Sample*Loc_Prim(PrimFound) Time_Sample*Loc_Prim(PrimFound)],...
        [-max(real(Channel_Signal)) max(real(Channel_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off
else 
    PrimFound = 0;
    fprintf('FAIL Find PSS\n');
    fprintf('Frame Detect in Point %d, which is not Point %d or %d\n',Loc_maxcorr,Loc_Prim(1),Loc_Prim(2));
end

%% AWGN Channel SSS
if PrimFound
    % xcorr
    [Corred_Channel_Signal_Sec1,idx1] = xcorr(Channel_Signal,IFFT_SecSync);
    [Corred_Channel_Signal_Sec2,idx2] = xcorr(Channel_Signal,IFFT_SecSync2);
    shifted_idx1=find(idx1==0);
    shifted_idx2=find(idx2==0);

    % PLOT
    [PLOTVECTOR(1),location_max_sec(1)] = max(abs(Corred_Channel_Signal_Sec1));
    [PLOTVECTOR(2),location_max_sec(2)] = max(abs(Corred_Channel_Signal_Sec2));
    figure(fig);fig=fig+1;
    subplot(3,1,1)
    stem(abs(Corred_Channel_Signal_Sec1),'Marker','None');
    axis([shifted_idx1 length(Channel_Signal)+shifted_idx1 0 max(abs(Corred_Channel_Signal_Sec1))*1.1])
    title('CORRED SIGNAL WITH SSS (1)')
    subplot(3,1,2)
    stem(abs(Corred_Channel_Signal_Sec2),'Marker','None');
    axis([shifted_idx2 length(Channel_Signal)+shifted_idx2 0 max(abs(Corred_Channel_Signal_Sec2))*1.1])
    title('CORRED SIGNAL WITH SSS (2)')
    subplot(3,1,3)
    plot(x_timeline,real(Channel_Signal));
    axis([0 ceil(max(x_timeline)) -max(real(Channel_Signal)) max(real(Channel_Signal))])
    title([num2str(BW),'M TD-LTE WAVE OVER AWGN CHANNEL(I) SNR =',num2str(SNR),'(dB)'])
    hold on
    line([Time_Sample*Loc_Prim(PrimFound) Time_Sample*Loc_Prim(PrimFound)],...
        [-max(real(Channel_Signal)) max(real(Channel_Signal))],...
        'Color',[1 0 0],'LineWidth',3);
    hold off

    Loc_Sec = (Nfft*(SecLoc-1) + AddCPlenSec) + (Nfft*Ns+AllCPlen).*SecSubLoc;
    Loc_maxcorr_SSS(1) = location_max_sec(1)-shifted_idx1;
    Loc_maxcorr_SSS(2) = location_max_sec(2)-shifted_idx2;
    % SSS
    for Secidx = [1 2]
        if abs(Loc_Sec(Secidx) - Loc_maxcorr_SSS(Secidx))<3 % 3*35ns=105ns(Offset)
            fprintf('SUCCESS Find SSS_1\n');

            % PLOT
            subplot(3,1,Secidx)
            hold on
            line([location_max_sec(Secidx) location_max_sec(Secidx)],[PLOTVECTOR(Secidx) PLOTVECTOR(Secidx)],...
                'Color',[0.5 0.5 0],'LineWidth',3,'Marker','o');
            hold off
            if PrimFound == Secidx
                subplot(3,1,3)
                hold on
                SecFound = Secidx;
                line([Time_Sample*Loc_Sec(SecFound) Time_Sample*Loc_Sec(SecFound)],...
                    [-max(real(Channel_Signal)) max(real(Channel_Signal))],...
                    'Color',[0.5 0.5 0],'LineWidth',3);
                hold off
            end
        else
            fprintf('FAIL Find SSS_%d\n',Secidx);
            fprintf('Frame Detect in Point %d, which is not Point %d\n',Loc_maxcorr_SSS(1),Loc_Sec(1));
        end
    end
end