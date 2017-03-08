clear;clc;close all;

%% INIT
load PrimSync
load SecSync
fig = 1;
BW = 5;
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


%% IFFTED SYNC SIGNALS
IFFT_PrimSync = Nfft*ifft([zeros(1,1);PrimSync(end/2+1:end);zeros(Nfft-Nprim-1,1);PrimSync(1:end/2)],Nfft);
IFFT_SecSync = Nfft*ifft([zeros(1,1);SecSync(end/2+1:end);zeros(Nfft-Nsec-1,1);SecSync(1:end/2)],Nfft);
IFFT_SecSync2 = Nfft*ifft([zeros(1,1);SecSync2(end/2+1:end);zeros(Nfft-Nsec-1,1);SecSync2(1:end/2)],Nfft);

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

%% OTHER GENERATED PARAMETERS
% CP increase the length of all symbols (From Bottom to Top)
if PrimLoc <= Ns/2
    AddCPlenPrim = CPlength(1)*1 + CPlength(2) * (PrimLoc-1);
else
    AddCPlenPrim = CPlength(1)*2 + CPlength(2) * (PrimLoc-2);
end
if SecLoc <= Ns/2
    AddCPlenSec = CPlength(1)*1 + CPlength(2) * (SecLoc-1);
else
    AddCPlenSec = CPlength(1)*2 + CPlength(2) * (SecLoc-2);
end
AllCPlen = CPlength(1)*2 + CPlength(2) * (Ns-2);



%% TEST LOOPS START FROM GENERATING THE SIGNALS
SNR_TEST = -20:5;
SNRmin = min(SNR_TEST);
MC_max = 1000;
Success_PSS = zeros(1,length(SNR_TEST));
Success_SSS = zeros(1,length(SNR_TEST));
h = waitbar(0,['CURRENT SNR = ',num2str(SNRmin)]);
for i = 1:length(SNR_TEST)
    waitbar((i-1)/length(SNR_TEST),h,['CURRENT SNR = ',num2str(i+SNRmin-1)]);
    for MC = 1:MC_max
        %% BPSK Signal With PSS/SSS fill banks
%                 BpskSignal = cell(Nframe,1);
%                 for SubF = 0:Nframe-1
%                     if any(PrimSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
%                         BpskSignalTemp = randi(2,Nsub*Ns-Nprim-PrimGILength*2,1);
%                     elseif any(SecSubLoc == SubF) && all(intersect(PrimSubLoc,SecSubLoc) ~= SubF) % not included in both Locations
%                         BpskSignalTemp = randi(2,Nsub*Ns-Nsec-SecGILength*2,1);
%                     elseif any(intersect(PrimSubLoc,SecSubLoc) == SubF) % included in both Locations
%                         BpskSignalTemp = randi(2,Nsub*Ns-Nprim-Nsec-PrimGILength*2-SecGILength*2,1);
%                     else
%                         BpskSignalTemp = randi(2,Nsub*Ns,1);
%                     end
%                     BpskSignal{SubF+1} = (BpskSignalTemp-1.5)*2; %Generated BPSK Signal (Easy Mode)
%                 end
%                 clear BpskSignalTemp;
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
        %x_timeline = linspace(0,Time_Sample*(length(Tx_Signal)-1),length(Tx_Signal));
        
        %% AWGN Channel PSS
        Channel_Signal = awgn(Tx_Signal,SNR_TEST(i),'measured');
        
        % xcorr
        [Corred_Channel_Signal,idx] = xcorr(Channel_Signal,IFFT_PrimSync);
        shifted_idx=find(idx==0);
        
        [PLOTVECTOR,location_max_ifft] = max(abs(Corred_Channel_Signal));
        
        Loc_Prim = (Nfft*(PrimLoc-1) + AddCPlenPrim) + (Nfft*Ns+AllCPlen).*PrimSubLoc;
        Loc_maxcorr = location_max_ifft-shifted_idx;
        if any(abs(Loc_Prim - Loc_maxcorr)<3) % 3*35ns=105ns(Offset)
            PrimFound = find(abs(Loc_Prim - Loc_maxcorr)<3);
        else
            PrimFound = 0;
        end
        
        %% AWGN Channel SSS
        if PrimFound
            % xcorr
            [Corred_Channel_Signal_Sec1,idx1] = xcorr(Channel_Signal,IFFT_SecSync);
            [Corred_Channel_Signal_Sec2,idx2] = xcorr(Channel_Signal,IFFT_SecSync2);
            shifted_idx1=find(idx1==0);
            shifted_idx2=find(idx2==0);
            
            [PLOTVECTOR(1),location_max_sec(1)] = max(abs(Corred_Channel_Signal_Sec1));
            [PLOTVECTOR(2),location_max_sec(2)] = max(abs(Corred_Channel_Signal_Sec2));
            
            Loc_Sec = (Nfft*(SecLoc-1) + AddCPlenSec) + (Nfft*Ns+AllCPlen).*SecSubLoc;
            Loc_maxcorr_SSS(1) = location_max_sec(1)-shifted_idx1;
            Loc_maxcorr_SSS(2) = location_max_sec(2)-shifted_idx2;
            % SSS
            for Secidx = [1 2]
                if abs(Loc_Sec(Secidx) - Loc_maxcorr_SSS(Secidx))<3 % 3*35ns=105ns(Offset)
                    if PrimFound == Secidx
                        SecFound = Secidx;
                    end
                else
                    SecFound = 0;
                end
            end
        end
        
        if PrimFound
            Success_PSS(i) = Success_PSS(i) + 1;
            if SecFound
                Success_SSS(i) = Success_SSS(i) + 1;
            end
        end
        
    end
end
delete(h);
Success_Rate_PSS = Success_PSS/MC_max;
Success_Rate_SSS = Success_SSS/MC_max;
plot(SNR_TEST,Success_Rate_PSS);
hold on
plot(SNR_TEST,Success_Rate_SSS,'r');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Successful Detect Rate');
legend('PSS','SSS');