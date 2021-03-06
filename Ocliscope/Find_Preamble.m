clear;clc;close all;
% init
[ Rx_AC,DC ] = dump_read( 'HWTRxDump_20150203_16qam.dmp' ); %在这里改名字
DATA_CARRIER = 32;
FFT_POINT = 128;%32.55ns
% START_POINT = 1;
SPEED = 2;
% 小信号检测阈值
Threshold_Find_Preamble = 1/5;
Preamble_Head = find(abs(Rx_AC)>max(abs(Rx_AC))*Threshold_Find_Preamble,1);
Preamble_Head_ori = Preamble_Head;
Preamble = Rx_AC(Preamble_Head:Preamble_Head+DATA_CARRIER*2-1);





hgcf = figure;
subplot(2,1,1);
time_line = Preamble_Head:DATA_CARRIER*10+Preamble_Head;
Preamble_time = Rx_AC(time_line);
time_gcf = plot(time_line,abs(Preamble_time),'g','Linewidth',1);
set(gca, 'color', [0 0 0]);
set(gca, 'xcolor', [1 1 1]);
set(gca, 'ycolor', [1 1 1]);
timeline_gcf = line([Preamble_Head Preamble_Head],[0 15000],'color',[1 1 1],'linewidth',2);
timeline_gcf2 = line([Preamble_Head+DATA_CARRIER-1 Preamble_Head+DATA_CARRIER-1],[0 15000],'color',[1 1 1],'linewidth',2);
xlabel('Time Line');ylabel('Waveform');

subplot(2,1,2);


FFT_N = 256;

Preamble_fft = 1/FFT_POINT*fftshift(fft(Preamble,FFT_N));
spectrum_max = max(abs(Preamble_fft));
Preamble_spectrum = 20*log10(abs(Preamble_fft)/spectrum_max);

Y_height = -50;
set(hgcf,'color',[0 0 0])
spectrum_line = linspace(-20,20,length(Preamble_spectrum));
spectrum_gcf = plot(spectrum_line,Preamble_spectrum,'g','Linewidth',2,'Marker','None');
set(gca, 'color', [0 0 0]);
set(gca, 'xcolor', [1 1 1]);
set(gca, 'ycolor', [1 1 1]);
% set(spectrum_gcf,'EraseMode','xor');
xlabel('Spectrum Line(MHz)');ylabel('Waveform');
axis([-15 15 Y_height 0]);
line([10 10],[Y_height 0],'color',[1 1 1],'linewidth',2);
line([-10 -10],[Y_height 0],'color',[1 1 1],'linewidth',2);

try %显示循环
    while(Preamble_Head-Preamble_Head_ori<DATA_CARRIER*9)
        Preamble_Head = Preamble_Head + SPEED;
        
        time_line = Preamble_Head:DATA_CARRIER*10+Preamble_Head;
        Preamble_time = Rx_AC(time_line);
        set(time_gcf,'xdata',time_line,'ydata',abs(Preamble_time));
        set(timeline_gcf,'xdata',[Preamble_Head Preamble_Head]);
        set(timeline_gcf2,'xdata',[Preamble_Head+DATA_CARRIER-1 Preamble_Head+DATA_CARRIER-1]);
        
        Preamble = Rx_AC(Preamble_Head:Preamble_Head+DATA_CARRIER*2-1);
        Preamble_fft = 1/FFT_POINT*fftshift(fft(Preamble,FFT_N));
        spectrum_max = max(spectrum_max,max(abs(Preamble_fft)));
        Preamble_spectrum = 20*log10(abs(Preamble_fft)/spectrum_max);
        
        set(spectrum_gcf,'xdata',spectrum_line,'ydata',Preamble_spectrum);%更新图像的坐标数据
        drawnow
        pause(0.05);
    end
catch %关掉图形会出错，正好退出
    fprintf('CLOSED PLOT, CURRENT DATA: %d\n',Preamble_Head);
end