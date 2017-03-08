clear;clc;close all;
% init
[ Rx_AC,DC ] = dump_read( 'dump_16qam.dmp' ); %在这里改名字
% Rx_AC = Rx_AC + DC;
DATA_CARRIER = 52;
FFT_POINT = 128;
SYMBOL_LEN = 160;
CP_HLEN = 16;
START_POINT = 1;
SPEED = SYMBOL_LEN/4;

% plot
hgcf = figure;
set(hgcf,'color',[0 0 0])

subplot(2,1,1);
time_line = START_POINT:SYMBOL_LEN*20+START_POINT;
DATA_time = Rx_AC(time_line);
time_gcf = plot(time_line,abs(DATA_time),'g','Linewidth',2);
set(gca, 'color', [0 0 0]);
set(gca, 'xcolor', [1 1 1]);
set(gca, 'ycolor', [1 1 1]);
timeline_gcf = line([START_POINT START_POINT],[0 max(abs(Rx_AC))],'color',[1 1 1],'linewidth',2);
xlabel('Time Line');ylabel('Waveform');

subplot(2,1,2);
FFT_N = 128; % 让图好看些
DATA_fft = 1/FFT_POINT*fftshift(fft(Rx_AC(1+START_POINT+CP_HLEN:SYMBOL_LEN-CP_HLEN+START_POINT),FFT_N));
spectrum_max = max(abs(DATA_fft));
DATA_spectrum = 20*log10(abs(DATA_fft)/spectrum_max);
spectrum_line = linspace(-20,20,length(DATA_spectrum));
spectrum_gcf = plot(spectrum_line,DATA_spectrum,'g','Linewidth',2,'Marker','None');
set(gca, 'color', [0 0 0]);
set(gca, 'xcolor', [1 1 1]);
set(gca, 'ycolor', [1 1 1]);
%set(spectrum_gcf,'EraseMode','xor');
xlabel('Spectrum Line(MHz)');ylabel('Waveform');
axis([-15 15 -80 0]);
line([10 10],[-80 0],'color',[1 1 1],'linewidth',2);
line([-10 -10],[-80 0],'color',[1 1 1],'linewidth',2);
% for i = 1:length(Rx_AC)-56
%     plot(abs(fft(Rx_AC(i:i+56),64)));
% end
%循环语句中更新坐标数据，一般使用for或者while
try %显示循环
    for i=START_POINT:SPEED:length(Rx_AC)-FFT_POINT 
        time_line = START_POINT+i:START_POINT+SYMBOL_LEN*20+i;
        DATA_time = Rx_AC(time_line);
        set(time_gcf,'xdata',time_line,'ydata',abs(DATA_time));
        set(timeline_gcf,'xdata',[START_POINT+i START_POINT+i]);
        
        DATA_fft = 1/FFT_POINT*fftshift(fft(Rx_AC(1+i+CP_HLEN:SYMBOL_LEN-CP_HLEN+i),FFT_N));
        spectrum_max = max(spectrum_max,max(abs(DATA_fft)));
        DATA_spectrum = 20*log10(abs(DATA_fft)/spectrum_max);
        
        set(spectrum_gcf,'xdata',spectrum_line,'ydata',DATA_spectrum);%更新图像的坐标数据
        drawnow
        pause(0.02);
    end
catch %关掉图形会出错，正好退出
    fprintf('CLOSED PLOT, CURRENT DATA: %d\n',i);
end