function [ Rx_complex,DC ] = dump_read16( file_path )
%SORA dump file reader
%2015-05-21 Li Songpeng

fid = fopen(file_path,'r');
% Data -> SORA int16 data, Count -> Number of data
[Data,Count]=fread(fid, inf,'int16'); 
fclose(fid);
% 
% Nsub = 64;
% %The first 8 data samples of every block are not real signal, drop them
% Ndrop = 8;
% signal = zeros(1,Count/Nsub*(Nsub-Ndrop)-1);
% for m=0:Count/64-1
%     signal((1:Nsub-Ndrop)+m*(Nsub-Ndrop)) = Data((1+Ndrop:Nsub)+m*Nsub); 
% end
% clear Data Count;
signal = Data;
Rx_I = signal(1:2:end-1);
Rx_Q = signal(2:2:end);
Rx_complex = Rx_I+1i*Rx_Q;

DC = mean(Rx_complex);
% DC = 0;
% Rx_AC = Rx_complex - DC;

end

