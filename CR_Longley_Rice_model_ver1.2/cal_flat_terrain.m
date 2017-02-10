function [ flat_terrain, terrain_tx_to_rx, lenth_of_flat ] = cal_flat_terrain( tx_x, tx_y, rx_x, rx_y, terrain, show )
%CAL_FLAT_TERRAIN Summary of this function goes here
%[ flat_terrain, terrain_tx_to_rx ] = cal_flat_terrain( tx_x, tx_y, rx_x, r
% x_y, terrain )
%   Detailed explanation goes here
%   tx_x, tx_y = Location x & y for Tx
%   rx_x, rx_y = Location x & y for Rx
%   terrain = the terrain data
%   out put flat_terrain = the terrain data on the flat from Location of Tx
%    to Location of Rx
%   out put terrain_tx_to_rx = the terrain data from Location of Tx to Loca
%    tion of Rx
%   show = 0 when the show is not inputed
if nargin < 6 
    show = 0;
end
%tic
%这一段的目的是为了能读取terrain的值 因为小值永远要在前面
Loc_x_s = min(tx_x+1,rx_x+1);
Loc_x_l = max(tx_x+1,rx_x+1);
Loc_y_s = min(tx_y+1,rx_y+1);
Loc_y_l = max(tx_y+1,rx_y+1);

terrain_temp = terrain(Loc_y_s:Loc_y_l,Loc_x_s:Loc_x_l); %坐标从0开始 数据从1开始 数据先x后y 读取则先y后x
                                                         %这是一个重要的改动！！！
[N,M] = size(terrain_temp);
%这一段的目的是让起点永远在小值 终点永远在大值
if tx_x<=rx_x && tx_y<=rx_y
    terrain_tx_to_rx = terrain_temp;
elseif tx_x>rx_x && tx_y<=rx_y   %x方向颠倒
    for i = 1:M
        terrain_tx_to_rx(1:N,i) = terrain_temp(1:N,M-i+1);
    end
elseif tx_x<=rx_x && tx_y>rx_y
    for i = 1:N
        terrain_tx_to_rx(i,1:M) = terrain_temp(N-i+1,1:M);
    end
elseif tx_x>rx_x && tx_y>rx_y
    for i = 1:N
        for j = 1:M
            terrain_tx_to_rx(i,j) = terrain_temp(N-i+1,M-j+1);
        end
    end
end

[N,M] = size(terrain_tx_to_rx);
if N == 1
    flat_terrain = terrain_tx_to_rx;
elseif M == 1
    flat_terrain = terrain_tx_to_rx';
else
    if max(N,M) == N
        for i = 1:N
            if int32(i*M/N) == 0
                flat_terrain(i) = terrain_tx_to_rx(i,1);
            else
                flat_terrain(i) = terrain_tx_to_rx(i,int32(i*M/N));
            end
        end
    else
        for i = 1:M
            if int32(i*N/M) == 0
                flat_terrain(i) = terrain_tx_to_rx(1,i);
            else
                flat_terrain(i) = terrain_tx_to_rx(int32(i*N/M),i);
            end
        end
    end
end

lenth_of_flat = sqrt(N^2 + M^2);

if show ~= 0
    %显示成图形
    figure;
    [x,y] = meshgrid(1:M,1:N);
    subplot(2,1,1);surf(x,y,terrain_tx_to_rx);
    subplot(2,1,2);plot(flat_terrain);
    axis([0 max(M,N) -max(flat_terrain)*0.2 max(flat_terrain)*1.2]);
end
%toc
end

