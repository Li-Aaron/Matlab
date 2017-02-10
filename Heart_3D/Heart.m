function Heart
close all;clear all;clc;
figure
f=@(x,y,z)(x.^2+(9/4)*y.^2+z.^2-1).^3-x.^2.*z.^3-(9/80)*y.^2.*z.^3;
love_value=implicitsurf(f,[-1.5 1.5],[-.8 .8],[-1.5 1.5],52.1);
set(love_value,'AmbientStrength',.5);
rotate3d on;
%axis off;
shading interp;	
title('ÀÏ¹«°®Ð¡Æß');