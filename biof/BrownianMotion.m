function [ dx,dy,dz ] = BrownianMotion( stepSize,T,eta,d,fluidVel )
%BROWNIANMOTION Summary of this function goes here
%   Detailed explanation goes here
kb=1.38e-23;
D=kb*T/3/pi/eta/d;
k=sqrt(2*D*stepSize);
dx=k*randn()+fluidVel(1)*stepSize;
dy=k*randn()+fluidVel(2)*stepSize;
dz=k*randn()+fluidVel(3)*stepSize;


end

