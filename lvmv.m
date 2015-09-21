%%%LCMV%%%%%%%%
clc;clear all;close all;
j=sqrt(-1);
N=4;                               
d_lamda=0.5;                  
theta=0:0.1:90;            
theta1=0;                      
theta2=20;                          
theta3=40;                        
theta_jam=70;                  
L=512;                              
%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:L;
    a1=0.0000001*randn(1);
    
    a2=0.127*randn(1);
    
    a3=0.127*randn(1);
    
    
    ajam=0.127*randn(1);
    
    an=1;
    s(:,k)=a1*exp(j*2*pi*d_lamda*sin(theta1*pi/180)*[0:N-1]')+...
            +a2*exp(j*2*pi*d_lamda*sin(theta2*pi/180)*[0:N-1]')+...
            +a3*exp(j*2*pi*d_lamda*sin(theta3*pi/180)*[0:N-1]');
    jam(:,k)=ajam*exp(j*2*pi*d_lamda*sin(theta_jam*pi/180)*[0:N-1]');
    n(:,k)=an*randn(N,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%


x=s+jam;
Rx=x*x'/L;                           
R=inv(Rx);                           
C=[1,0,0,0]';

Wopt=R*C/(C'*R*C);
%%%%%%%%%%%%%%%%%%%%%%%%%


for m=1:length(theta);
    a=exp(j*2*pi*d_lamda*sin(theta(m)*pi/180)*[0:N-1]');
    y(m)=Wopt'*a;
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
Y=20*log10(abs(y)/max(abs(y)));


plot(theta,Y);hold on;grid on;
axis([-90 90 -50 0]);
plot(theta1,-30:0,'.');
plot(theta2,-30:0,'.');
plot(theta3,-30:0,'.');
plot(theta_jam,-30:0,'.');
xlabel('\theta/o');
ylabel('Amplitude in dB');
title('LCMV-PI');