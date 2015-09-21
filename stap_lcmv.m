%%%LCMV%%%%%%%%
clc;
clear all;
close all;
j           = sqrt(-1);
M           = 4;                     %   node
lamda       = 0.01 ;                 %   lan
R           = 0.5*lamda ;            %   rada
N_snapshot  = 2048;                  %   snapshot  num

 N_bit      = 200;                   %   bit_num
 C          = 3e8;                   %   speed 
 fc         = 46.52e6;               %   fc
 fs         = 61.38e6;               %   fs
 Ts         = 1/fs ;                 %   TS
 
symbol_rate = 1.023e6;               %symbol_rate
 ca_rate    = 10.23e6;               %ca
 over_sample= fs/ca_rate;            %over_sample
 
 snr_db     = -20;                   %   snr  here fix niose
 sigma      = 1;
 snr        = 10^(snr_db/10);
 eb         = 2*snr^2*sigma;
 e_chip     = eb/1023;
 e_over     = e_chip/over_sample;
 
 
 isr_db0    = 60;
 isr0       = 10^(isr_db0/10);        % isr_one
 e_i0       = isr0*eb;
 a_coe0     = sqrt(e_i0);
 a_over0    = sqrt(e_i0/1023/over_sample);
 
 
 
 isr_db1    = 60;
 isr1       = 10^(isr_db1/10);        % isr_two
 e_i1       = isr1*eb;
 a_coe1     = sqrt(e_i1);
 a_over1    = sqrt(e_i1/1023/over_sample);
 
 
 
 isr_db2   = 60 ;
 isr2      = 10^(isr_db2/10);        % isr_three
 e_i2      = isr2*eb ;
 a_coe2    = sqrt(e_i2);
 a_over2   =  sqrt(e_i2/1023/over_sample);
 
 
 fj_0      =  45e6 ;
 fj_1      =  47e6 ;                 % fc_three
 fj_2      =  50e6 ; 
 
 bw_0      = 2.5e6;
 bw_1      = 2.5e6;                  % bw_three
 bw_2      = 2.5e6;  
 
 
 
 
 
 aj_0      = sqrt(e_i0/1023/over_sample/2);     % coeffcient three
 aj_1      = sqrt(e_i1/1023/over_sample/2);
 aj_2      = sqrt(e_i2/1023/over_sample/2);
 
 
 theta_wei = 0 ;
% theta_wei = pi*theta_wei/180;                  % satelite  angle
 phy_wei   = 0;
 phy_wei   = pi*(phy_wei/180);
 
 
 theta_0   = 20;
% theta_0   = pi*theta_0/180;                    % interference angle 0
 phy_0     = 100 ;
 phy_0     = pi*(phy_0/180);
 
 theta_1   = 40;
%  theta_1   = pi*theta_1/180;                    % interference angle 1
 phy_1     = 60 ;
 phy_1     = pi*(phy_1/180);
 
 theta_2   = 70;
% theta_2   = pi*theta_2/180;                    % interference angle 2
 phy_2     = 200 ;
 phy_2     = pi*(phy_2/180);
 
 
 
 x_index  = R*[1,-1,-1,1];                      % anita index
 y_index  = R*[1,1,-1,-1];
 
 angle_cell = [0:0.1:90];                       % angle cell
 
                            

%%%%%%%%%%%%%%%%%%%base_band%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bit_gen   =  randn(1,N_bit);
        bit_gen   =  bit_gen >0 ;
        bit_gen   =  2*bit_gen-1;                                          % source_bit_nz
        
        coeff_0   = [1,0,1,0,0,0,0,0,0,1];                                 % m_0 coefficient G1_poly=1+x^3+x^10 gps_c/a
        coeff_1   = [1,1,1,0,0,1,0,1,1,1];                                 % m_1 coefficient G2_poly=1+x^2+x^3+x^6+x^8+x^9+x^10 gps_c/a
        
        pn_code   = prn_code(coeff_0,coeff_1);
        pn_fix    = pn_code(3,:);
        pn_fix    = 2*pn_fix-1;                                            % pn_code_nz
        
        
        %spread_base = ds_mod(bit_gen ,pn_fix);
        
        %h_mod       = modem.pskmod('M',2);
        %mod_signal  = modulate(h_mod,spread_base);
        % ss         = real(mod_signal);
        
        
        ss = [];
        
         for i  =1:length(bit_gen)
             temp(1:1023)            = bit_gen(i);
             ss((i-1)*1023+1:i*1023) = temp.*pn_fix ;                       %  signal spread
            
            
        end
        
         for i= 1:length(ss)
            
            ss_os((i-1)*over_sample+1:i*over_sample)  = sqrt(e_over)*ss(i);  % oversample
            
         end
         
         
        n             =   [0:length(ss_os)-1];
      
      fc_ss           =   ss_os.*exp(j*2*pi*fc*n/fs);                   % fc
      
      ss_awgn         =   fc_ss+sigma*randn(1,length(fc_ss));                % add guass white noise
      
      interference_0  =  aj_0*cos(2*pi*fj_0*n/fs);                           % single    interference
      
   %  interference_1  =  aj_1*cos(2*pi*fj_1*n/fs+2*pi*delta_f1*(n/fs).^2);   % lfm       interference
    
      interference_1  =  aj_1*cos(2*pi*fj_1*n/fs);                          
      
      
      interference_2  =  aj_2*cos(2*pi*fj_2*n/fs);            
      
      
      ss_total      = ss_awgn + interference_0+interference_1+interference_2 ; % add   interference
      
      
      figure(1);
      
      fft_awgn     =  fft(ss_awgn);
      
      awgn_spec    = abs(fft_awgn).^2/length(fft_awgn);                       % ds   spectrum
      
      length_awgn  = [0:length(fft_awgn)-1]*fs/length(fft_awgn);
      
      awgn_spec    = awgn_spec(1:length(length_awgn));
      
      awgn_spec    = 10*log10(awgn_spec);
      
      
      plot(length_awgn,awgn_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      
      
      
      figure(2);
      
      ss_fft    =  fft(ss_total);
      ss_spec   = abs(ss_fft).^2/length(ss_fft);
      length_ss = [0:length(ss_fft)-1]*fs/length(ss_fft) ;
      
      ss_spec    = ss_spec(1:length(length_ss));
      ss_spec    = 10*log10(ss_spec);
      
      plot(length_ss,ss_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      
 
      
      for k=1:N_snapshot
          
    a1=fc_ss(k);
    
    a2=interference_0(k);
    
    a3=interference_1(k);
    
    
    ajam=interference_2(k);
    
   
    s(:,k)=a1*exp(j*2*pi*(R/lamda)*sin(theta_wei*pi/180)*[0:M-1]')+...
            +a2*exp(j*2*pi*(R/lamda)*sin(theta_0*pi/180)*[0:M-1]')+...
            +a3*exp(j*2*pi*(R/lamda)*sin(theta_1*pi/180)*[0:M-1]');
    jam(:,k)=ajam*exp(j*2*pi*(R/lamda)*sin(theta_2*pi/180)*[0:M-1]');
    
      end
      
      
x=s+jam;
Rx=x*x'/N_snapshot;                           
R=inv(Rx);                           
C=[1,0,0,0]';

Wopt=R*C/(C'*R*C);
%%%%%%%%%%%%%%%%%%%%%%%%%


for m=1:length(angle_cell);
    a=exp(j*2*pi*0.5*sin(angle_cell(m)*pi/180)*[0:M-1]');
    y(m)=Wopt'*a;
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
Y=20*log10(abs(y)/max(abs(y)));


plot(angle_cell,Y);hold on;grid on;
axis([-90 90 -50 0]);
plot(theta_wei,-30:0,'.');
plot(theta_0,-30:0,'.');
plot(theta_1,-30:0,'.');
plot(theta_2,-30:0,'.');
xlabel('\theta/o');
ylabel('Amplitude in dB');
title('LCMV-PI');
      
      
      
      
      
      angle_vector_wei  = [];
      
      
      for  i  =1:M
          
          x(i)         =  x_index(i);
          y(i)         =  y_index(i);
          z(i)         =  0 ;
          delay_wei(i) =  2*pi*(x(i)*sin(theta_wei)*cos(phy_wei)+y(i)*sin(theta_wei)*sin(phy_wei)+z(i)*cos(theta_wei))/lamda;   %standard
          
          angle_wei(i) =  exp(-j*delay_wei(i));
          angle_vector_wei  = [angle_vector_wei,angle_wei(i)];   
          
      end
 
 
      angle_vector_wei  = angle_vector_wei';
      
      
      
      angle_vector_0  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std0(i) = 2*pi*(x(i)*sin(theta_0)*cos(phy_0)+y(i)*sin(theta_0)*sin(phy_0)+z(i)*cos(theta_0))/lamda;  %standard
            angle_v_0(i)  = exp(-j*delay_std0(i));
            angle_vector_0= [angle_vector_0 ,angle_v_0(i)];
                
      end
 
     angle_vector_0    = angle_vector_0';
     
     
     
    angle_vector_1  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std1(i) = 2*pi*(x(i)*sin(theta_1)*cos(phy_1)+y(i)*sin(theta_1)*sin(phy_1)+z(i)*cos(theta_1))/lamda;  %standard
            angle_v_1(i)  = exp(-j*delay_std1(i));
            angle_vector_1= [angle_vector_1 ,angle_v_1(i)];
                
      end
 
     angle_vector_1    = angle_vector_1';
 
 
 
   angle_vector_2  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std2(i) = 2*pi*(x(i)*sin(theta_2)*cos(phy_2)+y(i)*sin(theta_2)*sin(phy_2)+z(i)*cos(theta_2))/lamda;  %standard
            angle_v_2(i)  = exp(-j*delay_std2(i));
            angle_vector_2= [angle_vector_2 ,angle_v_1(i)];
                
      end
 
     angle_vector_2    = angle_vector_2';
     









for k=1:N_snapshot;
    a1=1*randn(1);
    
    a2=10*randn(1);
    
    a3=10*randn(1);
    
    
    ajam=10*randn(1);
    
    an=1;
    s(:,k)=a1*exp(j*2*pi*d_lamda*sin(theta1*pi/180)*[0:N-1]')+...
            +a2*exp(j*2*pi*d_lamda*sin(theta2*pi/180)*[0:N-1]')+...
            +a3*exp(j*2*pi*d_lamda*sin(theta3*pi/180)*[0:N-1]');
    jam(:,k)=ajam*exp(j*2*pi*d_lamda*sin(theta_jam*pi/180)*[0:N-1]');
    n(:,k)=an*(randn(N,1)+j*randn(N,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%


x=s+jam+n;
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
plot(theta_wei,-30:0,'.');
plot(theta2,-30:0,'.');
plot(theta3,-30:0,'.');
plot(theta_jam,-30:0,'.');
xlabel('\theta/o');
ylabel('Amplitude in dB');
title('LCMV-PI');