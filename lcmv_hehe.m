%%%%%%%%????????????????????????-??????????????????????????????????????????%%%%%%%%%%%%%%%
  %%%%%%%author: wong %%%%%%%%%%%%%%%%%%%%
    %%%%LCMV-STAP  pipeline%%%%%%%%%%%%
    %%%%enhanced-adaptive-IIR%%%%%%%%
    %%%%%email:takeshineshiro@126.com%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 clc;
 clear  all ;
 
  M         = 4 ;                    %node
 lamda      = 0.01 ;                 %lan
 R          = 0.5*lamda;             %rada
 N_snapshot = 1024;                  %snapshot  num
 P          = 5 ;                    %delay_cell num
 N_bit      = 200;                   %bit_num
 C          = 3e8;                   %speed 
 fc         = 46.52e6;               %fc
 fs         = 61.38e6;               %fs
 Ts         = 1/fs ;                 %TS
 symbol_rate= 1.023e6;               %symbol_rate
 ca_rate    = 10.23e6;               %ca
 over_sample= fs/ca_rate;            %over_sample
 j          = sqrt(-1);              % j
 
 
 snr_db     = -20;                   % snr  here fix niose
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
 
 
 T0        = N_bit/1023;
 delta_f0  = bw_0/(2*T0);
 
 T1        = N_bit/symbol_rate;
 delta_f1  = bw_1/(2*T1);
 
 
 aj_0      = sqrt(e_i0/1023/over_sample/2);     % coeffcient three
 aj_1      = sqrt(e_i1/1023/over_sample/2);
 aj_2      = sqrt(e_i2/1023/over_sample/2);
 
 
 theta_wei = 10 ;
 theta_wei = pi*(theta_wei/180);                % satelite  angle
 phy_wei   = 0;
 phy_wei   = pi*(phy_wei/180);
 
 
 theta_0   = 20;
 theta_0   = pi*(theta_0/180);                  % interference angle 0
 phy_0     = 100 ;
 phy_0     = pi*(phy_0/180);
 
 theta_1   = 40;
 theta_1   = pi*(theta_1/180);                  % interference angle 1
 phy_1     = 60 ;
 phy_1     = pi*(phy_1/180);
 
 theta_2   = 70;
 theta_2   = pi*(theta_2/180);                  % interference angle 2
 phy_2     = 200 ;
 phy_2     = pi*(phy_2/180);
 
 
 
 x_index  = R*[1,-1,-1,1];                      % anita index
 y_index  = R*[1,1,-1,-1];
 
 
 
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
      
      fc_ss           =   pulse_base.*exp(j*2*pi*fc*n/fs);                   % fc
      
      ss_awgn         =   fc_ss+sigma*randn(1,length(fc_ss));                % add guass white noise
      
      interference_0  =  a_coe0*cos(2*pi*fj_0*n/fs);                           % single    interference
      
   %  interference_1  =  aj_1*cos(2*pi*fj_1*n/fs+2*pi*delta_f1*(n/fs).^2);   % lfm       interference
    
      interference_1  =  a_coe1*cos(2*pi*fj_1*n/fs);                          
      
      
      interference_2  =  a_coe2*cos(2*pi*fj_2*n/fs);            
      
      
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
      legend('??????????????????????????????');
      title('??????????????????????????????');
      
      
      figure(2);
      
      ss_fft    =  fft(ss_total);
      ss_spec   = abs(ss_fft).^2/length(ss_fft);
      length_ss = [0:length(ss_fft)-1]*fs/length(ss_fft) ;
      
      ss_spec    = ss_spec(1:length(length_ss));
      ss_spec    = 10*log10(ss_spec);
      
      plot(length_ss,ss_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('???????????????????????????');
      title ('???????????????????????????');
 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for   k  =  1:N_snapshot
               ss_mb     =  fc_ss(k);
               jam_mb_0  =  interference_0(k);
               jam_mb_1  =  interference_1(k);
               jam_mb_2  =  interference_2(k);
               
               
               
               source(:,k) = ss_mb*exp(j*2*pi*0.5*sin(theta_wei)*[0:M-1]');
               
               jam_0(:,k)  = jam_mb_0*exp(j*2*pi*0.5*sin(theta_0)*[0:M-1]');
               
               jam_1(:,k)  = jam_mb_1*exp(j*2*pi*0.5*sin(theta_1)*[0:M-1]');
               
               jam_2(:,k)  = jam_mb_2*exp(j*2*pi*0.5*sin(theta_2)*[0:M-1]');
               
               noise(:,k)  = sigma*randn(M,1);
       end
       
       
       xx_rev   =  source+jam_0+jam_1+noise;
       
       rr_xx    =  xx_rev*rr_rev'/N_snapshot;
       
       inv_rr   = inv(rr_xx);
       
       
       
       f_theta0 =  [1,0,0,0]';
       
       w_opt    =  inv_rr*f_theta0/(f_theta0'*inv_rr*f_theta0);
       
       
       theta_range  =  [0:0.01:90];
       
       
       for i  =  1:length(theta_range)
           
           angle_cell   =  exp(j*2*pi*0.5*sin(theta_range(i)*pi/180)*[0:M-1]');
           y_ww(i)      =  w_opt'*angle_cell;
           
       end
       
       y_av   =  20*log10(abs(y_ww)/max(abs(y_ww)));
       
       
       figure(3);
       
       plot(theta_range,y_av);
       
       title('az');
       xlabel('angle(degree)');
       ylabel('db');
       
       
       
       for   i  =  1:N_sapshot
           xx_i   =  xx_rev(:,i);
           
           y(i)   =  w_opt'*xx_i;
           
           
           
       end
       
       
         
  
         
      figure(3);
      
      ss_fft    =  fft(y);
      ss_spec   = abs(ss_fft).^2/length(ss_fft);
      length_ss = [0:length(ss_fft)-1]*fs/length(ss_fft) ;
      
      ss_spec    = ss_spec(1:length(length_ss));
      ss_spec    = 10*log10(ss_spec);
      
      plot(length_ss,ss_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('SMI-STAP ??????????????????????????????');
      title ('SMI-STAP ??????????????????????????????');
       
      
      
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 