%%%%%%%%GPS接收机抗干扰算法及实现研究%%%%%%%%%%%%%%%
  %%%%%%%author: wong %%%%%%%%%%%%%%%%%%%%
       %%%%stap algorithm %%%%%%%%%%%%
    %%%%%email:takeshineshiro@126.com%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clc;
 clear  all ;
 
  M         = 4 ;                    %node
 lamda      = 0.01 ;                 %lan
 R          = 0.5*lamda;             %rada
 N_snapshot = 4096;                  %snapshot  num
 P          = 5 ;                    %delay_cell  num
 N_bit      = 200;                   %bit_num
 C          = 3e8;                   %speed 
 fc         = 2.048e6;               %fc
 fs         = 4.096e6;               %fs
 Ts         = 1/fs ;                 %TS
 ca_rate    = 1.024e6;               %ca
 over_sample= fs/ca_rate;            %over_sample
 j          = sqrt(-1);              % j
 
 
 snr_db     = -20;                   % snr  here fix niose
 sigma      = 1;
 snr        = 10^(snr_db/10);
 eb         = 2*snr^2*sigma;
 e_chip     = eb/1023;
 e_over     = e_chip/over_sample;
 
 
 isr_db    = 60;
 isr       = 10^(isr_db/10);              % single  interference
 e_i       = isr*eb;
 a_coe     = sqrt(e_i);
 a_over    = sqrt(e_i/1023/over_sample);
 
 
 

 
 
 fj_0      =  0.5e6 ;
 fj_1      =  1e6 ;                            %two  interference
 
 aj_0      =  sqrt(e_i/1023/over_sample/2);
 
 aj_1      =  sqrt(e_i/1023/over_sample/2);
 
 
 
 
 
 aj_0      = sqrt(e_i0/1023/over_sample/2);     % coeffcient three
 aj_1      = sqrt(e_i1/1023/over_sample/2);
 aj_2      = sqrt(e_i2/1023/over_sample/2);
 
 
 theta_wei = 0 ;
 theta_wei = pi*(theta_wei/180);                % satelite  angle
 phy_wei   = 0;
 phy_wei   = pi*(phy_wei/180);
 
 
 theta_0   = 60;
 theta_0   = pi*(theta_0/180);                  % interference angle 0
 phy_0     = 100 ;
 phy_0     = pi*(phy_0/180);
 
 theta_1   = 70;
 theta_1   = pi*(theta_1/180);                  % interference angle 1
 phy_1     = 60 ;
 phy_1     = pi*(phy_1/180);
 
 theta_2   = 30;
 theta_2   = pi*(theta_2/180);                  % interference angle 2
 phy_2     = 200 ;
 phy_2     = pi*(phy_2/180);
 
 
 
 x_index  = R*[1,-1,-1,1];                      % anita index
 y_index  = R*[1,1,-1,-1];
 
 
 mu       =  0.12 ;                             %  step_factor >= 5
 
 % S_index(m)  =  [R*cos(2*pi*m/M),R*sin(2*pi*m/M),0]'         % space_index
 
 %  delay(m)   = sin(theta)*cos(fine-2*pi*m/M)*2*pi*R)/lamda   % node delay
 
 %  angle_vector(m)  =  exp(-j*delay(m));                      % angle

 
  %%%%%%%%%%%%%%%%%%%base_band%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bit_gen   =  randn(1,N_bit);
        bit_gen   =  bit_gen >0 ;
        bit_gen   =  2*bit_gen-1;                                          % source_bit_nz
        
        coeff_0   = [1,0,1,0,0,0,0,0,0,1];                                 % m_0 coefficient G1_poly=1+x^3+x^10 gps_c/a
        coeff_1   = [1,1,1,0,0,1,0,1,1,1];                                 % m_1 coefficient G2_poly=1+x^2+x^3+x^6+x^8+x^9+x^10 gps_c/a
        
        pn_code   = prn_code(coeff_0,coeff_1);
        pn_fix    = pn_code(1,:);
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
            
            ss_os((i-1)*over_sample+1:i*over_sample)  = sqrt(a_over)*ss(i);  % oversample
            
         end
         
         
        n             =   [0:length(ss_os)-1];
      
      fc_ss           =   pulse_base.*exp(j*2*pi*fc*n/fs);                   % fc
      
      ss_awgn         =   fc_ss+sigma*randn(1,length(fc_ss));                % add guass white noise
      
      interference_0  =  aj_0*cos(2*pi*fj_0*n/fs);                           % single    interference
      
   %  interference_1  =  aj_1*cos(2*pi*fj_1*n/fs+2*pi*delta_f1*(n/fs).^2);   % lfm       interference
    
      interference_1  =  aj_1*cos(2*pi*fj_1*n/fs);                          
      
      
    %  interference_2  =  aj_2*cos(2*pi*fj_2*n/fs);            
      
      
      ss_total      = ss_awgn + interference_0+interference_1;               % add   interference
      
      
      figure(1);
      
      fft_awgn     =  fft(ss_awgn);
      
      awgn_spec    = abs(fft_awgn).^2/length(fft_awgn);                       % ds   spectrum
      
      length_awgn  = [0:length(fft_awgn)-1]*fs/length(fft_awgn);
      
      awgn_spec    = awgn_spec(1:length(length_awgn));
      
      awgn_spec    = 10*log10(awgn_spec);
      
      
      plot(length_awgn,awgn_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('未加干扰时信号功率普');
      title('未加干扰时信号功率普');
      
      
      figure(2);
      
      ss_fft    =  fft(ss_total);
      ss_spec   = abs(ss_fft).^2/length(ss_fft);
      length_ss = [0:length(ss_fft)-1]*fs/length(ss_fft) ;
      
      ss_spec    = ss_spec(1:length(length_ss));
      ss_spec    = 10*log10(ss_spec);
      
      plot(length_ss,ss_spec);
      
      xlabel('Hz');
      ylabel('power spectrum(dB)');
      legend('加干扰时信号功率普');
      title ('加干扰时信号功率普');
 
      
      
      
%       angle_vector_wei  = [];
%       
%       
%       for  i  =1:M
%           
%           x(i)         =  x_index(i);
%           y(i)         =  y_index(i);
%           z(i)         =  0 ;
%           delay_wei(i) =  2*pi*(x(i)*sin(theta_wei)*cos(phy_wei)+y(i)*sin(theta_wei)*sin(phy_wei)+z(i)*cos(theta_wei))/lamda;   %standard
%           
%           angle_wei(i) =  exp(-j*delay_wei(i));
%           angle_vector_wei  = [angle_vector_wei,angle_wei(i)];   
%           
%       end
%  
%  
%       angle_vector_wei  = angle_vector_wei';
      
      
      
      angle_vector_0  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std0(i) = 2*pi*(x(i)*sin(theta_0)*cos(phy_0)+y(i)*sin(theta_0)*sin(phy_0))+z(i)*cos(theta_0))/lamda;  %standard
            delay_0(i)    =  sin(theta_0)*cos(phy_0-2*pi*i/M)*R*2*pi/lamda;  %  angle_vector_0  
            angle_v_0(i)  = exp(-j*delay_0(i));
            angle_vector_0= [angle_vector_0 ,angle_v_0(i)];
                
      end
 
     angle_vector_0    = angle_vector_0';
     
     
     
    angle_vector_1  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std1(i) = 2*pi*(x(i)*sin(theta_1)*cos(phy_1)+y(i)*sin(theta_1)*sin(phy_1))+z(i)*cos(theta_1))/lamda;  %standard
            delay_1(i)    =  sin(theta_1)*cos(phy_1-2*pi*i/M)*R*2*pi/lamda;  %  angle_vector_1 
            angle_v_1(i)  = exp(-j*delay_1(i));
            angle_vector_1= [angle_vector_1 ,angle_v_1(i)];
                
      end
 
     angle_vector_1    = angle_vector_1';
 
 
 
   angle_vector_2  = [];
      
      for  i  =1:M
            x(i)          = x_index(i);
            y(i)          = y_index(i);
            z(i)          = 0;
            delay_std2(i) = 2*pi*(x(i)*sin(theta_2)*cos(phy_2)+y(i)*sin(theta_2)*sin(phy_2))+z(i)*cos(theta_2))/lamda;  %standard
            delay_2(i)    =  sin(theta_2)*cos(phy_2-2*pi*i/M)*R*2*pi/lamda;  %  angle_vector_0  
             
            angle_v_2(i)  = exp(-j*delay_2(i));
            angle_vector_2= [angle_vector_2 ,angle_v_2(i)];
                
      end
 
     angle_vector_2    = angle_vector_2';
     
     
     angle_matrix      =  [angle_vector_wei, angle_vector_0,angle_vector_1, angle_vector_2];
     
     
     
      ss_vector_0     = fc_ss;
      %ss_vector_0    = fc_ss +sigma*randn(1,length(fc_ss));
      ss_vector_1     = interference_0;
      ss_vector_2     = interference_1;
   %   ss_vector_3     = interference_2;
      
      ss_matrix       = [ss_vector_0,ss_vector_1,ss_vector_2 ];
      
      
      vv_matrix       = [];
      
      
      for i=  1:M
          
          vv_vector(i,:)   =  sigma*randn(1,length(fc_ss));
          
      end
      
   
      
 
 
      xx_vector = [];
      
      for  i   = 1:N_snapshot
          
           xx_vector(:,i)  =angle_matrix(:,i)+vv_vector(:,i);
           %xx_vector(:,i) =  angle_matrix*ss_matrix(:,i);
          
          
      end
      
     
   xx_vector  =xx_vector;
   
   
   
   xx_rem   = [];
   
   
   
   for   i  = 1:size(ss_matrix,2)-N_snapshot
       
         xx_rem(:,i)  = angle_matrix*ss_matrix(:,i+N_snapshot);            % x(n)= A*s(n)+v(n)    n  = 0-P
       
   end
         
 
  
   ww_0       = [1,zeros(1,M-1)];
   ww_zero    = zeros(1,M);                                                % all node at  n
   ww_init    = [ww_0 ,ww_zero,ww_zero,ww_zero, ww_zero];
   ww_init_R  =  ww_init';
   
   
   
   xx_p_block  = [];                                                       % MP *1 
   
   ww_p_block  = [];                                                       % MP *1
   
   
   
   
    w0   =  [1,0,0,0]';                                                     %  w  initial value
   
   s   =   [1,0,0,0]';                                                     %  s
   
   s0  =   [0,0,0,0]';                                                     %  0
   
   s_matrix  = [s';s0';s0';s0'];                                           %  s*s^T/s^T*s
   
   eye_M    =  eye(M);                                                     %  I
   
   
   sub_M    = eye_M - s_matrix;                                            % I-s*s^T/s^T*s
   
   
   w        =  [];                                                         %  w  iterator
   
   
   
    for  i  =  1:N_snapshot
       
        xx_i       = xx_vector(:,i);
        xx_i_conj  = conj(xx_i);
        xx_i_T     = xx_i';
        sub_mul    = sub_M*xx_i_conj ;
        
        xx_mul     = sub_mul*xx_i_T ;
        
        coeff_xx   = 2*mu*xx_mul;
        
        if(i==1)
             w0_temp = coeff_xx*w0;
             w(:,i)  = w0-w0_temp;
             
        else
            w_temp  = coeff_xx*w(:,i-1);
            w(:,i)  = w(:,i-1)-w_temp;
            
        end
               y(i)  =  w(:,i)'*xx_i;
            
        
       
    end
   
   
   
      for  i  =1:size(xx_rem,2)
       
       y_rem(i)  = w(:,N_snapshot)'*xx_rem(:,i);
       
       
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
      legend('PI抗干扰后输出信号频谱');
      title ('PI 抗干扰后输出信号频谱');
       