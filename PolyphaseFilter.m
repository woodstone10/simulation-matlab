%-----------------------------------------------------------
%
%
%  Polyphase Filter  --- Last Version Fig.7.48 (Rabiner Book)
%  Multirate Digital Signal Processing
% 
%-----------------------------------------------------------
%
% January, 2000
%
% Created by Jonggil Nam
% https://www.linkedin.com/in/jonggil-nam-6099a162/ | https://github.com/woodstone10 | woodstone10@gmail.com | +82-10-8709-6299 
%-----------------------------------------------------------
clear all
close all

N=1; % number of PNCODE
F=1 ; % 1 is original
CC=-1;  % -1 is original
PN_length=128;  % length of PN code
OFF=0.625  % 0 or 0.625 , 16 band is 0.625
CUT_FREQ=1.25 % 1.25 is original
Xn = SSBsig(N,F,CC,PN_length,OFF,CUT_FREQ) ;  % SSB signal--(1 3 5 7)
%Want_FA=1;  Xn=ONEsig(Want_FA,OFF,CUT_FREQ);  % one FA--select Want_FA
%Xn = DSBsig(N,F) ; % DSB signal


%-------noise----------------%
%SNRindB=0 , NOISE=AWGN(PN_length,N,Xn,SNRindB);

%----- input signal -----------%
%Xn1 = (Xn+NOISE).'  ;
Xn1=Xn.' ;

$----- polyphase filter ------------%
NB = 16  % Subband
NC = 128  % Coefficient of Filter
KL  = NB*2  % Cutoff frequency of filter = 1/KL 
Bit_index=4

%  Initialization                        
BB= NC/NB+1 ;	                       % buffer size (default 16)
orig_X = zeros ( NB,PN_length*N*32/NB ) ;     %  original input 
anal_X = zeros ( NB,PN_length*N*32/NB);      % analysis input
syn_f = zeros ( NB, BB-1)  ;    % lowpass synthesis filter
X_BUF = zeros( NB, BB ) ;       % input buffer
Y_BUF = zeros( NB, BB ) ;
fft_X = zeros( NB, 1 ) ;              % FFT input 
ifft_X = zeros( NB, 1 ) ;             % FFT 
ifft_Y = zeros( NB, 1 ) ;             % IFFT 
anal_poly = zeros(NB, BB);      % analysis polyphase filter
syn_poly = zeros(NB, BB);	      % synthesis polyphase filter
syn_Y = zeros( NB, 1 ) ; 	
syn_X = zeros(1,PN_length*N*32) ;  
rea_ifft_X = zeros( NB, 1 ) ;  
      
% Polyphase filter of analysis 
BIT=2^Bit_index ;


%FIL_COEF =  fir1(NC-1, CUT_FREQ/40)*(BIT/311)   ;
%anal_proto=fix((FIL_COEF*1000))/1000; %*(BIT/311);
anal_proto=fir1(NC-1, CUT_FREQ/40)   ;
%anal_proto =  fir1(NC-1,1.25/40*11/12 ) ;
%M=NC-1;
%[anal_proto] = kaiserfil(M) ;

jm=[];
for int=1:NB
      Wm=[]; 
    for int2=1:BB
      Wm=[Wm,exp(j*pi*(int2-1)/2) ];
   end
   jm=[jm;Wm];
end
 
anal_poly=[];
for int = 1:NB
  anal_f=[];
  if int ==1 
     anal_f1 = [anal_proto( 1 : NB : NC  ),0] ;
     anal_f=[anal_f,anal_f1];
  else 
     anal_f1 =[0, anal_proto( NB+2-int : NB : NC )] ;
     anal_f=[anal_f,anal_f1];
  end
     anal_poly=[anal_poly ; anal_f];
end  
del_anal_poly=jm.*anal_poly ;
del_anal_poly(1,1)=0 ;

% Polyphase filter of synthesis
syn_proto=anal_proto;
for int = 1:NB
	syn_f(int , :) = syn_proto( int : NB : NC );   
    syn_poly(int,1:BB-1) = syn_f(int , :);
end
del_syn_poly=jm.*syn_poly;
     
     
%======== Routine ==============================================
ex_data=[];
save_fft_X=[];
for	ind = 1 : PN_length*N*32/NB	
    orig_X = reshape ( Xn1, NB, PN_length*N*32/NB ) ;	   
    anal_X(: , ind) = orig_X ( : , ind ) ;   
                
    X_BUF(: , 2 : BB) = X_BUF(: , 1 : BB-1) ;
    % buffer shifting
    X_BUF(:, 1 ) = anal_X(: , ind ) ;	
    %------polyphase filter--------%         
    for  t = 1 : NB 		
       fft_X1(t, 1 ) = X_BUF (t,1:BB)*del_anal_poly(t , :).'  ; 
    end
     %------WM-(M-1)/4----------%
    for  tt = 1:NB
       fft_X(tt,1)   = fft_X1(tt, 1 )* exp(-j*pi*(tt-1)/(2*NB));
    end
    save_fft_X=[save_fft_X , fft_X];
    %---------DFT---------------%
    ifft_X= fft ( fft_X ( : , 1 )) ;	 
    ex_data=[ex_data , ifft_X];
    %--------Re{}---------------%
    rea_ifft_X = real(ifft_X);
    %--- or ---%
    rea_ifft_X2(1,1)=0;
    rea_ifft_X2(2,1)=0;
    rea_ifft_X2(3,1)=rea_ifft_X(3,1);
    rea_ifft_X2(4,1)=rea_ifft_X(4,1);
    rea_ifft_X2(5,1)=rea_ifft_X(5,1);
    rea_ifft_X2(6,1)=rea_ifft_X(6,1);
    rea_ifft_X2(7:16,1)=0 ;     
    %--------IDFT---------------%        
    ifft_Y = ifft(rea_ifft_X2(:,1)) ;
    %--------Wm----------------%       
    for tt = 1:NB
       ifft_YY(tt,1) = ifft_Y(tt,1)*exp(j*pi*(tt-1)/(2*NB));
    end
    %------ polyphase filter ------%
    Y_BUF(:,2:BB) = Y_BUF(:,1:BB-1) ;	
    Y_BUF(:,1) = ifft_YY ;	
    for  t = 1 : NB 		
        syn_Y(t, 1) = Y_BUF(t,1:BB)*del_syn_poly(t, :).'  ;
    end
    %------- reshape------------%   
    syn_X((ind-1)*NB+1:ind*NB) = reshape(real(syn_Y(:,1)),1,NB) ;
end
%=============Routine===============================


%-------- Plot -----------------%
figure(1)
as1=fft(Xn1(PN_length*N*32-1999:PN_length*N*32-1000),1000);
max_as1=max(as1);
as2=fft(syn_X(PN_length*N*32-1999:PN_length*N*32-1000),1000);
max_as2=max(as2);
subplot(411);
plot(abs(as1))
title('Spectrum of Input Signal')
pos=1000/32;
line('xdata',[1000/32 1000/32 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*1 1000/32+pos*1 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*2 1000/32+pos*2 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*3 1000/32+pos*3 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*4 1000/32+pos*4 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*5 1000/32+pos*5 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*6 1000/32+pos*6 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*7 1000/32+pos*7 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*8 1000/32+pos*8 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*9 1000/32+pos*9 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*10 1000/32+pos*10 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*11 1000/32+pos*11 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*12 1000/32+pos*12 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*13 1000/32+pos*13 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*14 1000/32+pos*14 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*15 1000/32+pos*15 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*16 1000/32+pos*16 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*17 1000/32+pos*17 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*18 1000/32+pos*18 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*19 1000/32+pos*19 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*20 1000/32+pos*20 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*21 1000/32+pos*21 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*22 1000/32+pos*22 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*23 1000/32+pos*23 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*24 1000/32+pos*24 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*25 1000/32+pos*25 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*26 1000/32+pos*26 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*27 1000/32+pos*27 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*28 1000/32+pos*28 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*29 1000/32+pos*29 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*30 1000/32+pos*30 ],'ydata',[0 200],'linestyle',':')
line('xdata',[1000/32+pos*31 1000/32+pos*31 ],'ydata',[0 200],'linestyle',':')
subplot(412);
plot(abs(max_as1/max_as2*as2),'b')
title('Spectrum of Output Signal')
pos2=1000/8;
line('xdata',[pos2 pos2 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*2  pos2*2 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*3  pos2*3 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*4  pos2*4 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*5  pos2*5 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*6  pos2*6 ],'ydata',[0 200],'linestyle','-')
line('xdata',[pos2*7  pos2*7 ],'ydata',[0 200],'linestyle','-')
STA=3001;
END=3101;
max_syn = max(syn_X(STA+NC:END+NC)) ;
max_x  = max(Xn1(STA:END)) ;
subplot(413);
plot(Xn1(STA:END),'r');
title('Input Signal in Time Domain')
axis([0 END-STA -2 2])
grid on   
subplot(414); 
plot(max_x/max_syn*syn_X(STA+NC:END+NC),'b');
title('Output Signal in Time Domain')
axis([0 END-STA -2 2])
grid on
%gtext('# of coefficient in filter = 128')
%gtext('3,4,5,6th branch input')
%gtext('Using Kaiser Filter')
%gtext('9 bit')


figure(2)
A=0.5;
subplot(8,2,1)
plot(abs(real(ex_data(1,:))))
ylabel('1st')
axis([0 600  0 A])
subplot(8,2,3)
plot(abs(real(ex_data(2,:))))
ylabel('2nd')
axis([0 600  0 A])
subplot(8,2,5)
plot(abs(real(ex_data(3,:))))
ylabel('3rd')
axis([0 600  0 A])
subplot(8,2,7)
plot(abs(real(ex_data(4,:))))
ylabel('4th')
axis([0 600  0 A])
subplot(8,2,9)
plot(abs(real(ex_data(5,:))))
ylabel('5th')
axis([0 600  0 A])
subplot(8,2,11)
plot(abs(real(ex_data(6,:))))
ylabel('6th')
axis([0 600  0 A])
subplot(8,2,13)
plot(abs(real(ex_data(7,:))))
ylabel('7th')
axis([0 600  0 A])
subplot(8,2,15)
plot(abs(real(ex_data(8,:))))
ylabel('8th')
axis([0 600  0 A])
subplot(8,2,2)
plot(abs(real(ex_data(1+8,:))))
ylabel('9th')
axis([0 600  0 A])
subplot(8,2,4)
plot(abs(real(ex_data(2+8,:))))
ylabel('10th')
axis([0 600  0 A])
subplot(8,2,6)
plot(abs(real(ex_data(3+8,:))))
ylabel('11th')
axis([0 600  0 A])
subplot(8,2,8)
plot(abs(real(ex_data(4+8,:))))
ylabel('12th')
axis([0 600  0 A])
subplot(8,2,10)
plot(abs(real(ex_data(5+8,:))))
ylabel('13th')
axis([0 600  0 A])
subplot(8,2,12)
plot(abs(real(ex_data(6+8,:))))
ylabel('14th')
axis([0 600  0 A])
subplot(8,2,14)
plot(abs(real(ex_data(7+8,:))))
ylabel('15th')
axis([0 600  0 A])
subplot(8,2,16)
plot(abs(real(ex_data(8+8,:))))
ylabel('16th')
axis([0 600  0 A])




