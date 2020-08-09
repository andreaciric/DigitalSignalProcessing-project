function [b, N] = FIR_optim(Wp, Wa, Ap, Aa)

% fs = 48000;
dp = (10^(0.05*Ap)-1 )/(10^(0.05*Ap)+1);
da = 10^(-0.05*Aa);

Wa1 = Wa(1); fa1 = Wa1/pi;
Wa2 = Wa(2); fa2 = Wa2/pi;
Wa3 = Wa(3); fa3 = Wa3/pi;
Wa4 = Wa(4); fa4 = Wa4/pi;
Wa5 = Wa(5); fa5 = Wa5/pi;

Wp1 = Wp(1); fp1 = Wp1/pi;
Wp2 = Wp(2); fp2 = Wp2/pi;
Wp3 = Wp(3); fp3 = Wp3/pi;
Wp4 = Wp(4); fp4 = Wp4/pi;
Wp5 = Wp(5); fp5 = Wp5/pi;

B = [abs(Wa1-Wp1) abs(Wp2-Wa2) abs(Wa3-Wp3) abs(Wp4-Wa4) abs(Wa5-Wp5)];
Bt = min(B);

D = (0.01201*log10(dp)*log10(dp)+0.09664*log10(dp)-0.51325)*log10(da)+(0.00203*log10(dp)*log10(dp)-0.57054*log10(dp)-0.44314);
f = -16.9-14.6*(log10(dp)-log10(da));

M = ceil(2*pi*D/Bt - f*Bt/(2*pi) + 1)
N = M-1;

Hd = [1   1   0   0   1   1   0   0   1   1  0  0];
f1 = [0 fp1 fa1 fa2 fp2 fp3 fa3 fa4 fp4 fp5 fa5 1]; 
W = [1 dp/da 1 dp/da 1 dp/da];

while 1
    b = firpm(N, f1, Hd, W);
    
    Nfreqz = 15000; 
    [h, w] = freqz(b, 1, Nfreqz); 
    H = abs(h);

    ip1 = floor((Nfreqz*Wp1)/pi)+1;
    ip2 = ceil((Nfreqz*Wp2)/pi)+1;
    ip3 = floor((Nfreqz*Wp3)/pi)+1;
    ip4 = ceil((Nfreqz*Wp4)/pi)+1;
    ip5 = floor((Nfreqz*Wp5)/pi)+1;
    ia1 = ceil((Nfreqz*Wa1)/pi)+1;
    ia2 = floor((Nfreqz*Wa2)/pi)+1;
    ia3 = ceil((Nfreqz*Wa3)/pi)+1;
    ia4 = floor((Nfreqz*Wa4)/pi)+1;
    ia5 = ceil((Nfreqz*Wa5)/pi)+1;

    Hno = [H(ia1:ia2)' H(ia3:ia4)' H(ia5:Nfreqz)']; % nepropusni opsezi
    Hpo = [H(1:ip1)' H(ip2:ip3)' H(ip4:ip5)']; % propusni opsezi
    
    if((max(Hno) <= da) && (min(Hpo) >= (1-dp)) && ((max(Hpo) <= (1+dp))))
        break
    else
        N = N+1;
        M = M+1;
    end

end
end