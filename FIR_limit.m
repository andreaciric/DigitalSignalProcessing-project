function [b, N] = FIR_limit(Wp, Wa, Ap, Aa)

dp = (10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
da = 10^(-0.05*Aa);

delta = min(dp, da);
if delta ~= da
    Aa = -20*log10(delta);
end

if (Aa < 21)
    beta = 0;
else if (Aa <= 50)
        beta = 0.5842*(Aa-21)^0.4+0.07886*(Aa-21);
    else
        beta = 0.1102*(Aa-8.7);
    end
end

if Aa <= 21
    D = (Aa-7.95)/14.36; 
else
    D = (Aa-7.95)/14.36; 
end

Wa1 = Wa(1);
Wa2 = Wa(2);
Wa3 = Wa(3);
Wa4 = Wa(4);
Wa5 = Wa(5);

Wp1 = Wp(1);
Wp2 = Wp(2);
Wp3 = Wp(3);
Wp4 = Wp(4);
Wp5 = Wp(5);

%NO:
Bt1 = min(Wa1-Wp1, Wp2-Wa2);
Wc1 = Wp1+Bt1/2;
Wc2 = Wp2-Bt1/2;

Bt2 = min(Wa3-Wp3, Wp4-Wa4);
Wc3 = Wp3+Bt2/2;
Wc4 = Wp4-Bt2/2;

%NF:
Bt3 = Wa5-Wp5;
Wc5 = (Wa5+Wp5)/2;

Bt = min([Bt1 Bt2 Bt3]);

M = ceil(2*pi*D/Bt+1)
N = M-1; 

while 1
    window = kaiser(M, beta)';

    n = -(M-1)/2:(M-1)/2;
    b = (-sin(n*Wc1)+sin(n*Wc2)-sin(n*Wc3)+sin(n*Wc4)-sin(n*Wc5))./(n*pi);
    if (mod(M, 2) == 1)
        i = (M+1)/2;
        b(i) = (-Wc1+Wc2-Wc3+Wc4-Wc5)/pi;
    end

    b = b.*window;
    
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