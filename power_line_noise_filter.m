function [bd, ad] = power_line_noise_filter(fs, fc, Aa, Ap)

    Aav = Aa;
    Apv = Ap;
    Nfreqz = 15000;

    % odgovarajuce ucestanosti
    fp1 = fc - 2; 
    fp2 = fc + 2;
    fa1 = fc - 0.5; 
    fa2 = fc + 0.5;

    % diskretne granicne ucestanosti propusnog opsega
    Wpd1 = 2*pi*fp1/fs; 
    Wpd2 = 2*pi*fp2/fs;
    % diskretne granicne ucestanosti nepropusnog opsega
    Wad1 = 2*pi*fa1/fs; 
    Wad2 = 2*pi*fa2/fs;
    
    % predistorzovane granicne ucestanosti propusnog opsega
    wpa1 = 2*fs*tan(Wpd1/2); 
    wpa2 = 2*fs*tan(Wpd2/2);
    % predistorzovane granicne ucestanosti nepropusnog opsega
    waa1 = 2*fs*tan(Wad1/2); 
    waa2 = 2*fs*tan(Wad2/2);

    Ka = tan(Wpd2/2) - tan(Wpd1/2);
    Kb = tan(Wpd1/2) * tan(Wpd2/2);
    Kc = tan(Wad1/2) * tan(Wad2/2);
    
    B = wpa2 - wpa1;

    % zbog nelinearne transformacije, nakon predistorzije ucestanosti, nece
    % biti ispunjen uslov da je
    % w0 = sqrt(wp1_pred*wp2_pred) = sqrt(wa1_pred*wa2_pred),
    % zbog toga se uzima da je w0 = sqrt(wp1_pred*wp2_pred) i zadrzava se jos
    % jedna predistorzovana ucestanost, ili wa1_pred, ili wa2_pred, dok se
    % preostala ucestanost izracunava na osnovu ostale tri, kao sto je uradjeno
    w0 = sqrt(wpa1*wpa2);
    
    if (w0^2 < waa1*waa2)
        waa2 = w0^2/waa1;
    else if (w0^2 > waa1*waa2)
        waa1 = w0^2/waa2;
        end
    end
    
    wpp = 1; % granicna ucestanost propusnog opsega analognog prototipa
    wap = B*waa1/(w0^2-waa1^2); % granicna ucestanost nepropusnog opsega analognog prototipa

    k1 = (Ka*tan(Wad1/2))/(Kb - (tan(Wad1/2))^2);
    k2 = (Ka*tan(Wad2/2))/((tan(Wad2/2))^2 - Kb);
   
    if (Kb < Kc)
        k = 1/k2;
    else
        k = 1/k1;
    end
       
    %k = wpp/wap;
    k1 = sqrt(1-k^2);
    q0 = 1/2*(1-sqrt(k1))/(1+sqrt(k1));
    q = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13; 
  
while (1)

    D = (10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
   
    N = ceil((log10(16*D))/(log10(1/q))); % potreban red filtra, prema formuli
        
    [z,p,k] = ellipap(N,Apv,Aav);

    baN = k*poly(z);
    aaN = poly(p);
    
    %Transformacija normalizovanog prototipa u VF 
    [ba,aa] = lp2bs(baN,aaN,w0, B);
    
    % Diskretizacija bilinearnom transformacijom: nule i polovi analognog -> nule i polove digitalnog
    [bd, ad] = bilinear(ba,aa,fs); 
    
    [hd,wd] = freqz(bd,ad,Nfreqz);
    Hd = abs(hd);
        
   % provera gabarita
    
    df = (fs/2)/Nfreqz;

    NPOK = 0;
    POOK = 0;

    ia1 = floor(fa1/df)+1;
    ia2 = ceil(fa1/df)+1;
    ip1 = ceil(fp1/df)+1;
    ip2 = floor(fp2/df)+1;
    
    Hp=[Hd(1:ip1)' Hd(ip2:length(Hd))'];
    Ha=Hd(ia1:ia2);

    if(max(20*log10(Ha))>-Aa )
        Aav=Aav+0.1;
    else
        NPOK=1;
    end

    if(min(20*log10(Hp))<-Ap)
        Apv=Apv-0.1;
    else
        POOK=1;
    end   

    if((NPOK==1)&&(POOK==1))
        break
    end

end 
N
    
end