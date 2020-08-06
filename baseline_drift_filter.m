function [bd, ad] = baseline_drift_filter(fs, fa, fp, Aa, Ap)

    T = 1/fs;
    Nfreqz = 15000; % broj tacaka za izracunavanje frekvencijske karakteristike
        
    Aav = Aa;
    Apv = Ap;
    
    Wpd=2*pi*fp/fs; % diskretna granicna ucestanost propusnog opsega
    Wad=2*pi*fa/fs; % diskretna granicna ucestanost nepropusnog opsega
    
    % PREDISTORZIJA granicnih ucestanosti
    % posto ce se koristiti bilinearna transformacija za diskretizaciju
% analognog filtra, mora se izvrsiti predistorzija ucestanosti
    waa = 2/T*tan(Wad/2); % predistorzovana granicna ucestanost propusnog opsega
    wpa = 2/T*tan(Wpd/2); % predistorzovana granicna ucestanost nepropusog opsega
    
    wpp = 1; % normalizacija / % granicna ucestanost propusnog opsega analognog prototipa
    wap = wpa/waa; % granicna ucestanost nepropusnog opsega analognog prototipa

    % TEORIJA NF prototip -> VF     
    % wpN = a/wpa => a = 2*wpN*tan(Wpa/2)/T
    % waN = a/waa
    % k = wpN/waN = waa/wpa = tan(Waa/2)/tan(Wpa/2)

while(1)    
    
    D = (10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
    k = wpp/wap;
    
    N = ceil(acosh(sqrt(D))/acosh(1/k));    
   % N = uint8(N);
        
    [z,p,k]=cheb2ap(N,Aav);

    baN=k*poly(z);
    aaN=poly(p);
    
    %Transformacija normalizovanog prototipa u VF 
    [ba,aa]=lp2hp(baN,aaN,waa);
    
    % Diskretizacija bilinearnom transformacijom: nule i polovi analognog -> nule i polove digitalnog
    [bd,ad]=bilinear(ba,aa,fs); 
    
    [hd,wd] = freqz(bd,ad,Nfreqz);
    Hd = abs(hd);
      
   % provera gabarita
    
    df=(fs/2)/Nfreqz;

    NPOK = 0;
    POOK = 0;

    ia=ceil(fa/df)+1;
    ip=floor(fp/df)+1;
   
    Hp=Hd(ip:length(Hd));
    Ha=Hd(1:ia);

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