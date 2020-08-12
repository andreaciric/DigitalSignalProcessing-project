function [bd, ad] = baseline_drift_filter(fs, fa, fp, Aa, Ap)

    T = 1/fs;
    Nfreqz = 15000; 
        
    Aav = Aa;
    Apv = Ap;
    
    Wpd=2*pi*fp/fs; 
    Wad=2*pi*fa/fs; 
    
    waa = 2/T*tan(Wad/2); 
    wpa = 2/T*tan(Wpd/2); 
    
    wpp = 1; 
    wap = wpa/waa; 

while(1)    
    
    D = (10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
    k = wpp/wap;
    
    N = ceil(acosh(sqrt(D))/acosh(1/k));    
        
    [z,p,k]=cheb2ap(N,Aav);

    baN=k*poly(z);
    aaN=poly(p);
    
    [ba,aa]=lp2hp(baN,aaN,waa);
    
    [bd,ad]=bilinear(ba,aa,fs); 
    
    [hd,wd] = freqz(bd,ad,Nfreqz);
    Hd = abs(hd);
 
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