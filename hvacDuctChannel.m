function channel = hvacDuctChannel(frequency,WGlen)
radius=0.305/2;
%radius=0.125/2;
%WGlen = WGlen*1e-2;

if str2double(frequency) == 60
%     freq = 58.2E9:2E6:60.2E9;
        freq = 59E9:1E8:60E9;
else
    freq = 2.4E9:0.3125E6:2.48E9;
end
Zo=50;
eta = 377;
l = 0.035;
%l = 0.001;

load('besDerZerMat.mat');
load('besZerMat.mat');

sigma = 10^6;
mu = 4*pi*10^-7;
R = ((2*pi*freq*mu)./(2*sigma)).^0.5;
lightSpeed = 2.98E8;
k = 2*pi*freq./lightSpeed;
antresTE = zeros(1,length(freq));
antresTM = zeros(1,length(freq));
antreacTE = zeros(1,length(freq));
antreacTM = zeros(1,length(freq));
antimpTE = zeros(1,length(freq));
antimpTM = zeros(1,length(freq));
channel = zeros(1,length(freq));

for fi=1:length(freq)
    n_TE1=[];
    m_TE1=[];
    fc_TE1=[];
    coWnTE1=[];
    for m=1:1001
        for n=1:1000
            fc_TE_temp=(lightSpeed/(2*pi*radius))*besDerZerMat(m,n);
            if fc_TE_temp <= freq(fi)
                n_TE1 = [n_TE1; n];
                m_TE1 = [m_TE1; m-1];
                fc_TE1 = [fc_TE1; fc_TE_temp];
                coWnTE1 = [coWnTE1; besDerZerMat(m,n)/radius];
            end
        end
    end
    n_TE = n_TE1;
    m_TE = m_TE1;
    fc_TE = fc_TE1; 
    coWnTE = coWnTE1; 
    n_TM1=[];
    m_TM1=[];
    fc_TM1=[];
    coWnTM1=[];
    for m=1:1001
        for n=1:1000
            fc_TM_temp=(lightSpeed/(2*pi*radius))*besZerMat(m,n);
            if fc_TM_temp <= freq(fi)
                n_TM1 = [n_TM1; n];
                m_TM1 = [m_TM1; m-1];
                fc_TM1 = [fc_TM1; fc_TM_temp];
                coWnTM1 = [coWnTM1; (besZerMat(m,n))/radius];
            end
        end
    end
    n_TM = n_TM1;
    m_TM = m_TM1;
    fc_TM = fc_TM1; 
    coWnTM = coWnTM1;
end
gammaTE = zeros(length(freq), length(n_TE));
gammaTM = zeros(length(freq), length(n_TM));
radresTE = zeros(length(freq), length(n_TE));
radresTM = zeros(length(freq), length(n_TM));
TEmodeimp = zeros(length(freq),length(n_TE));
TMmodeimp = zeros(length(freq),length(n_TM));

for fr=1:length(freq)
    for n = 1:length(n_TE)
        gTESq = (coWnTE(n))^2;
        betaTE = (sqrt(k(fr)^2 - gTESq));
        alphaTempTE = (R(fr)*k(fr))/(eta*radius*betaTE);
        alphaTE=8.85*(alphaTempTE*(((m_TE(n)^2)/(((radius*(coWnTE(n)))^2) ...
            - ((m_TE(n))^2)))+(((lightSpeed*coWnTE(n))/(2*pi*freq(fr)))^2)));
        gammaTE(fr,n)=alphaTE+1i*betaTE;
        constRadTE = eta*(k(fr))*((m_TE(n))^2)/(betaTE*pi*((sin(l*(k(fr))* ...
            pi/180))^2)*((besselj(m_TE(n),radius*coWnTE(n)))^4)*...
            ((((radius*coWnTE(n))^2)-((m_TE(n))^2))^2));
        radTEfunc=@(zeta) ((besselj(m_TE(n),radius.*zeta.*coWnTE(n)).*...
            sin((k(fr).*radius.*((l/radius)-1+zeta))*pi/180))./zeta);
        radTEnum = (integral(radTEfunc,1-(l/radius),1)).^2;
        radresTE(fr,n)= constRadTE*radTEnum;
        if imag(radresTE(fr,n)) ~= 0
            radresTE(fr,n) =0;
            gammaTE(fr,n) = 0;
        end
    end
    for n = 1:length(n_TM)
        gTMSq = (coWnTM(n))^2;
        betaTM = (sqrt(k(fr)^2 - gTMSq));%*1E-3;
        alphaconstTM = (R(fr)*k(fr))/(eta*radius*betaTM);
        alphaTM=8.85*alphaconstTM;
        gammaTM(fr,n)=alphaTM+1i*betaTM;
        constRadTM = eta*betaTM/(pi*((sin((k(fr))*l*pi/180))^2)*(k(fr))*...
            ((0.5*((besselj(m_TM(n)-1,radius*coWnTM(n)))-...
            (besselj(m_TM(n)+1,radius*coWnTM(n)))))^2));
        
        radTMfunc=@(zeta) 0.5*((besselj((m_TM(n))-1,zeta.*(coWnTM(n))*radius))...
            -(besselj(m_TM(n)+1,zeta.*(coWnTM(n))*radius))).*...
            sin((k(fr).*radius.*((l./radius)-1+zeta))*pi/180);
        radTMnum =(integral(radTMfunc,1-(l./radius),1)).^2;
        radresTM(fr,n)= constRadTM*radTMnum;
        if imag(radresTM(fr,n)) ~= 0
            radresTM(fr,n) =0;
            gammaTM(fr,n) = 0;
        end
    end
end

radreacTE = imag(hilbert(radresTE));
radreacTM = imag(hilbert(radresTM));

for fi=1:length(freq)
    antresTE(fi) = sum(radresTE(fi,:));
    antresTM(fi) = sum(radresTM(fi,:));
    antreacTE(fi) = sum(radreacTE(fi,:));
    antreacTM(fi) = sum(radreacTM(fi,:));
    antimpTE(fi) = antresTE(fi) +1i*antreacTE(fi);
    antimpTM(fi) = antresTM(fi) +1i*antreacTM(fi);
    for n = 1:length(n_TE)
        TEmodeimp(fi,n)=radresTE(fi,n)+1i*radreacTE(fi,n);
    end
    for n = 1:length(n_TM)
        TMmodeimp(fi,n)=radresTM(fi,n)+1i*radreacTM(fi,n);
    end   
    TsTE = diag(exp(-1*gammaTE(fi,:)*WGlen));
    chTEmode = TEmodeimp(fi,:)*TsTE;
    TsTM = diag(exp(-1*gammaTM(fi,:)*WGlen));
    chTMmode = TMmodeimp(fi,:)*TsTM;
    channel(fi) = ((2*Zo)/(abs(antimpTM(fi) + Zo + antimpTE(fi))^2))*...
        (sum(chTEmode)+sum(chTMmode));
    if isnan(channel(fi)) 
        channel(fi) = 0;
    end
end

ghu = (ifft(fftshift(channel)))./sum(ifft(fftshift(channel)));
%[bb,aa] = invfreqs(channel,freq,4,5);

