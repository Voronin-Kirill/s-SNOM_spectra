clear
format long

Eps0=8.85418782*10^(-12);   % Vacuum permittivity

R=25;                       % Spheroid apex curvature radius

a=300;                      % a=L/2 is the spheroid major semi-axis
c=sqrt(a*(a-R));

HminR=2;                    % Minimum shperoid-reference material distance
HminS=2;                    % Minimum shperoid-sample distance
A=30;                       % Oscillation amplitude

H0=A+HminR;

E0=1*4*pi*Eps0;             % External electric field
phi=45/180*pi;              % Angle of incidence
                   
EpsAu=-5000+800i;           % Dielectric permittivity of gold (used as a reference material)
EpsT=-10^10;                % Dielectric permittivity of the tip

% 4HSiC Lorentz parameters
Eb_par=6.78;
Eb_per=6.56;
wTO_par=782;
wTO_per=797;
wLO_par=967;
wLO_per=971;
gam_par=3.3;
gam_per=3.3;
wp_par=220;
wp_per=275;
gamp_par=450;
gamp_per=450;

% 6HSiC Lorentz parameters
% Eb_par=6.72;
% Eb_per=6.56;
% wTO_par=788;
% wTO_per=797;
% wLO_par=964;
% wLO_per=970;
% gam_par=5.5;
% gam_per=5.9;
% wp_par=120;
% wp_per=230;
% gamp_par=250;
% gamp_per=500;

% SiO2 Lorentz parameters
% NR=7;
% wLO_par=[1223, 1222, 794, 540, 514, 559, 413];
% wTO_par=[1220, 1080, 778, 539, 509, 495, 364];
% gam_par=[183, 7.5, 78, 22, 7.1, 4.5, 5.1];
% 
% wLO_per=[1229, 1165, 1215, 815, 700, 522, 421];
% wTO_per=[1227, 1163, 1072, 797, 697, 450, 394];
% gam_per=[135, 7, 7.6, 7.2, 8.4, 4, 2.8];

M=101;                      % Number of point over the frequency range
om=linspace(850,1050,M);    % Frequency range

N=50;                       % Spheroid oscillation discretization (number of spheroid-sample distances)
m=20;                       % Matrix size

% Arrays initialization
J=zeros(m,m,N); C=zeros(m,1); H=zeros(1,N); psi=zeros(1,N);
pz_ref=zeros(1,N); pz_sample=zeros(1,N);
pz2_sample=zeros(1,M); pz3_sample=zeros(1,M); pz4_sample=zeros(1,M);
Phase2=zeros(1,M); Phase3=zeros(1,M); Phase4=zeros(1,M);

% Integrals calculation
for k = 1:(N+1)
    psi(k)=pi*(k-1)/N;
    H(k)=H0-A*cos(psi(k));
    for i = 1:m
        for j = 1:m
            f=@(x) legendrePf(i,x).*legendrePf(j,(2*a/c+2*H(k)/c-a/c*x)./sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)+sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)).^2-4*(2*a/c+2*H(k)/c-a/c*x).^2))/2)).*legendreQ(j,sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)+sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)).^2-4*(2*a/c+2*H(k)/c-a/c*x).^2))/2));
            J(i,j,k)=integral(f,-1,1);
        end
    end
end

pz2_ref=0; pz3_ref=0; pz4_ref=0;
Eps2=EpsAu;
% Choose include or not far-feild coefficient
% kz1=cos(phi);
% kz2=sqrt(Eps2-(sin(phi))^2);
% rp=(Eps2/kz2-1/kz1)/(Eps2/kz2+1/kz1);
% E=E0*(1+rp)^2*sin(phi);
E=E0;
beta=(Eps2-1)/(Eps2+1);

% Reference material response
for k = 1:(N+1)
    % Matrix and free term creation
    Mat=eye(m);
    for i = 1:m
        C(i,1)=-beta*(EpsT-1)^2*c*E/4*J(1,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2))/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
        for j = 1:m
            Mat(i,j)=Mat(i,j)-beta*(EpsT-1)*(2*j+1)/2*J(j,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
        end
    end
    % System solution
    U=Mat\C;
    % Dipole moment calculation
    pz_ref(k)=2*a*c*U(1)-(EpsT-1)*a*c^2*E/3/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2));
    % Fourier transorms
    if (k==1 || k==(N+1))
        pz2_ref=pz2_ref+pz_ref(k)*exp(-1i*pi*(k-1)*2/N)/(2*N);
    end
    if (k>1 && k<(N+1))
        pz2_ref=pz2_ref+pz_ref(k)*cos(pi*(k-1)*2/N)/N;
    end
    if (k==1 || k==(N+1))
        pz3_ref=pz3_ref+pz_ref(k)*exp(-1i*pi*(k-1)*3/N)/(2*N);
    end
    if (k>1 && k<(N+1))
        pz3_ref=pz3_ref+pz_ref(k)*cos(pi*(k-1)*3/N)/N;
    end
    if (k==1 || k==(N+1))
        pz4_ref=pz4_ref+pz_ref(k)*exp(-1i*pi*(k-1)*4/N)/(2*N);
    end
    if (k>1 && k<(N+1))
        pz4_ref=pz4_ref+pz_ref(k)*cos(pi*(k-1)*4/N)/N;
    end
end

% Integrals for the sample (needed only if HminS and HminR are different)

H0=HminS+A;
for k = 1:(N+1)
    psi(k)=pi*(k-1)/N;
    H(k)=H0-A*cos(psi(k));
    for i = 1:m
        for j = 1:m
            f=@(x) legendrePf(i,x).*legendrePf(j,(2*a/c+2*H(k)/c-a/c*x)./sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)+sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)).^2-4*(2*a/c+2*H(k)/c-a/c*x).^2))/2)).*legendreQ(j,sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)+sqrt((1+(2*a/c+2*H(k)/c-a/c*x).^2+(a^2/c^2-1)*(1-x.^2)).^2-4*(2*a/c+2*H(k)/c-a/c*x).^2))/2));
            J(i,j,k)=integral(f,-1,1);
        end
    end
end

% Spectrum calculation
Eps_par=zeros(1,M); Eps_per=zeros(1,M); rp=zeros(1,M);
for q = 1:M
    % Dielectric permittivity of 4H/6HSiC
    Eps_per(q)=Eb_per*(1+(wLO_per^2-wTO_per^2)/(wTO_per^2-om(q)^2-1i*gam_per*om(q))-wp_per^2/(om(q)^2+1i*om(q)*gamp_per));
    Eps_par(q)=Eb_par*(1+(wLO_par^2-wTO_par^2)/(wTO_par^2-om(q)^2-1i*gam_par*om(q))-wp_par^2/(om(q)^2+1i*om(q)*gamp_par));

    % Dielectric permittivity of SiO2
    % Eps_par(q)=1;
    % Eps_per(q)=1;
    % for o = 1:NR
    % Eps_par(q)=Eps_par(q)+(wLO_par(o)^2-wTO_par(o)^2)/(wTO_par(o)^2-om(q)^2-1i*gam_par(o)*om(q));
    % Eps_per(q)=Eps_per(q)+(wLO_per(o)^2-wTO_per(o)^2)/(wTO_per(o)^2-om(q)^2-1i*gam_per(o)*om(q));
    % end
    % Eps_par(q)=Eps_par(q)*2.383;
    % Eps_per(q)=Eps_per(q)*2.356;
    
    if (imag(sqrt(Eps_par(q)*Eps_per(q)))<0)
        Eps2=-sqrt(Eps_par(q)*Eps_per(q));
    end
    if (imag(sqrt(Eps_par(q)*Eps_per(q)))>0)
        Eps2=sqrt(Eps_par(q)*Eps_per(q));
    end
    beta=(Eps2-1)/(Eps2+1);
    % Choose include or not far-feild coefficient
    % kz1=cos(phi);
    % kz2=sqrt(Eps_per(q)-Eps_per(q)/Eps_par(q)*(sin(phi))^2);
    % rp(q)=(Eps_per(q)/kz2-1/kz1)/(Eps_per(q)/kz2+1/kz1);
    % E=E0*(1+rp(q))^2*sin(phi);
    E=E0;
    for k = 1:(N+1)
        Mat=eye(m);
        for i = 1:m
            C(i,1)=-beta*(EpsT-1)^2*c*E/4*J(1,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2))/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
            for j = 1:m
                Mat(i,j)=Mat(i,j)-beta*(EpsT-1)*(2*j+1)/2*J(j,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
            end
        end
        U=Mat\C;
        pz_sample(k)=2*a*c*U(1)-(EpsT-1)*a*c^2*E/3/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2));
        if (k==1 || k==(N+1))
            pz2_sample(q)=pz2_sample(q)+pz_sample(k)*exp(-1i*pi*(k-1)*2/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz2_sample(q)=pz2_sample(q)+pz_sample(k)*cos(pi*(k-1)*2/N)/N;
        end
        if (k==1 || k==(N+1))
            pz3_sample(q)=pz3_sample(q)+pz_sample(k)*exp(-1i*pi*(k-1)*3/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz3_sample(q)=pz3_sample(q)+pz_sample(k)*cos(pi*(k-1)*3/N)/N;
        end
        if (k==1 || k==(N+1))
            pz4_sample(q)=pz4_sample(q)+pz_sample(k)*exp(-1i*pi*(k-1)*4/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz4_sample(q)=pz4_sample(q)+pz_sample(k)*cos(pi*(k-1)*4/N)/N;
        end
    end
    % Normalized phase
    Phase2(q)=angle(pz2_sample(q)/pz2_ref);
    Phase3(q)=angle(pz3_sample(q)/pz3_ref);
    Phase4(q)=angle(pz4_sample(q)/pz4_ref);
    if (Phase2(q)<0)
        Phase2(q)=Phase2(q)+2*pi;
    end
    if (Phase3(q)<0)
        Phase3(q)=Phase3(q)+2*pi;
    end
    if (Phase4(q)<0)
        Phase4(q)=Phase4(q)+2*pi;
    end
end

% Normalized amplitude
Ampl2=abs(pz2_sample/pz2_ref);
Ampl3=abs(pz3_sample/pz3_ref);
Ampl4=abs(pz4_sample/pz4_ref);

% Phase correction algorithm. To remove 2pi jumps

for i = 1:(M-1)
    Delta=Phase2(i)-Phase2(i+1);
    if (abs(Delta)>3)
        for k = (i+1):M
            Phase2(k)=Phase2(k)+sign(Delta)*2*pi;
        end
    end
end
for i = 1:(M-1)
    Delta=Phase3(i)-Phase3(i+1);
    if (abs(Delta)>3)
        for k = (i+1):M
            Phase3(k)=Phase3(k)+sign(Delta)*2*pi;
        end
    end
end
for i = 1:(M-1)
    Delta=Phase4(i)-Phase4(i+1);
    if (abs(Delta)>3)
        for k = (i+1):M
            Phase4(k)=Phase4(k)+sign(Delta)*2*pi;
        end
    end
end

% Plot the spectrum

figure;
plot(om, Ampl2, 'red', om, Ampl3, 'green', om, Ampl4, 'blue');
xlabel('\omega, cm^{-1}');
ylabel('|\sigma_n|');
title('Amplitude spectrum');
legend('n=2', 'n=3', 'n=4');

figure;
plot(om, Phase2, 'red', om, Phase3, 'green', om, Phase4, 'blue');
xlabel('\omega, cm^{-1}');
ylabel('arg(\sigma_n), rad');
title('Phase spectrum');
legend('n=2', 'n=3', 'n=4');

% Legendre functions calculation

function pn=legendrePf(n,x)

if (n==0)
    pn=1;
end

if (n==1)
    pn=x;
end

p0=1;
p1=x;

if (n>1)
    for i = 1:n-1
        pn=(2*i+1)/(i+1)*x.*p1-i/(i+1)*p0;
        p0=p1;
        p1=pn;
    end
end

end

function qn=legendreQ(n,x)

if (n==0)
    qn=log(abs((x+1)./(x-1)))/2;
end

if (n==1)
    qn=log(abs((x+1)./(x-1))).*x/2-1;
end

p0=1;
p1=x;

if (n>1)
    for i = 1:n-1
        pn=(2*i+1)/(i+1)*x.*p1-i/(i+1)*p0;
        p0=p1;
        p1=pn;
    end
    esp=10^(-15);
    t=10^(-30);
    Hn=t;
    Cj1=Hn;
    Dj1=0;
    f=1;
    j=0;
    aj=1;
    bj=(2+1/n)*x;
    Dj=bj+aj*Dj1;
    if (Dj==0)
        Dj=t;
    end
    Cj=bj+aj/Cj1;
    if (Cj==0)
        Cj=t;
    end
    Dj=1./Dj;
    Del=Cj.*Dj;
    Hn=Hn.*Del;
    Dj1=Dj;
    Cj1=Cj;
    while (f==1)
        aj=-(1+1/(n+j));
        bj=(2+1/(n+j+1))*x;
        Dj=bj+aj*Dj1;
        if (Dj==0)
            Dj=t;
        end
        Cj=bj+aj./Cj1;
        if (Cj==0)
            Cj=t;
        end
        Dj=1./Dj;
        Del=Cj.*Dj;
        Hn=Hn.*Del;
        if (abs(Del-1)<esp)
            f=0;
        end
        Dj1=Dj;
        Cj1=Cj;
        j=j+1;
    end
    qn=Hn/n./(pn-Hn.*p0);
end

end

