function Spec = BulkSpectrum(EpsS1, EpsS2, R, L, A, Hmin, Hmin0, EpsT1, EpsT2, phi, const, n, om_min, om_max, Num, FF, flag, Path)

Eps=EpsS1+1i*EpsS2;
EpsT=EpsT1+1i*EpsT2;

Spec=zeros(Num,5);

a=L/2;
c=sqrt(a*(a-R));
H0=Hmin0+A;

EpsRef=Eps;

M=Num;
om=linspace(om_min,om_max,M);

switch flag
   case 1
        Epsilon=Path;
        fMat=Epsilon(:,1);
        Eps_Mat=Epsilon(:,2)+1i*Epsilon(:,3);
   case 2
        EpsRM=load('real_eps_PMMA.txt');
        EpsIM=load('imag_eps_PMMA.txt');
        fMat=EpsRM(:,1);
        Eps_Mat=EpsRM(:,2)+1i*EpsIM(:,2);
    case 3
        fMat=om;
        NR=7;
        Eps_par=1;
        Eps_per=1;
        wLO_par=[1223, 1222, 794, 540, 514, 559, 413];
        wTO_par=[1220, 1080, 778, 539, 509, 495, 364];
        gam_par=[183, 7.5, 78, 22, 7.1, 4.5, 5.1];

        wLO_per=[1229, 1165, 1215, 815, 700, 522, 421];
        wTO_per=[1227, 1163, 1072, 797, 697, 450, 394];
        gam_per=[135, 7, 7.6, 7.2, 8.4, 4, 2.8];

        for i = 1:NR
            Eps_par=Eps_par+(wLO_par(i)^2-wTO_par(i)^2)./(wTO_par(i)^2-om.^2-1i*gam_par(i)*om);
            Eps_per=Eps_per+(wLO_per(i)^2-wTO_per(i)^2)./(wTO_per(i)^2-om.^2-1i*gam_per(i)*om);
        end
        Eps_par=Eps_par*2.383;
        Eps_per=Eps_per*2.356;
        Eps_Mat=zeros(1,M);
        for q = 1:M
            if (imag(sqrt(Eps_par(q)*Eps_per(q)))<0)
                Eps_Mat(q)=-sqrt(Eps_par(q)*Eps_per(q));
            else
                Eps_Mat(q)=sqrt(Eps_par(q)*Eps_per(q));
            end
        end
    case 4
        fMat=om;

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
        Eps_Mat=zeros(1,M);
        Eps_par=zeros(1,M);
        Eps_per=zeros(1,M);
        for q = 1:M
            Eps_par(q)=Eb_par*(1+(wLO_par^2-wTO_par^2)/(wTO_par^2-om(q)^2-1i*gam_par*om(q))-wp_par^2/(om(q)^2+1i*om(q)*gamp_par));
            Eps_per(q)=Eb_per*(1+(wLO_per^2-wTO_per^2)/(wTO_per^2-om(q)^2-1i*gam_per*om(q))-wp_per^2/(om(q)^2+1i*om(q)*gamp_per));
            if (imag(sqrt(Eps_par(q)*Eps_per(q)))<0)
                Eps_Mat(q)=-sqrt(Eps_par(q)*Eps_per(q));
            else
                Eps_Mat(q)=sqrt(Eps_par(q)*Eps_per(q));
            end
        end
    case 5
        fMat=om;

        Eb_par=6.72;
        Eb_per=6.56;
        wTO_par=788;
        wTO_per=797;
        wLO_par=964;
        wLO_per=970;
        gam_par=5.5;
        gam_per=5.9;
        wp_par=120;
        wp_per=230;
        gamp_par=250;
        gamp_per=500;
        Eps_Mat=zeros(1,M);
        Eps_par=zeros(1,M);
        Eps_per=zeros(1,M);
        for q = 1:M
            Eps_par(q)=Eb_par*(1+(wLO_par^2-wTO_par^2)/(wTO_par^2-om(q)^2-1i*gam_par*om(q))-wp_par^2/(om(q)^2+1i*om(q)*gamp_par));
            Eps_per(q)=Eb_per*(1+(wLO_per^2-wTO_per^2)/(wTO_per^2-om(q)^2-1i*gam_per*om(q))-wp_per^2/(om(q)^2+1i*om(q)*gamp_per));
            if (imag(sqrt(Eps_par(q)*Eps_per(q)))<0)
                Eps_Mat(q)=-sqrt(Eps_par(q)*Eps_per(q));
            else
                Eps_Mat(q)=sqrt(Eps_par(q)*Eps_per(q));
            end
        end
end

N=10;
m=10;
J=zeros(m,m,N); C=zeros(m,1); H=zeros(1,N); psi=zeros(1,N);
pz_PMMA=zeros(1,N); pzn_PMMA=zeros(1,M); Phase=zeros(1,M);

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
pzn_sub=0; pz_sub=zeros(1,N);
Eps2=EpsRef;
Eps1=1;
kz1=sqrt(Eps1)*cos(phi);
kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
rp=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
E=(1+const*rp)^(2*FF);
beta=(Eps2-Eps1)/(Eps2+Eps1);
for k = 1:(N+1)
    Mat=eye(m);
    for i = 1:m
        C(i,1)=-beta*(EpsT-1)^2*c*E/4*J(1,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2))/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
        for j = 1:m
            Mat(i,j)=Mat(i,j)-beta*(EpsT-1)*(2*j+1)/2*J(j,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
        end
    end
    U=Mat^(-1)*C;
    pz_sub(k)=2*a*c*U(1)-(EpsT-1)*a*c^2*E/3/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2));
    if (k==1 || k==(N+1))
        pzn_sub=pzn_sub+pz_sub(k)*exp(-1i*pi*(k-1)*n/N)/(2*N);
    end
    if (k>1 && k<(N+1))
        pzn_sub=pzn_sub+pz_sub(k)*cos(pi*(k-1)*n/N)/N;
    end
end

H0=Hmin+A;
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

for q = 1:M
    Eps2=interp1(fMat,real(Eps_Mat),om(q),'spline')+1i*interp1(fMat,imag(Eps_Mat),om(q),'spline');
    beta=(Eps2-Eps1)/(Eps2+Eps1);
    kz1=sqrt(Eps1)*cos(phi);
    kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
    rp=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
    E=(1+const*rp)^(2*FF);
    for k = 1:(N+1)
        Mat=eye(m);
        for i = 1:m
            C(i,1)=-beta*(EpsT-1)^2*c*E/4*J(1,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2))/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
            for j = 1:m
                Mat(i,j)=Mat(i,j)-beta*(EpsT-1)*(2*j+1)/2*J(j,i,k)/legendrePf(i,a/c)/(EpsT*legendreQ(i,a/c)/legendrePf(i,a/c)-(a/c*legendreQ(i,a/c)-legendreQ(i-1,a/c))/(a/c*legendrePf(i,a/c)-legendrePf(i-1,a/c)));
            end
        end
        U=Mat^(-1)*C;
        pz_PMMA(k)=2*a*c*U(1)-(EpsT-1)*a*c^2*E/3/(EpsT*legendreQ(1,a/c)-a/c*legendreQ(0,a/c)+a^2/(a^2-c^2));
        if (k==1 || k==(N+1))
            pzn_PMMA(q)=pzn_PMMA(q)+pz_PMMA(k)*exp(-1i*pi*(k-1)*n/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pzn_PMMA(q)=pzn_PMMA(q)+pz_PMMA(k)*cos(pi*(k-1)*n/N)/N;
        end
    end
    Phase(q)=angle(pzn_PMMA(q));
    if (Phase(q)<0)
        Phase(q)=Phase(q)+2*pi;
    end
    Spec(q,1)=om(q);
    Spec(q,2)=abs(pzn_PMMA(q)/pzn_sub);
    Spec(q,3)=Phase(q)-pi;
    Spec(q,4)=real(Eps2);
    Spec(q,5)=imag(Eps2);
end

for i = 1:(M-1)
    Delta=Spec(i,3)-Spec(i+1,3);
    if (abs(Delta)>3)
        for k = (i+1):M
            Spec(k,3)=Spec(k,3)+sign(Delta)*2*pi;
        end
    end
end

end

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
