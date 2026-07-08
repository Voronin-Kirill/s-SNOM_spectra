function Spec = DoubleLayerSpectrum(EpsS1, EpsS2, R, L, A, Hmin, Hmin0, EpsT1, EpsT2, phi, const, n, om_min, om_max, Num, FF, flag1, Path_layer1, d1, flag2, Path_layer2, d2, flag3, Path_substr, N_acc, M_acc, K)

Eps=EpsS1+1i*EpsS2;
EpsT=EpsT1+1i*EpsT2;

Spec=zeros(Num,7);

a=L/2;
c=sqrt(a*(a-R));

M=Num;
om=linspace(om_min,om_max,M);
wl=1./om*10^(7);
ko=2*pi./wl;
Eps0=8.85418782*10^(-12);   % Vacuum permittivity
E0=1*4*pi*Eps0;             % External electric field

switch flag1
   case 1
        EpsilonL=Path_layer1;
        fMat1=EpsilonL(:,1);
        Eps_Mat1=EpsilonL(:,2)+1i*EpsilonL(:,3);
   case 2
        EpsRM=load('real_eps_PMMA.txt');
        EpsIM=load('imag_eps_PMMA.txt');
        fMat1=EpsRM(:,1);
        Eps_Mat1=EpsRM(:,2)+1i*EpsIM(:,2);
    case 3
        Eps_SiO2=load('Eps_SiO2.txt');
        fMat1=Eps_SiO2(:,1);
        Eps_Mat1=Eps_SiO2(:,2)+1i*Eps_SiO2(:,3);
    case 4
        fMat1=om;
        Einf=6.6;
        wTO=797;
        wLO=973;
        gam=6;
        Eps_Mat1=Einf*(om.^2-wLO^2+1i*gam*om)./(om.^2-wTO^2+1i*gam*om);
end

switch flag2
   case 1
        EpsilonL=Path_layer2;
        fMat2=EpsilonL(:,1);
        Eps_Mat2=EpsilonL(:,2)+1i*EpsilonL(:,3);
   case 2
        EpsRM=load('real_eps_PMMA.txt');
        EpsIM=load('imag_eps_PMMA.txt');
        fMat2=EpsRM(:,1);
        Eps_Mat2=EpsRM(:,2)+1i*EpsIM(:,2);
    case 3
        Eps_SiO2=load('Eps_SiO2.txt');
        fMat2=Eps_SiO2(:,1);
        Eps_Mat2=Eps_SiO2(:,2)+1i*Eps_SiO2(:,3);
    case 4
        fMat2=om;
        Einf=6.6;
        wTO=797;
        wLO=973;
        gam=6;
        Eps_Mat2=Einf*(om.^2-wLO^2+1i*gam*om)./(om.^2-wTO^2+1i*gam*om);
end

switch flag3
   case 1
        EpsilonS=Path_substr;
        fSub=EpsilonS(:,1);
        Eps_Sub=EpsilonS(:,2)+1i*EpsilonS(:,3);
   case 2
        fSub=om;
        Eps_Sub=ones(1, M)*12;
    case 3
        Eps_SiO2=load('Eps_SiO2.txt');
        fSub=Eps_SiO2(:,1);
        Eps_Sub=Eps_SiO2(:,2)+1i*Eps_SiO2(:,3);
    case 4
        fSub=om;
        Einf=6.6;
        wTO=797;
        wLO=973;
        gam=6;
        Eps_Sub=Einf*(om.^2-wLO^2+1i*gam*om)./(om.^2-wTO^2+1i*gam*om);
end

m=30;
N=30;
quadN = 48;

psi=zeros(1,N+1);
H=zeros(1,N+1);
J=zeros(m,m,N+1);

% Reference material

Eps1=1;
H0=Hmin0+A;
Eps2=Eps;
beta=(Eps2-Eps1)/(Eps2+Eps1);
kz1=sqrt(Eps1)*cos(phi);
kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
rp=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
E=E0*(1+const*rp)^(2*FF);                     % Far-field reflection coefficient

% Fixed quadrature nodes on [-1,1]
[xq,wq] = gaussLegendreRule(quadN);
xq = xq(:).';
wq = wq(:).';

% P_i(xq), i=1..m, used for all heights
P_x_all = legendreP_all_stable(m,xq);
Px = P_x_all(2:m+1,:);     % rows correspond to degrees 1..m

alpha = a/c;
P_alpha_all = legendreP_all_stable(m,alpha);
Q_alpha_all = legendreQ_all_stable(m,alpha);

P_i_alpha = P_alpha_all(2:m+1).';
Q_i_alpha = Q_alpha_all(2:m+1).';
Q_im1_alpha = Q_alpha_all(1:m).';
P_im1_alpha = P_alpha_all(1:m).';

den_i = EpsT.*Q_i_alpha./P_i_alpha -(alpha.*Q_i_alpha - Q_im1_alpha) ./ (alpha.*P_i_alpha - P_im1_alpha);

row_scale = 1 ./ (P_i_alpha .* den_i);
col_scale = ((2*(1:m)+1)/2).';

common_C_denom = EpsT*Q_alpha_all(2) - alpha*Q_alpha_all(1) + a^2/(a^2-c^2);
const_C = (EpsT-1)^2*c*E/4 / common_C_denom;
const_pz = (EpsT-1)*a*c^2*E/3 / common_C_denom;

% Integral tables.

for k = 1:(N+1)
    psi(k)=pi*(k-1)/N;
    H(k)=H0-A*cos(psi(k));

    J(:,:,k) = integral_table_for_shift(Px, xq, wq, m, a, c, H(k), 0);
end

pzn_ref=0; pz=zeros(1,N+1);
for k = 1:(N+1)
        Wtab = -beta * J(:,:,k);

        C = const_C * (Wtab(1,:).' .* row_scale);
        Mat = eye(m) + (EpsT-1) * ((Wtab.' .* col_scale.') .* row_scale);

        U = Mat\C;
        pz(k)=2*a*c*U(1)-const_pz;

        % Harmonic calculation.
        if (k==1 || k==(N+1))
            pzn_ref=pzn_ref+pz(k)*exp(-1i*pi*(k-1)*n/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pzn_ref=pzn_ref+pz(k)*cos(pi*(k-1)*n/N)/N;
        end
end

% Sample

N=N_acc;
m=M_acc;
quadN = 128;

psi=zeros(1,N+1);
H=zeros(1,N+1);
J=zeros(m,m,N+1);
I=zeros(m,m,N+1,K);

H0=A+Hmin;
% Fixed quadrature nodes on [-1,1]
[xq,wq] = gaussLegendreRule(quadN);
xq = xq(:).';
wq = wq(:).';

% P_i(xq), i=1..m, used for all heights
P_x_all = legendreP_all_stable(m,xq);
Px = P_x_all(2:m+1,:);     % rows correspond to degrees 1..m

alpha = a/c;
P_alpha_all = legendreP_all_stable(m,alpha);
Q_alpha_all = legendreQ_all_stable(m,alpha);

P_i_alpha = P_alpha_all(2:m+1).';
Q_i_alpha = Q_alpha_all(2:m+1).';
Q_im1_alpha = Q_alpha_all(1:m).';
P_im1_alpha = P_alpha_all(1:m).';

den_i = EpsT.*Q_i_alpha./P_i_alpha -(alpha.*Q_i_alpha - Q_im1_alpha) ./ (alpha.*P_i_alpha - P_im1_alpha);

row_scale = 1 ./ (P_i_alpha .* den_i);
col_scale = ((2*(1:m)+1)/2).';

common_C_denom = EpsT*Q_alpha_all(2) - alpha*Q_alpha_all(1) + a^2/(a^2-c^2);
const_C = (EpsT-1)^2*c/4 / common_C_denom;
const_pz = (EpsT-1)*a*c^2/3 / common_C_denom;

pzn=zeros(1,M);
Phasen=zeros(1,M);

% Integral tables.
di=[2*d1, 2*d1+2*min(d1,d2)*linspace(1,K-1,K-1)];
for k = 1:(N+1)
    psi(k)=pi*(k-1)/N;
    H(k)=H0-A*cos(psi(k));

    J(:,:,k) = integral_table_for_shift(Px, xq, wq, m, a, c, H(k), 0);

    for l = 1:K
        I(:,:,k,l) = integral_table_for_shift(Px, xq, wq, m, a, c, H(k), di(l));
    end
end

for i = 1:M
    Eps4=interp1(fSub,real(Eps_Sub),om(i),'spline')+1i*interp1(fSub,imag(Eps_Sub),om(i),'spline');
    Eps3=interp1(fMat2,real(Eps_Mat2),om(i),'spline')+1i*interp1(fMat2,imag(Eps_Mat2),om(i),'spline');
    Eps2=interp1(fMat1,real(Eps_Mat1),om(i),'spline')+1i*interp1(fMat1,imag(Eps_Mat1),om(i),'spline');
    beta_1=(Eps2-Eps1)/(Eps2+Eps1);
    CA=charge(Eps1, Eps2, Eps3, Eps4, d1, d2, K);
    kz4=sqrt(Eps4-Eps1*(sin(phi))^2);
    kz3=sqrt(Eps3-Eps1*(sin(phi))^2);
    kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
    kz1=sqrt(Eps1)*cos(phi);
    rp1=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
    rp2=(Eps3/kz3-Eps2/kz2)/(Eps3/kz3+Eps2/kz2);
    rp3=(Eps4/kz4-Eps3/kz3)/(Eps4/kz4+Eps3/kz3);
    rp=(rp2+rp3*exp(2i*kz3*d2*ko(i)))/(1+rp2*rp3*exp(2i*kz3*d2*ko(i)));
    rp=(rp1+rp*exp(2i*kz2*d1*ko(i)))/(1+rp1*rp*exp(2i*kz2*d1*ko(i)));
    E=E0*(1+const*rp)^(2*FF);                         % Far-field reflection coefficient

    pz=zeros(1,N+1);
    for k = 1:(N+1)
        Wtab = -beta_1 * J(:,:,k);
        for l = 1:K
            Wtab = Wtab + CA(l)*I(:,:,k,l);
        end

        C = const_C*E*(Wtab(1,:).' .* row_scale);
        Mat = eye(m) + (EpsT-1) * ((Wtab.' .* col_scale.') .* row_scale);

        U = Mat\C;
        pz(k)=2*a*c*U(1)-const_pz*E;

        % Harmonic calculation
        if (k==1 || k==(N+1))
            pzn(i)=pzn(i)+pz(k)*exp(-1i*pi*(k-1)*n/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pzn(i)=pzn(i)+pz(k)*cos(pi*(k-1)*n/N)/N;
        end
    end
    Phasen(i)=angle(pzn(i)/pzn_ref);
    % Phase correction (making phase positive)
    if (Phasen(i)<-pi/4)
        Phasen(i)=Phasen(i)+2*pi;
    end
    Spec(i,1)=om(i);
    Spec(i,2)=abs(pzn(i)/pzn_ref);
    Spec(i,3)=Phasen(i);
    Spec(i,4)=real(Eps2);
    Spec(i,5)=imag(Eps2);
    Spec(i,6)=real(Eps3);
    Spec(i,7)=imag(Eps3);
    Spec(i,8)=real(Eps4);
    Spec(i,9)=imag(Eps4);
end

% Phase correction (removing phase jumps)

for i = 1:(M-1)
    Delta=Spec(i,3)-Spec(i+1,3);
    if (abs(Delta)>3)
        for k = (i+1):M
            Spec(k,3)=Spec(k,3)+sign(Delta)*2*pi;
        end
    end
end

end

% Local functions

function T = integral_table_for_shift(Px, x, w, m, a, c, H, imageShift)
    % Computes T(i,j) = int_{-1}^1 P_i(x)*P_j(y(x))*Q_j(eta(x)) dx
    % for i,j = 1..m, using fixed quadrature.

    s = 2*a/c + 2*H/c + imageShift/c - a/c*x;
    Aexpr = 1 + s.^2 + (a^2/c^2 - 1)*(1 - x.^2);
    rad = sqrt(Aexpr.^2 - 4*s.^2);
    eta = sqrt((Aexpr + rad)/2);
    y = s ./ eta;

    P_y_all = legendreP_all_stable(m,y);
    Q_eta_all = legendreQ_all_stable(m,eta);

    Py = P_y_all(2:m+1,:);       % degrees 1..m
    Qe = Q_eta_all(2:m+1,:);     % degrees 1..m

    B = Py .* Qe .* w;           % implicit expansion over rows
    T = Px * B.';
end

function P = legendreP_all_stable(n,x)
    % P(degree+1, pointIndex) = P_degree(x)
    x = x(:).';
    P = zeros(n+1,numel(x));
    P(1,:) = 1;
    if n >= 1
        P(2,:) = x;
    end
    for k = 1:n-1
        P(k+2,:) = ((2*k+1)/(k+1))*x.*P(k+1,:) - (k/(k+1))*P(k,:);
    end
end

function Q = legendreQ_all_stable(n,x)
    % Stable evaluation of Legendre functions Q_0..Q_n for x > 1.
    % For degrees >=2 this follows the same continued-fraction/Wronskian
    % construction as the user's original legendreQ.m, vectorized over x.

    x = x(:).';
    Q = zeros(n+1,numel(x));

    Q(1,:) = 0.5*log(abs((x+1)./(x-1)));
    if n >= 1
        Q(2,:) = x.*Q(1,:) - 1;
    end

    if n <= 1
        return;
    end

    P = legendreP_all_stable(n,x);
    for degree = 2:n
        pn = P(degree+1,:);
        p0 = P(degree,:);        % P_{degree-1}
        Hn = continued_fraction_Hn(degree,x);
        Q(degree+1,:) = Hn/degree ./ (pn - Hn.*p0);
    end
end

function Hn = continued_fraction_Hn(n,x)
    % Vectorized version of the continued fraction part of legendreQ.m.
    esp = 1e-15;
    tiny = 1e-30;
    maxIter = 10000;

    Hn = tiny + zeros(size(x));
    Cj1 = Hn;
    Dj1 = zeros(size(x));

    aj = 1;
    bj = (2 + 1/n) * x;

    Dj = bj + aj*Dj1;
    Dj(abs(Dj) < tiny) = tiny;

    Cj = bj + aj./Cj1;
    Cj(abs(Cj) < tiny) = tiny;

    Dj = 1./Dj;
    Del = Cj.*Dj;
    Hn = Hn.*Del;
    Dj1 = Dj;
    Cj1 = Cj;

    converged = abs(Del - 1) < esp;
    j = 0;

    while any(~converged) && j < maxIter
        aj = -(1 + 1/(n+j));
        bj = (2 + 1/(n+j+1)) * x;

        Dj = bj + aj*Dj1;
        Dj(abs(Dj) < tiny) = tiny;

        Cj = bj + aj./Cj1;
        Cj(abs(Cj) < tiny) = tiny;

        Dj = 1./Dj;
        Del = Cj.*Dj;
        Hn = Hn.*Del;

        converged = converged | (abs(Del - 1) < esp);
        Dj1 = Dj;
        Cj1 = Cj;
        j = j + 1;
    end

    if j >= maxIter
        warning('continued_fraction_Hn:notConverged', ...
            'Continued fraction did not converge for n=%d at all quadrature nodes.', n);
    end
end

function [x,w] = gaussLegendreRule(n)
    % Golub-Welsch nodes and weights for Gauss-Legendre quadrature on [-1,1].
    beta = (1:n-1) ./ sqrt(4*(1:n-1).^2 - 1);
    J = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);
    [x,idx] = sort(diag(D));
    V = V(:,idx);
    w = 2*(V(1,:).^2).';
end

function C = charge(Eps1, Eps2, Eps3, Eps4, d1, d2, M)

    beta1=(Eps2-Eps1)/(Eps2+Eps1);
    beta2=(Eps3-Eps2)/(Eps3+Eps2);
    beta3 = (Eps4-Eps3)/(Eps4+Eps3);
    di=[2*d1, 2*d1+2*min(d1,d2)*linspace(1,M-1,M-1)];
    
    N = 1000;
    kn = linspace(0, 10/min(d1,d2), N);
    dk = kn(2) - kn(1);
    rn = ref(beta1, beta2, beta3, d1, d2, kn);
    
    z=exp(-di*dk);
    
    Mat = zeros(M);
    F = zeros(M,1);
    
    for i = 1:N
        F(i) = rn(i);
        for j = 1:M
            Mat(i,j) = z(j)^(i-1);
        end
    end

    C = Mat\F;

end

function r=ref(beta1, beta2, beta3, d1, d2, k)
    r1=-(beta2 + beta3.*exp(-2*k*d2))./(1 + beta2.*beta3.*exp(-2*k*d2));
    r =beta1-(beta1-r1.*exp(-2*k*d1))./(1 - beta1.*r1.*exp(-2*k*d1));
end
