clear
format long

R=20;                       % Spheroid apex curvature radius
a=300;                     % a=L/2 is the spheroid major semi-axis
c=sqrt(a*(a-R));

HminR=2;                    % Minimum spheroid-reference material distance
HminS=2;                    % Minimum spheroid-sample distance
A=80;                       % Oscillation amplitude

Eps0=8.85418782*10^(-12);   % Vacuum permittivity
E0=1*4*pi*Eps0;             % External electric field
phi=30/180*pi;              % Angle of incidence
zet=1;                      % Far-field factor

Eps1=1;                     % Dielectric permittivity of the air
EpsT=-10^10;                % Dielectric permittivity of the tip

M=501;
om=linspace(700,1200,M);    % Array of frequencies
wl=1./om*10^(7);
ko=2*pi./wl;

% Dielectric permittivity of the layer
Einf=6.6;
wTO=797;
wLO=973;
gam=6;
Eps=Einf*(om.^2-wLO^2+1i*gam*om)./(om.^2-wTO^2+1i*gam*om);

Epsi=15;

wp=845;
gamma=125;
EpsInAs=Epsi*(1-wp^2./(om.^2+1i*om*gamma));

N=30;                       % Height discretization
m=30;                       % Matrix size
quadN = 48;                 % Gauss-Legendre quadrature order; try 32, 48, 64, 96

psi=zeros(1,N+1);
H=zeros(1,N+1);
J=zeros(m,m,N+1);

% Reference material

Eps2=Epsi;
beta=(Eps2-Eps1)/(Eps2+Eps1);
kz1=sqrt(Eps1)*cos(phi);
kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
rp=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
E=E0*(1+zet*rp)^2*sin(phi);                     % Far-field reflection coefficient included
% E=E0;                                     % Far-field reflection coefficient not included

H0=A+HminR;
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

pz=zeros(1,N+1); pz2_ref=0; pz3_ref=0; pz4_ref=0;
for k = 1:(N+1)
        Wtab = -beta * J(:,:,k);

        C = const_C * (Wtab(1,:).' .* row_scale);
        Mat = eye(m) + (EpsT-1) * ((Wtab.' .* col_scale.') .* row_scale);

        U = Mat\C;
        pz(k)=2*a*c*U(1)-const_pz;

        % Harmonics calculation.
        if (k==1 || k==(N+1))
            pz2_ref=pz2_ref+pz(k)*exp(-1i*pi*(k-1)*2/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz2_ref=pz2_ref+pz(k)*cos(pi*(k-1)*2/N)/N;
        end
        if (k==1 || k==(N+1))
            pz3_ref=pz3_ref+pz(k)*exp(-1i*pi*(k-1)*3/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz3_ref=pz3_ref+pz(k)*cos(pi*(k-1)*3/N)/N;
        end
        if (k==1 || k==(N+1))
            pz4_ref=pz4_ref+pz(k)*exp(-1i*pi*(k-1)*4/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz4_ref=pz4_ref+pz(k)*cos(pi*(k-1)*4/N)/N;
        end
end

% Sample

H0=A+HminS;
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

pz2=zeros(1,M); pz3=zeros(1,M); pz4=zeros(1,M); rp=zeros(1,M); 
Phase2=zeros(1,M); Phase3=zeros(1,M); Phase4=zeros(1,M);

% Integral tables.

for k = 1:(N+1)
    psi(k)=pi*(k-1)/N;
    H(k)=H0-A*cos(psi(k));

    J(:,:,k) = integral_table_for_shift(Px, xq, wq, m, a, c, H(k), 0);
end

for i = 1:M
    Eps2=EpsInAs(i);
    beta=(Eps2-Eps1)/(Eps2+Eps1);
    kz1=sqrt(Eps1)*cos(phi);
    kz2=sqrt(Eps2-Eps1*(sin(phi))^2);
    rp(i)=(Eps2/kz2-Eps1/kz1)/(Eps2/kz2+Eps1/kz1);
    E=E0*(1+zet*rp(i))^2*sin(phi);                    % Far-field reflection coefficient included
    % E=E0;                                           % Far-field reflection coefficient not included

    pz=zeros(1,N+1);
    for k = 1:(N+1)
        Wtab = -beta * J(:,:,k);

        C = const_C*E*(Wtab(1,:).' .* row_scale);
        Mat = eye(m) + (EpsT-1) * ((Wtab.' .* col_scale.') .* row_scale);

        U = Mat\C;
        pz(k)=2*a*c*U(1)-const_pz*E;

        % Harmonics calculation.
        if (k==1 || k==(N+1))
            pz2(i)=pz2(i)+pz(k)*exp(-1i*pi*(k-1)*2/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz2(i)=pz2(i)+pz(k)*cos(pi*(k-1)*2/N)/N;
        end
        if (k==1 || k==(N+1))
            pz3(i)=pz3(i)+pz(k)*exp(-1i*pi*(k-1)*3/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz3(i)=pz3(i)+pz(k)*cos(pi*(k-1)*3/N)/N;
        end
        if (k==1 || k==(N+1))
            pz4(i)=pz4(i)+pz(k)*exp(-1i*pi*(k-1)*4/N)/(2*N);
        end
        if (k>1 && k<(N+1))
            pz4(i)=pz4(i)+pz(k)*cos(pi*(k-1)*4/N)/N;
        end
    end
    Phase2(i)=angle(pz2(i)/pz2_ref);
    Phase3(i)=angle(pz3(i)/pz3_ref);
    Phase4(i)=angle(pz4(i)/pz4_ref);

    % Phase correction (making phase positive)
    if (Phase2(i)<0)
        Phase2(i)=Phase2(i)+2*pi;
    end
    if (Phase3(i)<0)
        Phase3(i)=Phase3(i)+2*pi;
    end
    if (Phase4(i)<0)
        Phase4(i)=Phase4(i)+2*pi;
    end
end

Ampl2=abs(pz2/pz2_ref);
Ampl3=abs(pz3/pz3_ref);
Ampl4=abs(pz4/pz4_ref);

% Phase correction (removing phase jumps)

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

figure;
plot(om, Ampl2, om, Ampl3, om, Ampl4);
xlim([700 1200]);
figure;
plot(om, Phase2, om, Phase3, om, Phase4);
xlim([700 1200]);

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
