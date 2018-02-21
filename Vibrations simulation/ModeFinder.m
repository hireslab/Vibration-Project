
% This code prepares the normalized vibrational eigenmodes
% The results are stored in:
% x (the x-axis)
% X_array, dX_array, ddX_array, dddX_array
% X2_array, dX2_array, ddX2_array, dddX2_array
% Z_array, p_array, nu_array
% Z2_array, p2_array, nu2_array

% % tic

Npoints = 10000; % number of points for the gross scan
beta_max =  30000; %20000; %3000; %20000;
red_factor = 100000;

FreqLimit = 1000e3; %1000e3; % true frequency (Hz)

zpoints = 1000; %5000; %100 % points for each section (to left and right from point of force)
xpoints = zpoints;

% --- Calculation parameters ---
% --- SI UNITS ---
L = 18e-3; % length from origin to whisker base, m
r0 = 37e-6; % radius at whisker base, m
x0 = 0.05*L; % TRUNCATED TIP POSITION, m
c = 0.30*L; % POLE POSITION, m
A = r0/L;
rho = 1.0e3; % whisker density, kg/m^3
E = 3.00e9; % Young's modulus, Pascal
alpha = 2 * 215; % damping constant, rad/sec

% ------- FINDING GROSS ZEROS -------

beta = transpose(vec(0,beta_max,Npoints));
M = MatrixGen2(x0,c,L,beta,4); % truncated conical w/supprot
M2 = MatrixGen2(x0,c,L,beta,2); % truncated conical w/o support

D = zeros(Npoints,1);
D2 = zeros(Npoints,1);
for j=1:1:Npoints
    D(j,1) = det(M(:,:,j)/red_factor);
    D2(j,1) = det(M2(:,:,j)); %/red_factor);
end

Z = FindMinima(beta,abs(D));
NZ = length(Z(:,1));
Z_array = zeros(NZ,6);
Z_array(:,2) = Z(:,1);
for j=1:1:NZ
    Z_gross = Z_array(j,2);
    ind = find(beta(:,1)==Z_gross);
    Z_left = beta(ind-1,1);
    Z_right = beta(ind+1,1);
    Z_array(j,3:4) = [ Z_left, Z_right ];
    Z_array(j,6) = abs(D(ind,1));
end

Z2 = FindMinima(beta,abs(D2));
NZ2 = length(Z2(:,1));
Z2_array = zeros(NZ2,6);
Z2_array(:,2) = Z2(:,1);
for j=1:1:NZ2
    Z2_gross = Z2_array(j,2);
    ind = find(beta(:,1)==Z2_gross);
    Z2_left = beta(ind-1,1);
    Z2_right = beta(ind+1,1);
    Z2_array(j,3:4) = [ Z2_left, Z2_right ];
    Z2_array(j,6) = abs(D2(ind,1));
end

% ------- FINDING PRECISE ZEROS -------

n_iter = 15; % number of iterations
ppi = 10; % points per iteration
% 1: Truncatred conical whisker with simple support added
for j=1:1:NZ
    Z_left = Z_array(j,3);
    Z0 = Z_array(j,2);
    Z_right = Z_array(j,4);
    beta_iter = transpose(vec(Z_left,Z_right,ppi));
    D_iter = zeros(ppi,1);
    for iter=1:1:n_iter
        M = MatrixGen2(x0,c,L,beta_iter,4);
        for w=1:1:ppi
            D_iter(w,1) = det(M(:,:,w)/red_factor);
        end
        for w=2:1:ppi-1
            D_left = abs(D_iter(w-1,1));
            D0 = abs(D_iter(w,1));
            D_right = abs(D_iter(w+1,1));
            if D_left>D0 && D0<D_right
                Z_left = beta_iter(w-1,1);
                Z0 = beta_iter(w,1);
                Z_right = beta_iter(w+1,1);
                beta_iter = transpose(vec(Z_left,Z_right,ppi));
                break
            end
        end
        disp(['1. Root: ',num2str(j),'   Iteration: ',num2str(iter),'   Value: ',num2str(Z0,'%1.15e')])
    end % end of iteration
    Z_array(j,1) = Z0;
    M = MatrixGen2(x0,c,L,Z0,4);
    Z_array(j,5) = abs(det(M/red_factor));
end % end of root

% 2: Truncatred conical whisker w/o middle support
for j=1:1:NZ2
    Z_left = Z2_array(j,3);
    Z0 = Z2_array(j,2);
    Z_right = Z2_array(j,4);
    beta_iter = transpose(vec(Z_left,Z_right,ppi));
    D_iter = zeros(ppi,1);
    for iter=1:1:n_iter
        M = MatrixGen2(x0,c,L,beta_iter,2);
        for w=1:1:ppi
            D_iter(w,1) = det(M(:,:,w)); %/red_factor);
        end
        for w=2:1:ppi-1
            D_left = abs(D_iter(w-1,1));
            D0 = abs(D_iter(w,1));
            D_right = abs(D_iter(w+1,1));
            if D_left>D0 && D0<D_right
                Z_left = beta_iter(w-1,1);
                Z0 = beta_iter(w,1);
                Z_right = beta_iter(w+1,1);
                beta_iter = transpose(vec(Z_left,Z_right,ppi));
                break
            end
        end
        disp(['2. Root: ',num2str(j),'   Iteration: ',num2str(iter),'   Value: ',num2str(Z0,'%1.15e')])
    end % end of iteration
    Z2_array(j,1) = Z0;
    M = MatrixGen2(x0,c,L,Z0,2);
    Z2_array(j,5) = abs(det(M)); %/red_factor));
end % end of root


figure(1); clf
semilogy(beta,abs(D))
hold on
for j=1:1:NZ
    line([Z_array(j,3) Z_array(j,1) Z_array(j,4)],[Z_array(j,6) Z_array(j,5) Z_array(j,6)],'color','r')
end
title('Truncated conical whisker w/simple support in the middle')
xlabel('dimensionless frequency, (\rho/E)^{1/2}(2L/A)p')
ylabel('det M')

figure(2); clf
semilogy(beta,abs(D2))
hold on
for j=1:1:NZ2
    line([Z2_array(j,3) Z2_array(j,1) Z2_array(j,4)],[Z2_array(j,6) Z2_array(j,5) Z2_array(j,6)],'color','r')
end
title('Truncated conical whisker w/o additional support')
xlabel('dimensionless frequency, (\rho/E)^{1/2}(2L/A)p')
ylabel('det M_2')


% ------- CONSTRUCTING THE EIGENMODES -------

p_array = Z_array(:,1) * sqrt(E/rho)*A/(2*L); % angular frequency (rad/sec)
nu_array = p_array/(2*pi); % real frequency (1/sec)

p2_array = Z2_array(:,1) * sqrt(E/rho)*A/(2*L); % angular frequency (rad/sec)
nu2_array = p2_array/(2*pi); % real frequency (1/sec)

for j=1:1:NZ
    if nu_array(j,1) > FreqLimit
        NZ = j-1;
        break
    end
end
for j=1:1:NZ2
    if nu2_array(j,1) > FreqLimit
        NZ2 = j-1;
        break
    end
end
if FreqLimit~=inf
    p_array(NZ+1:end)=[];
    nu_array(NZ+1:end)=[];
    p2_array(NZ2+1:end)=[];
    nu2_array(NZ2+1:end)=[];
end

% 1: Truncatred conical whisker with simple support added
chi_array = zeros(2*zpoints,NZ);
dchi_array = zeros(2*zpoints,NZ);
ddchi_array = zeros(2*zpoints,NZ);
dddchi_array = zeros(2*zpoints,NZ);
dX_array = zeros(2*zpoints,NZ);
ddX_array = zeros(2*zpoints,NZ);
dddX_array = zeros(2*zpoints,NZ);
for j=1:1:NZ
    
    beta_mode = Z_array(j,1);
    M = MatrixGen2(x0,c,L,beta_mode,4);
    M_small = M(1:7,1:7);
    b = - M(1:7,8);
    x = M_small\b;
    C = [ x; 1 ];
    
    q = beta_mode/L;
    
    z_left = transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints));
    z_right = transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints));
    
    chi_left = C(1)*besselj(2,2*sqrt(z_left))./z_left + C(2)*bessely(2,2*sqrt(z_left))./z_left + C(3)*besseli(2,2*sqrt(z_left))./z_left + C(4)*besselk(2,2*sqrt(z_left))./z_left;
    chi_right = C(5)*besselj(2,2*sqrt(z_right))./z_right + C(6)*bessely(2,2*sqrt(z_right))./z_right + C(7)*besseli(2,2*sqrt(z_right))./z_right + C(8)*besselk(2,2*sqrt(z_right))./z_right;
    CN = 1/sqrt( (1./q^3) * trapz([z_left; z_right],[z_left; z_right].^2.*[chi_left; chi_right].^2));
    chi_left = CN*chi_left;
    chi_right = CN*chi_right;
    C = CN*C;

    dchi_left = -C(1)*besselj(3,2*sqrt(z_left))./z_left.^(3/2) - C(2)*bessely(3,2*sqrt(z_left))./z_left.^(3/2) + C(3)*besseli(3,2*sqrt(z_left))./z_left.^(3/2) - C(4)*besselk(3,2*sqrt(z_left))./z_left.^(3/2);
    dchi_right = -C(5)*besselj(3,2*sqrt(z_right))./z_right.^(3/2) - C(6)*bessely(3,2*sqrt(z_right))./z_right.^(3/2) + C(7)*besseli(3,2*sqrt(z_right))./z_right.^(3/2) - C(8)*besselk(3,2*sqrt(z_right))./z_right.^(3/2);
    ddchi_left = C(1)*besselj(4,2*sqrt(z_left))./z_left.^(4/2) + C(2)*bessely(4,2*sqrt(z_left))./z_left.^(4/2) + C(3)*besseli(4,2*sqrt(z_left))./z_left.^(4/2) + C(4)*besselk(4,2*sqrt(z_left))./z_left.^(4/2);
    ddchi_right = C(5)*besselj(4,2*sqrt(z_right))./z_right.^(4/2) + C(6)*bessely(4,2*sqrt(z_right))./z_right.^(4/2) + C(7)*besseli(4,2*sqrt(z_right))./z_right.^(4/2) + C(8)*besselk(4,2*sqrt(z_right))./z_right.^(4/2);
    dddchi_left = -C(1)*besselj(5,2*sqrt(z_left))./z_left.^(5/2) - C(2)*bessely(5,2*sqrt(z_left))./z_left.^(5/2) + C(3)*besseli(5,2*sqrt(z_left))./z_left.^(5/2) - C(4)*besselk(5,2*sqrt(z_left))./z_left.^(5/2);
    dddchi_right = -C(5)*besselj(5,2*sqrt(z_right))./z_right.^(5/2) - C(6)*bessely(5,2*sqrt(z_right))./z_right.^(5/2) + C(7)*besseli(5,2*sqrt(z_right))./z_right.^(5/2) - C(8)*besselk(5,2*sqrt(z_right))./z_right.^(5/2);
    
    chi_array(:,j) = [ chi_left; chi_right ];
    dchi_array(:,j) = [ dchi_left; dchi_right ];
    ddchi_array(:,j) = [ ddchi_left; ddchi_right ];
    dddchi_array(:,j) = [ dddchi_left; dddchi_right ];

    dX_array(:,j) = q * dchi_array(:,j);
    ddX_array(:,j) = q^2 * ddchi_array(:,j);
    dddX_array(:,j) = q^3 * dddchi_array(:,j);
    
end

% % p_array = Z_array(:,1) * sqrt(E/rho)*A/(2*L); % angular frequency (rad/sec)
% % nu_array = p_array/(2*pi); % real frequency (1/sec)

% 2: Truncatred conical whisker w/o middle support
chi2_array = zeros(2*zpoints,NZ2);
dchi2_array = zeros(2*zpoints,NZ2);
ddchi2_array = zeros(2*zpoints,NZ2);
dddchi2_array = zeros(2*zpoints,NZ2);
dX2_array = zeros(2*zpoints,NZ2);
ddX2_array = zeros(2*zpoints,NZ2);
dddX2_array = zeros(2*zpoints,NZ2);
for j=1:1:NZ2
    
    beta_mode = Z2_array(j,1);
    M = MatrixGen2(x0,c,L,beta_mode,2);
    M_small = M(1:3,1:3);
    b = - M(1:3,4);
    x = M_small\b;
    C = [ x; 1 ];

    q = beta_mode/L;
    
    z_left = transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints));
    z_right = transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints));
    
    z = [ z_left; z_right ];
    
    chi = C(1)*besselj(2,2*sqrt(z))./z + C(2)*bessely(2,2*sqrt(z))./z + C(3)*besseli(2,2*sqrt(z))./z + C(4)*besselk(2,2*sqrt(z))./z;
    CN = 1/sqrt( (1./q^3) * trapz(z,z.^2.*chi.^2));
    chi = CN*chi;
    C = CN*C;

    dchi = -C(1)*besselj(3,2*sqrt(z))./z.^(3/2) - C(2)*bessely(3,2*sqrt(z))./z.^(3/2) + C(3)*besseli(3,2*sqrt(z))./z.^(3/2) - C(4)*besselk(3,2*sqrt(z))./z.^(3/2);
    ddchi = C(1)*besselj(4,2*sqrt(z))./z.^(4/2) + C(2)*bessely(4,2*sqrt(z))./z.^(4/2) + C(3)*besseli(4,2*sqrt(z))./z.^(4/2) + C(4)*besselk(4,2*sqrt(z))./z.^(4/2);
    dddchi = -C(1)*besselj(5,2*sqrt(z))./z.^(5/2) - C(2)*bessely(5,2*sqrt(z))./z.^(5/2) + C(3)*besseli(5,2*sqrt(z))./z.^(5/2) - C(4)*besselk(5,2*sqrt(z))./z.^(5/2);
    
    chi2_array(:,j) = chi;
    dchi2_array(:,j) = dchi;
    ddchi2_array(:,j) = ddchi;
    dddchi2_array(:,j) = dddchi;

    dX2_array(:,j) = q * dchi2_array(:,j);
    ddX2_array(:,j) = q^2 * ddchi2_array(:,j);
    dddX2_array(:,j) = q^3 * dddchi2_array(:,j);
    
end

% % p2_array = Z2_array(:,1) * sqrt(E/rho)*A/(2*L); % angular frequency (rad/sec)
% % nu2_array = p2_array/(2*pi); % real frequency (1/sec)

x = [ transpose(vec(x0,c,zpoints)); transpose(vec(c,L,zpoints)) ];

X_array = chi_array;
X2_array = chi2_array;

    xL = transpose(vec(x0,c,xpoints));
    xR = transpose(vec(c,L,xpoints));
    tilde_ys_right = 4/(pi*E*A^4) * ( (3*xR-c)./(6*xR.^2)+((3*L-2*c)*xR+3*L*(c-2*L))/(6*L^3) );
    tilde_ys_c = 4/(pi*E*A^4) * (1/(3*c)+((3*L-2*c)*c+3*L*(c-2*L))/(6*L^3));
    tilde_ys_prime_c = 4/(pi*E*A^4) * (-1/(6*c^2)+(3*L-2*c)/(6*L^3));
    tilde_ys_left = tilde_ys_prime_c*(xL-c) + tilde_ys_c;
    tilde_ys = [ tilde_ys_left; tilde_ys_right ];

    I_array = zeros(NZ,1);
    omega_array = p_array;
    for W=1:1:NZ
        I_array(W,1) = trapz(x,tilde_ys.*x.^2.*X_array(:,W));
    end
    Amp_array = abs(I_array)./omega_array;
    MAmp_array = (pi/4)*E*A^4*L^4*abs(transpose(ddX_array(end,:)).*I_array)./omega_array;

FM = FindMinima(x,abs(X_array(:,end)));
dmin = FM(2,1) - FM(1,1); % minimum distance between two adjacent nodes of the fastest eigenmode
FM2 = FindMinima(x,abs(X2_array(:,end)));
dmin2 = FM2(2,1) - FM2(1,1); % minimum distance between two adjacent nodes of the fastest eigenmode

% --- PLOTTING MODES W/SUPPORT ---
figure(3); clf
subplot(2,2,1)
box on; hold on
title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   NZ = ',num2str(NZ)])
xlabel('x (meters)')
ylabel('X_j(x) (meters)')
norm_dev = 1;
for j=1:1:NZ
    plot(x,X_array(:,j))
    norm_check = trapz(x,x.^2.*X_array(:,j).^2);
    if abs(norm_dev-1)<abs(norm_check-1)
        norm_dev=norm_check;
    end
    disp(num2str(norm_check,'%1.16f'))
end
% % ScaleAxis(0.05)
line(xlim,[0 0],'linestyle',':','color','k')
line([c c],ylim,'linestyle',':','color','k')

subplot(2,2,2); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^\prime_j(z)')
for j=1:1:NZ
    beta_mode = Z_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,dchi_array(:,j))
end
% % ScaleAxis(0.05)
line(xlim,[0 0],'linestyle',':','color','k')
line([c c]/L,ylim,'linestyle',':','color','k')

subplot(2,2,3); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^{\prime\prime}_j(z)')
for j=1:1:NZ
    beta_mode = Z_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,ddchi_array(:,j))
end
line(xlim,[0 0],'linestyle',':','color','k')
line([c c]/L,ylim,'linestyle',':','color','k')

subplot(2,2,4); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^{\prime\prime\prime}_j(z)')
for j=1:1:NZ
    beta_mode = Z_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,dddchi_array(:,j))
end
line(xlim,[0 0],'linestyle',':','color','k')
line([c c]/L,ylim,'linestyle',':','color','k')


% --- PLOTTING MODES W/O SUPPORT ---
figure(4); clf
subplot(2,2,1); box on; hold on
title(['x_0/L = ',num2str(x0/L),'   NZ2 = ',num2str(NZ2)])
xlabel('x (meters)')
ylabel('X_j(x) (meters)')
norm_dev2 = 1;
for j=1:1:NZ2
    plot(x,X2_array(:,j))
    norm_check = trapz(x,x.^2.*X2_array(:,j).^2);
    if abs(norm_dev2-1)<abs(norm_check-1)
        norm_dev2=norm_check;
    end
    disp(num2str(norm_check,'%1.16f'))
end
line(xlim,[0 0],'linestyle',':','color','k')
% % line([c c],ylim,'linestyle',':','color','k')

subplot(2,2,2); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^\prime_j(z)')
for j=1:1:NZ2
    beta_mode = Z2_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,dchi2_array(:,j))
end
line(xlim,[0 0],'linestyle',':','color','k')
% % line([c c]/L,ylim,'linestyle',':','color','k')

subplot(2,2,3); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^{\prime\prime}_j(z)')
for j=1:1:NZ2
    beta_mode = Z2_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,ddchi2_array(:,j))
end
line(xlim,[0 0],'linestyle',':','color','k')
% % line([c c]/L,ylim,'linestyle',':','color','k')

subplot(2,2,4); box on; hold on
xlabel('z/\beta_j')
ylabel('\chi^{\prime\prime\prime}_j(z)')
for j=1:1:NZ2
    beta_mode = Z2_array(j,1);
    z = [ transpose(vec(x0/L*beta_mode,c/L*beta_mode,zpoints)); transpose(vec(c/L*beta_mode,L/L*beta_mode,zpoints)) ];
    plot(z/beta_mode,dddchi2_array(:,j))
end
line(xlim,[0 0],'linestyle',':','color','k')
% % line([c c]/L,ylim,'linestyle',':','color','k')





% ------- 1. BCs check w/supprt -------
disp(' ')
disp('1: Vanishing of X at base')
max_dev1 = 0;
for j=1:1:NZ
    Xmax = max(abs(X_array(:,j)));
    Xbase = X_array(end,j);
    rel_dev = Xbase/Xmax;
    if max_dev1<abs(rel_dev)
        max_dev1 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('2: Vanishing of X'' at base')
max_dev2 = 0;
for j=1:1:NZ
    dXmax = max(abs(dX_array(:,j)));
    dXbase = dX_array(end,j);
    rel_dev = dXbase/dXmax;
    if max_dev2<abs(rel_dev)
        max_dev2 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('3,4: Vanishing of X at pole (on left and right)')
max_dev3=0;
max_dev4=0;
for j=1:1:NZ
    Xmax = max(abs(X_array(:,j)));
    XL = X_array(xpoints,j);
    XR = X_array(xpoints+1,j);
    rel_dev3 = XL/Xmax;
    rel_dev4 = XR/Xmax;
    if max_dev3<abs(rel_dev3)
        max_dev3 = abs(rel_dev3);
    end
    if max_dev4<abs(rel_dev4)
        max_dev4 = abs(rel_dev4);
    end
    disp([num2str(rel_dev3,'%1.16e'),'   ',num2str(rel_dev4,'%1.16e')])
end

disp(' ')
disp('5: Continuity of X'' at pole')
max_dev5=0;
for j=1:1:NZ
    dXpoleL = dX_array(xpoints,j);
    dXpoleR = dX_array(xpoints+1,j);
    rel_dev = (dXpoleL-dXpoleR)/dXpoleL;
    if max_dev5<abs(rel_dev)
        max_dev5 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('6: Continuity of X'''' at pole')
max_dev6 = 0;
for j=1:1:NZ
    ddXpoleL = ddX_array(xpoints,j);
    ddXpoleR = ddX_array(xpoints+1,j);
    rel_dev = (ddXpoleL-ddXpoleR)/ddXpoleL;
    if max_dev6<abs(rel_dev)
        max_dev6 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('7: Vanishing of X'''' at tip')
max_dev7=0;
for j=1:1:NZ
    ddXmax = max(abs(ddX_array(:,j)));
    ddXtip = ddX_array(1,j);
    rel_dev = ddXtip/ddXmax;
    if max_dev7<abs(rel_dev)
        max_dev7 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('8: Vanishing of X'''''' at tip')
max_dev8=0;
for j=1:1:NZ
    dddXmax = max(abs(dddX_array(:,j)));
    dddXtip = dddX_array(1,j);
    rel_dev = dddXtip/dddXmax;
    if max_dev8<abs(rel_dev)
        max_dev8 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end


% ------- 2. BCs check w/o supprt -------
disp(' ')
disp('1: Vanishing of X at base')
max_dev9 = 0;
for j=1:1:NZ2
    Xmax = max(abs(X2_array(:,j)));
    Xbase = X2_array(end,j);
    rel_dev = Xbase/Xmax;
    if max_dev9<abs(rel_dev)
        max_dev9 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('2: Vanishing of X'' at base')
max_dev10 = 0;
for j=1:1:NZ2
    dXmax = max(abs(dX2_array(:,j)));
    dXbase = dX2_array(end,j);
    rel_dev = dXbase/dXmax;
    if max_dev10<abs(rel_dev)
        max_dev10 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('3: Vanishing of X'''' at tip')
max_dev11=0;
for j=1:1:NZ2
    ddXmax = max(abs(ddX2_array(:,j)));
    ddXtip = ddX2_array(1,j);
    rel_dev = ddXtip/ddXmax;
    if max_dev11<abs(rel_dev)
        max_dev11 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end

disp(' ')
disp('4: Vanishing of X'''''' at tip')
max_dev12=0;
for j=1:1:NZ2
    dddXmax = max(abs(dddX2_array(:,j)));
    dddXtip = dddX2_array(1,j);
    rel_dev = dddXtip/dddXmax;
    if max_dev12<abs(rel_dev)
        max_dev12 = abs(rel_dev);
    end
    disp(num2str(rel_dev,'%1.16e'))
end



disp('EIGENMODES W/SUPPORT. Largest BC relative deviations:')
disp(['BC 1: ',num2str(max_dev1,'%1.16e'),'   BC 2: ',num2str(max_dev2,'%1.16e'),'   BC 3: ',num2str(max_dev3,'%1.16e'),'   BC 4: ',num2str(max_dev4,'%1.16e')])
disp(['BC 5: ',num2str(max_dev5,'%1.16e'),'   BC 6: ',num2str(max_dev6,'%1.16e'),'   BC 7: ',num2str(max_dev7,'%1.16e'),'   BC 8: ',num2str(max_dev8,'%1.16e')])
disp(['Most deviating normalization: ',num2str(norm_dev,'%1.16f')])
disp(['Eigenmodes found: ',num2str(NZ)])

disp('EIGENMODES W/O SUPPORT. Largest BC relative deviations:')
disp(['BC 1: ',num2str(max_dev9,'%1.16e'),'   BC 2: ',num2str(max_dev10,'%1.16e'),'   BC 3: ',num2str(max_dev11,'%1.16e'),'   BC 4: ',num2str(max_dev12,'%1.16e')])
disp(['Most deviating normalization: ',num2str(norm_dev2,'%1.16f')])
disp(['Eigenmodes found: ',num2str(NZ2)])

disp(['x-points per fastest half-period... W/SUPPORT: ',num2str(dmin/(L-x0)*2*xpoints),'   W/O SUPPORT: ',num2str(dmin2/(L-x0)*2*xpoints)])


% % toc

% % % --- PLOTS FOR PAPER ---
% % figure(5); clf; box on; hold on
% % %title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   NZ = ',num2str(NZ)])
% % xlabel('Normalized position along whisker, x/L')
% % ylabel('Amplitude, X_j(x) (arb. units)')
% % norm_dev = 1;
% % maxM = max(max(abs(X_array(:,:))));
% % for j=1:1:4 %NZ
% %     plot(x/L,X_array(:,j)/(maxM/2),'linewidth',1.5)
% %     norm_check = trapz(x,x.^2.*X_array(:,j).^2);
% %     if abs(norm_dev-1)<abs(norm_check-1)
% %         norm_dev=norm_check;
% %     end
% %     disp(num2str(norm_check,'%1.16f'))
% % end
% % legend('Mode 1','2','3','4')
% % ScaleAxis(0.05)
% % line(xlim,[0 0],'linestyle',':','color','k','linewidth',1.0)
% % line([c c]/L,ylim,'linestyle',':','color','k','linewidth',1.0)
% % set(gca,'fontsize',12)
% % 
% % figure(6); clf; box on; hold on
% % %title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   NZ = ',num2str(NZ)])
% % xlabel('Normalized position along whisker, x/L')
% % ylabel('Amplitude, Z_j(x) (arb. units)')
% % norm_dev = 1;
% % maxM = max(max(abs(X2_array(:,:))));
% % for j=1:1:4 %NZ
% %     plot(x/L,X2_array(:,j)/(maxM/2),'linewidth',1.5)
% %     norm_check = trapz(x,x.^2.*X2_array(:,j).^2);
% %     if abs(norm_dev-1)<abs(norm_check-1)
% %         norm_dev=norm_check;
% %     end
% %     disp(num2str(norm_check,'%1.16f'))
% % end
% % legend('Mode 1','2','3','4')
% % ScaleAxis(0.05)
% % line(xlim,[0 0],'linestyle',':','color','k','linewidth',1.0)
% % % % line([c c]/L,ylim,'linestyle',':','color','k','linewidth',1.0)
% % set(gca,'fontsize',12)

