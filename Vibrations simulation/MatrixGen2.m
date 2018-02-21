
% Modes:
% m = 1 : full conical beam
% m = 2 : truncated conical beam
% m = 3 : full conical beam with additional simple support
% m = 4 : truncated conical beam with additional simple support

function M = MatrixGen2(x0,c,L,beta,m)

if length(beta(:,1))==1
    beta = transpose(beta);
end

N = length(beta(:,1));

if m==1 % full conical beam
    
elseif m==2 % truncated conical beam
    
    gamma = (x0/L)*beta;
    
    M = zeros(4,4,N);
    
    s2beta = zeros(1,1,N);
    s2gamma = zeros(1,1,N);
    
    s2beta(1,1,:) = 2*sqrt(beta);
    s2gamma(1,1,:) = 2*sqrt(gamma);
    
    J2beta = besselj(2,s2beta); Y2beta = bessely(2,s2beta); I2beta = besseli(2,s2beta); K2beta = besselk(2,s2beta);
    J3beta = besselj(3,s2beta); Y3beta = bessely(3,s2beta); I3beta = besseli(3,s2beta); K3beta = besselk(3,s2beta);
    J4gamma = besselj(4,s2gamma); Y4gamma = bessely(4,s2gamma); I4gamma = besseli(4,s2gamma); K4gamma = besselk(4,s2gamma);
    J3gamma = besselj(3,s2gamma); Y3gamma = bessely(3,s2gamma); I3gamma = besseli(3,s2gamma); K3gamma = besselk(3,s2gamma);
    
    M(1,:,:) = [ J2beta, Y2beta, I2beta, K2beta ];
    M(2,:,:) = [ J3beta, Y3beta, -I3beta, K3beta ];
    M(3,:,:) = [ J4gamma, Y4gamma, I4gamma, K4gamma ];
    M(4,:,:) = [ J3gamma, Y3gamma, I3gamma, -K3gamma ];
    
elseif m==3 % full conical beam with additional simple support (6x6 matrix)

    eta = (c/L)*beta;
    
    M = zeros(6,6,N);
    
    s2eta = zeros(1,1,N);
    s2beta = zeros(1,1,N);
    Z = zeros(1,1,N);
    
    s2eta(1,1,:) = 2*sqrt(eta);
    s2beta(1,1,:) = 2*sqrt(beta);
    
    J2beta = besselj(2,s2beta); Y2beta = bessely(2,s2beta); I2beta = besseli(2,s2beta); K2beta = besselk(2,s2beta);
    J3beta = besselj(3,s2beta); Y3beta = bessely(3,s2beta); I3beta = besseli(3,s2beta); K3beta = besselk(3,s2beta);
    
    J2eta = besselj(2,s2eta); Y2eta = bessely(2,s2eta); I2eta = besseli(2,s2eta); K2eta = besselk(2,s2eta);
    J3eta = besselj(3,s2eta); Y3eta = bessely(3,s2eta); I3eta = besseli(3,s2eta); K3eta = besselk(3,s2eta);
    J4eta = besselj(4,s2eta); Y4eta = bessely(4,s2eta); I4eta = besseli(4,s2eta); K4eta = besselk(4,s2eta);
    
    M(1,:,:) = [ Z, Z, J2beta, Y2beta, I2beta, K2beta ];
    M(2,:,:) = [ Z, Z, -J3beta, -Y3beta, I3beta, -K3beta ];
    M(3,:,:) = [ J2eta, I2eta, Z, Z, Z, Z ];
    M(4,:,:) = [ Z, Z, J2eta, Y2eta, I2eta, K2eta ];
    M(5,:,:) = [ -J3eta, I3eta, J3eta, Y3eta, -I3eta, K3eta ];
    M(6,:,:) = [ J4eta, I4eta, -J4eta, -Y4eta, -I4eta, -K4eta ];

elseif m==4 % truncated conical beam with additional simple support (8x8 matrix)
    
    gamma = (x0/L)*beta;
    eta = (c/L)*beta;
    
    M = zeros(8,8,N);
    
    s2gamma = zeros(1,1,N);
    s2eta = zeros(1,1,N);
    s2beta = zeros(1,1,N);
    Z = zeros(1,1,N);
    
    s2gamma(1,1,:) = 2*sqrt(gamma);
    s2eta(1,1,:) = 2*sqrt(eta);
    s2beta(1,1,:) = 2*sqrt(beta);
    
    J2beta = besselj(2,s2beta); Y2beta = bessely(2,s2beta); I2beta = besseli(2,s2beta); K2beta = besselk(2,s2beta);
    J3beta = besselj(3,s2beta); Y3beta = bessely(3,s2beta); I3beta = besseli(3,s2beta); K3beta = besselk(3,s2beta);
    
    J2eta = besselj(2,s2eta); Y2eta = bessely(2,s2eta); I2eta = besseli(2,s2eta); K2eta = besselk(2,s2eta);
    J3eta = besselj(3,s2eta); Y3eta = bessely(3,s2eta); I3eta = besseli(3,s2eta); K3eta = besselk(3,s2eta);
    J4eta = besselj(4,s2eta); Y4eta = bessely(4,s2eta); I4eta = besseli(4,s2eta); K4eta = besselk(4,s2eta);
    
    J4gamma = besselj(4,s2gamma); Y4gamma = bessely(4,s2gamma); I4gamma = besseli(4,s2gamma); K4gamma = besselk(4,s2gamma);
    J5gamma = besselj(5,s2gamma); Y5gamma = bessely(5,s2gamma); I5gamma = besseli(5,s2gamma); K5gamma = besselk(5,s2gamma);

    M(1,:,:) = [ Z, Z, Z, Z, J2beta, Y2beta, I2beta, K2beta ];
    M(2,:,:) = [ Z, Z, Z, Z, -J3beta, -Y3beta, I3beta, -K3beta ];
    M(3,:,:) = [ J2eta, Y2eta, I2eta, K2eta, Z, Z, Z, Z ];
    M(4,:,:) = [ Z, Z, Z, Z, J2eta, Y2eta, I2eta, K2eta ];
    M(5,:,:) = [ -J3eta, -Y3eta, I3eta, -K3eta, J3eta, Y3eta, -I3eta, K3eta ];
    M(6,:,:) = [ J4eta, Y4eta, I4eta, K4eta, -J4eta, -Y4eta, -I4eta, -K4eta ];
    M(7,:,:) = [ J4gamma, Y4gamma, I4gamma, K4gamma, Z, Z, Z, Z ];
    M(8,:,:) = [ -J5gamma, -Y5gamma, I5gamma, -K5gamma, Z, Z, Z, Z ];
    
    
end




