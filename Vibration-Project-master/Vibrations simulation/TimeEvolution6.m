
% This script calculates the time evolution
% Should be used after ModeFinder.m

tic

PLOT_MODE = 0;

if ~PLOT_MODE

    t1 = 0; % start of contact
    t2 = 10e-3; % end of contact
    t3 = 30e-3; % end of calculation
% %     dt = 0.05e-3; 
% %     tsteps = (t2-t1)/dt + 1;
    tsteps = 100; %201; %500; %500; %1000; %400;

    xL = transpose(vec(x0,c,xpoints));
    xR = transpose(vec(c,L,xpoints));
    tilde_ys_right = 4/(pi*E*A^4) * ( (3*xR-c)./(6*xR.^2)+((3*L-2*c)*xR+3*L*(c-2*L))/(6*L^3) );
    tilde_ys_c = 4/(pi*E*A^4) * (1/(3*c)+((3*L-2*c)*c+3*L*(c-2*L))/(6*L^3));
    tilde_ys_prime_c = 4/(pi*E*A^4) * (-1/(6*c^2)+(3*L-2*c)/(6*L^3));
    tilde_ys_left = tilde_ys_prime_c*(xL-c) + tilde_ys_c;
    tilde_ys = [ tilde_ys_left; tilde_ys_right ];

    % Calculating the I_j integrals
    I_array = zeros(NZ,1);
    for j=1:1:NZ
        I_array(j,1) = trapz(x,tilde_ys.*x.^2.*X_array(:,j));
    end

    % --- Calculation of the time-dependent expansion coefficients ---
    % --- GAUSSIAN FORCE, smooth onset, sharp end ---
    tp_points = 40000; %100000; %30000;
    tau = 0.1e-3; % smooth onset duration
    W = t2; % force duration
    Fmax = -5e-6; %-20e-6; % peak force
    C = 0.5; % Gaussian cut-off
    zeta = -0.1; % if set zeta<0, then the standard constant alpha is used

    % --- Search for the Gaussian's center, a, and width, b ---
    amin = W/2;
    amax = (W+tau)/2;
    asteps = 100000;
    a_vec = transpose(vec(amin,amax,asteps));
    b = (W-a_vec)/(sqrt(log(1/C)));
    dF_tau = -2*Fmax./(b.^2*(1-C)).*(tau-a_vec).*exp(-((tau-a_vec)./b).^2);
    ddF_tau = 2*Fmax./(b.^2*(1-C)).*(2./b.^2.*(tau-a_vec).^2-1).*exp(-((tau-a_vec)./b).^2);
    F_R_tau = Fmax/(1-C)*(exp(-((tau-a_vec)./b).^2)-C);
    F_L_tau = tau/2*dF_tau - tau^2/12*ddF_tau;
    delta_F = F_R_tau - F_L_tau;
    Ma = FindMinima(a_vec,abs(delta_F));
    a = Ma(1,1);
    b = (W-a)/(sqrt(log(1/C)));
    % --- end of search ---

    dF_tau = -2*Fmax/(b*(1-C))*((tau-a)/b)*exp(-((tau-a)/b)^2);
    ddF_tau = 2*Fmax/(b^2*(1-C))*(2*((tau-a)/b)^2-1)*exp(-((tau-a)/b)^2);
    F_R_tau = Fmax/(1-C)*(exp(-((tau-a)/b)^2)-C);
    F_L_tau = tau/2*dF_tau-tau^2/12*ddF_tau;

    disp(['Force... Relative left condition: ',num2str(abs((F_R_tau-F_L_tau)/Fmax)),',   Rel. right condition : ',num2str(1/(1-C)*(exp(-((W-a)/b).^2)-C))])

    F_L = @(t)( tau/3*(3*dF_tau-tau*ddF_tau)*(t/tau).^3 + tau/4*(tau*ddF_tau-2*dF_tau)*(t/tau).^4 );
    dF_L = @(t)( (3*dF_tau-tau*ddF_tau)*(t/tau).^2 + (tau*ddF_tau-2*dF_tau)*(t/tau).^3 );
    ddF_L = @(t)( 2/tau*(3*dF_tau-tau*ddF_tau)*(t/tau) + (3/tau)*(tau*ddF_tau-2*dF_tau)*(t/tau).^2 );
    F_R = @(t)( Fmax/(1-C)*(exp(-((t-a)/b).^2)-C) );
    dF_R = @(t)( -2*Fmax/(b^2*(1-C))*(t-a).*exp(-((t-a)/b).^2) );
    ddF_R = @(t)( 2*Fmax/(b^2*(1-C))*(2/b^2*(t-a).^2-1).*exp(-((t-a)/b).^2) );

    phi_array = zeros(NZ,tsteps); % time-dependent expansion coefficients
    t_array = transpose(vec(t1,t2,tsteps)); % time points to calculate
    ys_array = zeros(2*xpoints,tsteps); % time-dependent steady-state profile
    F_array = zeros(tsteps,1); % force
    dF = zeros(tp_points,1);
    ddF = zeros(tp_points,1);
    final_tp_step = (t2-t1)/(tp_points-1);
    disp(['t-points per tau onset interval at t=t2: ',num2str(tau/final_tp_step)])
    disp(['t-points per sine period at t=t2: ',num2str(tp_points/(p_array(NZ)*t2/(2*pi)))])
    % % disp(['x-points per sine period at t=t2: ',num2str(tp_points/(p_array(end)*t2/(2*pi)))])
    phi_f = zeros(NZ,1); % phi at t=t2
    dphi_f = zeros(NZ,1); % time-derivative of phi at t=t2
    for qt=2:1:tsteps % Beginning from qt=2. At qt=1 phi is zero anyway
        t = t_array(qt,1);
        tp = transpose(vec(0,t,tp_points));
        if t<=tau % t<=tau
            F = F_L(t);
        elseif t>tau
            F = F_R(t);
        end
        F_array(qt,1) = F;
        ys_array(:,qt) = F * tilde_ys;
        for u=1:1:tp_points
            if tp(u,1)<=tau % t'<=tau
                dF(u,1) = dF_L(tp(u,1));
                ddF(u,1) = ddF_L(tp(u,1));
            elseif tp(u,1)>tau 
                dF(u,1) = dF_R(tp(u,1));
                ddF(u,1) = ddF_R(tp(u,1));
            end
        end
        for j=1:1:NZ
            p = p_array(j,1);
            if zeta>=0
                alpha_eff = 2*zeta*p;
            elseif zeta<0
                alpha_eff = alpha;
            end
            d = sqrt(p^2-(alpha_eff/2)^2);
            arg = -(ddF+alpha_eff*dF) .* exp(-alpha_eff/2*(t-tp)).*sin(d*(t-tp));
            phi_array(j,qt) = I_array(j,1)/d * trapz(tp,arg);
% %                 figure(50); clf
% %                 plot(tp,arg)
% %                 pause

            if qt==tsteps % collecting the values of phi and dphi at the last moment of contact (for initial conditions after the detachment)
                Q = -(ddF+alpha_eff*dF)*I_array(j,1);
                arg = Q .* exp(-alpha_eff/2*(t-tp)) .* ( cos(d*(t-tp)) - alpha_eff/(2*d)*sin(d*(t-tp)) );
                dphi_f(j,1) = trapz(tp,arg);
                phi_f(j,1) = phi_array(j,qt);
            end

        end
        
    end
    % ------- END OF THE COEFFICIENT CALCULATION (Also: force profile and steady-state whisker profile) -------


    % ------- FINDING THE VIBRATIONAL PART OF DEFLECTION DURING CONTACT -------
    ye_array = zeros(2*xpoints,tsteps);
    ye_f = zeros(2*xpoints,1); % this will be the vibrational profile at t=t2
    dye_f = zeros(2*xpoints,1); % this will be the time-derivative of vibrational profile at t=t2
    conv_sample = 20; % number of points along the whisker to sample for convergence
    conv_ind = round(transpose(vec(1,2*xpoints,conv_sample)));
    conv_time = round(tsteps*rand);
    conv_array = zeros(conv_sample,NZ);
    ddXL = zeros(tsteps,1);
    dddXL = zeros(tsteps,1);
    conv_array_ddXL = zeros(NZ,1);
    conv_array_dddXL = zeros(NZ,1);
    for qt=1:1:tsteps
        for j=1:1:NZ
            ye_array(:,qt) = ye_array(:,qt) + phi_array(j,qt)*X_array(:,j);
            ddXL(qt,1) = ddXL(qt,1) + phi_array(j,qt)*ddX_array(end,j);
            dddXL(qt,1) = dddXL(qt,1) + phi_array(j,qt)*dddX_array(end,j);
            if qt==conv_time
                conv_array(:,j) = ye_array(conv_ind,qt);
                conv_array_ddXL(j,1) = ddXL(qt,1);
                conv_array_dddXL(j,1) = dddXL(qt,1);
            end
            if qt==tsteps % time t2 reached
                ye_f(:,1) = ye_f(:,1) + phi_f(j,1)*X_array(:,j);
                dye_f(:,1) = dye_f(:,1) + dphi_f(j,1)*X_array(:,j);
            end
        end
    end

    % This is the FULL initial profile and its time-derivative 
    % used as input for time-evolution calculation beyond t=W=t2
    y_f = 0*tilde_ys + ye_f;
    dy_f = dF_R(t2)*tilde_ys + dye_f;

    J_array = zeros(NZ2,1);
    K_array = zeros(NZ2,1);
    for j=1:1:NZ2
        J_array(j,1) = trapz(x,X2_array(:,j).*x.^2.*y_f);
        K_array(j,1) = trapz(x,X2_array(:,j).*x.^2.*dy_f);
    end

    Ms = F_array*(L-c);
    Me = (pi/4)*E*A^4*L^4 * ddXL;

    Vs = F_array;
    Ve = (pi/4)*E*A^4*L^4 * ( 4/L*ddXL + dddXL );

    % Calculation of the mode activation hisgtogram 
    FT_array = zeros(NZ,1);
    for j=1:1:NZ
        FT_array(j,1) = mean(abs(phi_array(j,:)));
    end


    % Finding the coefficients after the whisker detachment
    varphi_array_after = zeros(NZ2,tsteps);
    t_array_after = transpose(vec(t2,t3,tsteps));
    t = t_array_after;
    for j=1:1:NZ2
        p = p2_array(j,1);
        if zeta>=0
            alpha_eff = 2*zeta*p;
        elseif zeta<0
            alpha_eff = alpha;
        end
        d = sqrt(p^2-(alpha_eff/2)^2);
        varphi_array_after(j,:) = exp(-alpha_eff/2*(t-t2)) .* ( J_array(j,1) * ( alpha_eff/(2*d)*sin(d*(t-t2)) + cos(d*(t-t2)) ) + K_array(j,1)/d*sin(d*(t-t2)) );
    end

    % ------- FINDING THE VIBRATIONAL PART AFTER THE DETACHMENT -------
    ye_array_after = zeros(2*xpoints,tsteps);
    ddXL_after = zeros(tsteps,1);
    dddXL_after = zeros(tsteps,1);
    conv_array2 = zeros(conv_sample,NZ2);
    conv_array2_ddXL = zeros(NZ2,1);
    conv_array2_dddXL = zeros(NZ2,1);
    for qt=1:1:tsteps
        for j=1:1:NZ2
            ye_array_after(:,qt) = ye_array_after(:,qt) + varphi_array_after(j,qt)*X2_array(:,j);
            ddXL_after(qt,1) = ddXL_after(qt,1) + varphi_array_after(j,qt)*ddX2_array(end,j);
            dddXL_after(qt,1) = dddXL_after(qt,1) + varphi_array_after(j,qt)*dddX2_array(end,j);
            if qt==conv_time
                conv_array2(:,j) = ye_array_after(conv_ind,qt);
                conv_array2_ddXL(j,1) = ddXL_after(qt,1);
                conv_array2_dddXL(j,1) = dddXL_after(qt,1);
            end
        end
    end

    Ms_after = zeros(tsteps,1);
    Me_after = (pi/4)*E*A^4*L^4 * ddXL_after;

    Vs_after = zeros(tsteps,1);
    Ve_after = (pi/4)*E*A^4*L^4 * ( 4/L*ddXL_after + dddXL_after );

    % % vel_sample = 3000; % 1 = whisker tip
    % % velocity = [ diff(transpose(ys_array(1,:)))./diff(t_array), diff(transpose(ye_array(1,:)))./diff(t_array) ];
    trace_points = 5;
    trace_ind = round(vec(1,2*xpoints,trace_points));
    y_trace = [ ys_array(trace_ind,:)+ye_array(trace_ind,:), ye_array_after(trace_ind,:) ];
    dy_trace_num = zeros(trace_points,2*tsteps-1);
    for j=1:1:trace_points
        dy_trace_num(j,:) = diff(y_trace(j,:))./transpose(diff([t_array(:,1); t_array_after(1:end,1)]));
    end
    % ------- END OF ALL CALCULATIONS -------

end % end of PLOT_MODE


figure(5); clf; 
subplot(1,2,1); box on; hold on
    title([num2str(conv_sample),' points sampled at t=',num2str(t_array(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e(x,t) during summation')
    for w=1:1:conv_sample
        plot(1:1:NZ,conv_array(w,:),'o-','markersize',2)
    end
    ScaleAxis(0.05)
subplot(1,2,2); box on; hold on
    title([num2str(conv_sample),' points sampled at t=',num2str(t_array_after(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e(x,t) during summation')
    for w=1:1:conv_sample
        plot(1:1:NZ2,conv_array2(w,:),'o-','markersize',2)
    end
    ScaleAxis(0.05)

figure(6); clf
subplot(1,2,1); box on; hold on
    title(['One point (x=L) sampled at t=',num2str(t_array(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e''''(x=L,t) during summation')
    plot(1:1:NZ,conv_array_ddXL,'o-','markersize',2)
subplot(1,2,2); box on; hold on
    title(['One point (x=L) sampled at t=',num2str(t_array_after(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e''''(x=L,t) during summation')
    plot(1:1:NZ2,conv_array2_ddXL,'o-','markersize',2)

figure(7); clf;
subplot(1,2,1); box on; hold on
    title(['One point (x=L) sampled at t=',num2str(t_array(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e''''''(x=L,t) during summation')
    plot(1:1:NZ,conv_array_dddXL,'o-','markersize',2)
subplot(1,2,2); box on; hold on
    title(['One point (x=L) sampled at t=',num2str(t_array_after(conv_time)*1e3),' ms (frame ',num2str(conv_time),'/',num2str(tsteps),')'])
    xlabel('eigenmode summation index, j')
    ylabel('value of y_e''''''(x=L,t) during summation')
    plot(1:1:NZ2,conv_array2_dddXL,'o-','markersize',2)




figure(8); clf; box on; hold on
xlabel('t (ms)')
ylabel('\phi_j(t)')
for j=1:1:NZ
    plot([t_array; t_array_after]*1e3,[phi_array(j,:), varphi_array_after(j,:)])
end

figure(8); clf; box on; hold on
plot([t_array(:,1); t_array_after(:,1)]*1e3,[Ms(:,1)*1e9; zeros(tsteps,1)])
title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   \tau = ',num2str(tau*1e3),' ms'])
xlabel('Time (ms)')
ylabel('Static bending moment at base (nN*m)')

figure(9); clf; box on; hold on
plot([t_array(:,1); t_array_after(:,1)]*1e3,[Me(:,1); Me_after(:,1)]*1e9)
title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   \tau = ',num2str(tau*1e3),' ms'])
xlabel('Time (ms)')
ylabel('Vibrational bending moment at base (nN*m)')
    disp(['Max vibrational bending moment at base = ',num2str(max(abs(Me(:,1)))*1e9),' nN*m'])
    disp(['Max static bending moment at base = ',num2str(max(abs(Ms(:,1)))*1e9),' nN*m'])
    disp(['Average vibrational bending moment at base = ',num2str(mean(abs(Me(:,1)))*1e9),' nN*m'])
    disp(['Average static bending moment at base = ',num2str(mean(abs(Ms(:,1)))*1e9),' nN*m'])

figure(10); clf; box on; hold on
plot([t_array(:,1); t_array_after(:,1)]*1e3,[(Ms(:,1)+Me(:,1)); Me_after(:,1)]*1e9)
title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   \tau = ',num2str(tau*1e3),' ms'])
xlabel('Time (ms)')
ylabel('Total bending moment at base (nN*m)')
set(gca,'fontsize',16)

figure(18); clf
subplot(1,3,1); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Vs(:,1); Vs_after(:,1)]*1e6)
    xlabel('Time (ms)')
    ylabel('Steady-state shear force at base (\muN)')
    title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   \tau = ',num2str(tau*1e3),' ms'])
    ScaleAxis(0.1)
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')
subplot(1,3,2); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Ve(:,1); Ve_after(:,1)]*1e6)
    xlabel('Time (ms)')
    ylabel('Vibrational shear force at base (\muN)')
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')
subplot(1,3,3); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Vs(:,1)+Ve(:,1); Vs_after(:,1)+Ve_after(:,1)]*1e6)
    xlabel('Time (ms)')
    ylabel('Total shear force at base (\muN)')
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')
disp(['Max vibrational shear force at base = ',num2str(max(abs(Ve(:,1)))*1e6),' uN'])
disp(['Max static shear force at base = ',num2str(max(abs(Vs(:,1)))*1e6),' uN'])
disp(['Average vibrational shear force at base = ',num2str(mean(abs(Ve(:,1)))*1e6),' uN'])
disp(['Average static shear force at base = ',num2str(mean(abs(Vs(:,1)))*1e6),' uN'])

figure(17); clf
subplot(1,3,1); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Ms(:,1); Ms_after(:,1)]*1e9)
    xlabel('Time (ms)')
    ylabel('Steady-state bending moment at base (nN*m)')
    title(['x_0/L = ',num2str(x0/L),'   c/L = ',num2str(c/L),'   \tau = ',num2str(tau*1e3),' ms'])
    ScaleAxis(0.1)
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')
subplot(1,3,2); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Me(:,1); Me_after(:,1)]*1e9)
    xlabel('Time (ms)')
    ylabel('Vibrational bending moment at base (nN*m)')
% %     ScaleAxis(0.1)
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')
subplot(1,3,3); box on; hold on
    plot([t_array(:,1); t_array_after(:,1)]*1e3,[ Ms(:,1)+Me(:,1); Ms_after(:,1)+Me_after(:,1)]*1e6)
    xlabel('Time (ms)')
    ylabel('Total bending moment at base (nN*m)')
% %     ScaleAxis(0.1)
    line([t2 t2]*1e3,ylim,'linestyle',':','color','k')

% % % Trajectory of a point close to the base
% % figure(11); clf; box on; hold on
% % plot(t_array(:,1)*1e3,ye_array(3990,:))

figure(12); clf; box on; hold on
xlabel('x (m)')
ylabel('y_s (m)')
for qt=1:1:tsteps
    plot(x,ys_array(:,qt))
end

% --- Plot of second derivative squared of y_s (small-angle approximation test) ---
figure(14); clf; box on; hold on
xlabel('x (m)')
ylabel('( dy_s/dx )^2 (dimensionless)')
for qt=1:1:tsteps
    dys = diff(ys_array(:,qt))./diff(x);
    plot(x(1:end-1),dys.^2)
end

figure(16); clf; box on
semilogy(1:1:NZ,FT_array(:,1),'o-')
xlabel('Eignemode #')
ylabel('Amplitude')
title(['\tau = ',num2str(tau*1e3),' ms. Amplitude averaged over ',num2str(t2*1e3),' ms'])



tL=transpose(vec(0,tau,100));
tR=transpose(vec(tau,W,1000));
fL = F_L(tL); %tau/3*(3*dF_tau-tau*ddF_tau) * (tL/tau).^3 + tau/4*(tau*ddF_tau-2*dF_tau) * (tL/tau).^4;
fR = F_R(tR); %tau/3*(3*dF_tau-tau*ddF_tau) * ((W-tR)/tau).^3 + tau/4*(tau*ddF_tau-2*dF_tau) * ((W-tR)/tau).^4;
dfL = dF_L(tL); %(3*dF_tau-tau*ddF_tau) * (tL/tau).^2 + (tau*ddF_tau-2*dF_tau) * (tL/tau).^3;
dfR = dF_R(tR); %- (3*dF_tau-tau*ddF_tau) * ((W-tR)/tau).^2 - (tau*ddF_tau-2*dF_tau) * ((W-tR)/tau).^3;
ddfL = ddF_L(tL); %2/tau*(3*dF_tau-tau*ddF_tau) * (tL/tau) + 3/tau*(tau*ddF_tau-2*dF_tau) * (tL/tau).^2;
ddfR = ddF_R(tR); %2/tau*(3*dF_tau-tau*ddF_tau) * ((W-tR)/tau) + 3/tau*(tau*ddF_tau-2*dF_tau) * ((W-tR)/tau).^2;

figure(20); clf
subplot(3,1,1); box on; hold on
    plot(tL*1e3,fL*1e6,'linewidth',2)
    plot(tR*1e3,fR*1e6,'linewidth',2)
    ScaleAxis(0.1)
    title(['Contact time = ',num2str(t2*1e3),' ms,   \tau = ',num2str(tau*1e3),' ms'])
    xlabel('Time (ms)')
    ylabel('Force, F(t) (\muN)')
subplot(3,1,2); box on; hold on
    plot(tL*1e3,dfL,'linewidth',2)
    plot(tR*1e3,dfR,'linewidth',2)
    ScaleAxis(0.1)
    xlabel('Time (ms)')
    ylabel('dF/dt')
subplot(3,1,3); box on; hold on
    plot(tL*1e3,ddfL,'linewidth',2)
    plot(tR*1e3,ddfR,'linewidth',2)
    ScaleAxis(0.1)
    xlabel('Time (ms)')
    ylabel('d^2F/dt^2')
figure(21); clf; box on; hold on
t_array_after = transpose(vec(t2,t3,tsteps));
plot([t_array; t_array_after]*1e3,[F_array*1e6; zeros(tsteps,1)],'linewidth',2)
line([1 1]*tau*1e3,[F_array(end) F_array(1)]*1e6,'linestyle','--','color','k')
title(['Contact time = ',num2str(t2*1e3),' ms,   \tau = ',num2str(tau*1e3),' ms'])
xlabel('time (ms)')
ylabel('Force (\muN)')
axis tight
ScaleAxis(0.1)

figure(21); clf; box on; hold on
plot([t_array(:,1); t_array_after(:,1)]*1e3,y_trace)
xlabel('time (ms)')
ylabel('y(x,t) traced at several points along the whisker')

figure(22); clf; box on; hold on
plot([t_array(:,1); t_array_after(1:end-1,1)]*1e3,dy_trace_num)
xlabel('time (ms)')
ylabel('time-derivative of y(x,t) traced at several points along the whisker')


toc


% --- THE ANIMATION SECTION ---
figure(15); clf; box on
clear f
f(tsteps) = struct('cdata',[],'colormap',[]);
% % Ymin = min(min(ye_array(:,:)));
% % Ymax = max(max(ye_array(:,:)));
Ymin = min(min(ys_array(:,:)));
Ymax = max(max(ys_array(:,:)));
Ymax = +max(abs([Ymin Ymax]));
Ymin = - Ymax;
pre_frames = 0;
% %     pause
% % for j=1:1:pre_frames
% %     plot(c,Ymax*(1-j/pre_frames),'o')
% %     line([x0 L],[0 0])    
% % % %     line([0 L],[0 0],'linestyle',':','color','k')
% %     line([c c],[Ymax Ymin],'linestyle',':','color','k')
% %     axis([0 L Ymin Ymax])
% %     xlabel('x (m)')
% %     ylabel('y_e(x,t) (m)')
% %     f(j) = getframe;
% % end
for qt=1:1:2*tsteps
    if qt<=tsteps
        plot(x,ye_array(:,qt)+ys_array(:,qt),'linewidth',1.5)
        line([0 L],[0 0],'linestyle',':','color','k')
        line([c c],[Ymin Ymax],'linestyle',':','color','k')
        axis([0 L Ymin Ymax ])
        xlabel('x (m)')
        ylabel('y_e(x,t) (m)')
        t = t_array(qt);
        text(L*0.80,Ymax*0.80,['t=',num2str(t*1e3,'%1.3f'),'ms'])
    else
        plot(x,ye_array_after(:,qt-tsteps),'linewidth',1.5)
        line([0 L],[0 0],'linestyle',':','color','k')
        line([c c],[Ymin Ymax],'linestyle',':','color','k')
        axis([0 L Ymin Ymax ])
        xlabel('x (m)')
        ylabel('y_e(x,t) (m)')
        t = t_array_after(qt-tsteps);
        text(L*0.80,Ymax*0.80,['t=',num2str(t*1e3,'%1.3f'),'ms'])
    end
    f(qt+pre_frames) = getframe;
    if qt==1
% %         pause
    end
end

movie(f,10)
% --- END OF ANIMATION SECTION ---


