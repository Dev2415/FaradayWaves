
%% Read parameters
Re = 1;
We = 0.1;
f0 = 1;
rho1 = 1;
rho2 = 1e-2;
eta1 = 1/Re;
eta2 = 1e-2 * eta1;
h1 = 1;
h2 = 1;
sigma = 1/We;
kmin = 100e4;
kmax = 100e5;
f =  f0;
mm = 4;
ll = 5;
chi = 60*pi/180;
g = 9.8066;
amax = 15.0*g;
CASE = 3;
omega = 2*pi*f;
warning('off');
phi = 0;
ephiP = exp(complex(0, phi));
ephiM = exp(complex(0, -phi));
%% Subharmonic
alpha = 0.5*omega;
subharmonicData = [];
N = 20;
mu = 0;
for zz = kmin:1000:kmax
    k = zz;
    A = zeros(2*(N+1),2*(N+1));
    for n = 0:N
        floquetExponent = complex(mu, alpha+n*omega);
        q1 = sqrt(k^2+rho1/eta1*floquetExponent);
        q2 = sqrt(k^2+rho2/eta2*floquetExponent);
        b = zeros(8,1);
        b(1) = floquetExponent;
        if CASE == 2
        boundaries = [1,1,1,1,0,0,0,0;
            1,1,1,1,-1,-1,-1,-1;
            k, -k, q1, -q1, -k, k, -q2, q2;
            2*eta1*k^2, 2*eta1*k^2, eta1*(q1^2+k^2), eta1*(q1^2+k^2),-2*eta2*k^2, -2*eta2*k^2, -eta2*(q2^2+k^2), -eta2*(q2^2+k^2);
            exp(-k*h1), exp(k*h1), exp(-q1*h1), exp(q1*h1), 0, 0, 0, 0;
            k*exp(-k*h1), -k*exp(k*h1), q1*exp(-q1*h1), -q1*exp(q1*h1), 0, 0, 0, 0;
            0, 0, 0, 0, exp(k*h2), exp(-k*h2), exp(q2*h2), exp(-q2*h2);
            0, 0, 0, 0, k*exp(k*h2), -k*exp(-k*h2), q2*exp(q2*h2), -q2*exp(-q2*h2)];
        elseif CASE == 1
            boundaries = [1,1,1,1,0,0,0,0;
            1,1,1,1,-1,-1,-1,-1;
            k, -k, q1, -q1, -k, k, -q2, q2;
            2*eta1*k^2, 2*eta1*k^2, eta1*(q1^2+k^2), eta1*(q1^2+k^2),-2*eta2*k^2, -2*eta2*k^2, -eta2*(q2^2+k^2), -eta2*(q2^2+k^2);
            0, 1, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 1, 0];

        else
            boundaries = [1,1,1,1,0,0,0,0;
            1,1,1,1,-1,-1,-1,-1;
            k, -k, q1, -q1, -k, k, -q2, q2;
            2*eta1*k^2, 2*eta1*k^2, eta1*(q1^2+k^2), eta1*(q1^2+k^2),-2*eta2*k^2, -2*eta2*k^2, -eta2*(q2^2+k^2), -eta2*(q2^2+k^2);
            exp(-k*h1), exp(k*h1), exp(-q1*h1), exp(q1*h1), 0, 0, 0, 0;
            k*exp(-k*h1), -k*exp(k*h1), q1*exp(-q1*h1), -q1*exp(q1*h1), 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 1, 0];
        end

        coefficients = boundaries\b;

        dw1 = coefficients(1)*k-coefficients(2)*k+coefficients(3)*q1-coefficients(4)*q1;
        dwww1 = coefficients(1)*k^3-coefficients(2)*k^3+coefficients(3)*q1^3-coefficients(4)*q1^3;
        dw2 = coefficients(5)*k-coefficients(6)*k+coefficients(7)*q1-coefficients(8)*q1;
        dwww2 = coefficients(5)*k^3-coefficients(6)*k^3+coefficients(7)*q2^3-coefficients(8)*q2^3;

        An = (rho2*floquetExponent+3*eta2*k^2)*dw2-eta2*dwww2-...
            ((rho1*floquetExponent+3*eta1*k^2)*dw1-eta1*dwww1)+ ...
            ((rho2-rho1)*g-sigma*k^2)*k^2;

        A(2*n+1, 2*n+1) = real(An);
        A(2*n+1, 2*n+2) = -imag(An);
        A(2*n+2, 2*n+1) = imag(An);
        A(2*n+2,2*n+2) = real(An);
    end


%     B = zeros(2*(N+1),2*(N+1));
%     B(1,1) = 1; B(2,2) = -1;
%     for II = 1 : 2*(N+1)-2
%         B(II, II+2) = 1;
%         B(II+2, II) = 1;
%     end
    
    % B general
    B2 = zeros(2*(N+1), 2*(N+1));
    for II = 0 : N
    
        if II+mm <= N
            B2(idr(II), idr(II+mm)) = B2(idr(II), idr(II+mm))+ cos(chi);
            B2(idi(II), idi(II+mm)) = B2(idi(II), idi(II+mm))+ cos(chi);
        end

        if II-mm >= 0
            B2(idr(II), idr(II-mm)) = B2(idr(II),idr(II-mm)) + cos(chi);
            B2(idi(II), idi(II-mm)) = B2(idi(II),idi(II-mm)) + cos(chi);
        else
            B2(idr(II), idr(-(II-mm)-1)) = B2(idr(II), idr(-(II-mm)-1)) + cos(chi);
            B2(idi(II), idi(-(II-mm)-1)) = B2(idi(II), idi(-(II-mm)-1)) - cos(chi); 
        end

        if II+ll <= N
            B2(idr(II), idr(II+ll)) = B2(idr(II), idr(II+ll))+ sin(chi)*ephiM;
            B2(idi(II), idi(II+ll)) = B2(idi(II), idi(II+ll))+ sin(chi)*ephiM;
        end

        if II-ll >= 0
            B2(idr(II), idr(II-ll)) = B2(idr(II),idr(II-ll)) + sin(chi)*ephiP;
            B2(idi(II), idi(II-ll)) = B2(idi(II),idi(II-ll)) + sin(chi)*ephiP;
        else
            B2(idr(II), idr(-(II-ll)-1)) = B2(idr(II), idr(-(II-ll)-1)) + sin(chi)*ephiP;
            B2(idi(II), idi(-(II-ll)-1)) = B2(idi(II), idi(-(II-ll)-1)) - sin(chi)*ephiP; 
        end
    end

    B2 = B2*0.5*(rho2-rho1)*k^2;
    AB = A\B2;
    EV = eig(AB);
    amplitude = 1./EV;
    amplitude = amplitude(amplitude>0);
    amplitude = amplitude(amplitude<amax);
    amplitude = amplitude(imag(amplitude)==0);
    for i = 1 : length(amplitude)
        subharmonicData = [subharmonicData; amplitude(i)/g, k/1000];
    end
end

%% Harmonic
alpha = 0;
harmonicData = [];
for zz = kmin:1000:kmax
    k = zz;
    A = zeros(2*(N+1),2*(N+1));
    for n = 0:N
        floquetExponent = complex(mu, alpha+n*omega);
        q1 = sqrt(k^2+rho1/eta1*floquetExponent);
        q2 = sqrt(k^2+rho2/eta2*floquetExponent);
        b = zeros(8,1);
        b(1) = floquetExponent;
        if CASE == 2
        boundaries = [1,1,1,1,0,0,0,0;
            1,1,1,1,-1,-1,-1,-1;
            k, -k, q1, -q1, -k, k, -q2, q2;
            2*eta1*k^2, 2*eta1*k^2, eta1*(q1^2+k^2), eta1*(q1^2+k^2),-2*eta2*k^2, -2*eta2*k^2, -eta2*(q2^2+k^2), -eta2*(q2^2+k^2);
            exp(-k*h1), exp(k*h1), exp(-q1*h1), exp(q1*h1), 0, 0, 0, 0;
            k*exp(-k*h1), -k*exp(k*h1), q1*exp(-q1*h1), -q1*exp(q1*h1), 0, 0, 0, 0;
            0, 0, 0, 0, exp(k*h2), exp(-k*h2), exp(q2*h2), exp(-q2*h2);
            0, 0, 0, 0, k*exp(k*h2), -k*exp(-k*h2), q2*exp(q2*h2), -q2*exp(-q2*h2)];
        else
        boundaries = [1,1,1,1,0,0,0,0;
            1,1,1,1,-1,-1,-1,-1;
            k, -k, q1, -q1, -k, k, -q2, q2;
            2*eta1*k^2, 2*eta1*k^2, eta1*(q1^2+k^2), eta1*(q1^2+k^2),-2*eta2*k^2, -2*eta2*k^2, -eta2*(q2^2+k^2), -eta2*(q2^2+k^2);
            0, 1, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0, 0, 0;
            0, 0, 0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 1, 0];
        end
        if n==0
            coefficients = zeros(8,1);
        else
            coefficients = boundaries\b;
        end
        
        dw1 = coefficients(1)*k-coefficients(2)*k+coefficients(3)*q1-coefficients(4)*q1;
        dwww1 = coefficients(1)*k^3-coefficients(2)*k^3+coefficients(3)*q1^3-coefficients(4)*q1^3;
        dw2 = coefficients(5)*k-coefficients(6)*k+coefficients(7)*q1-coefficients(8)*q1;
        dwww2 = coefficients(5)*k^3-coefficients(6)*k^3+coefficients(7)*q2^3-coefficients(8)*q2^3;

        An = (rho2*floquetExponent+3*eta2*k^2)*dw2-eta2*dwww2-...
            ((rho1*floquetExponent+3*eta1*k^2)*dw1-eta1*dwww1)+ ...
            ((rho2-rho1)*g-sigma*k^2)*k^2;

        A(2*n+1, 2*n+1) = real(An);
        A(2*n+1, 2*n+2) = -imag(An);
        A(2*n+2, 2*n+1) = imag(An);
        A(2*n+2,2*n+2) = real(An);
    end
    B = zeros(2*(N+1),2*(N+1));
    
    for II = 1 : 2*(N+1)-2
        B(II, II+2) = 1;
        B(II+2, II) = 1;
    end
    B(1,3) = 2; B(2,4) = 0;

    % B general
    B4 = zeros(2*(N+1),2*(N+1));

    for II = 0 : N

     if II+mm <= N
         B4(idr(II), idr(II+mm)) = B4(idr(II), idr(II+mm))+ cos(chi);
         B4(idi(II), idi(II+mm)) = B4(idi(II), idi(II+mm))+ cos(chi);
     end

     if II-mm >= 0
         B4(idr(II), idr(II-mm)) = B4(idr(II),idr(II-mm)) + cos(chi);
         B4(idi(II), idi(II-mm)) = B4(idi(II),idi(II-mm)) + cos(chi);
     else
         B4(idr(II), idr(-(II-mm))) = B4(idr(II), idr(-(II-mm))) + cos(chi);
         B4(idi(II), idi(-(II-mm))) = B4(idi(II), idi(-(II-mm))) - cos(chi);
     end
    
     if II+ll <= N
         B4(idr(II), idr(II+ll)) = B4(idr(II), idr(II+ll))+ sin(chi)*ephiM;
         B4(idi(II), idi(II+ll)) = B4(idi(II), idi(II+ll))+ sin(chi)*ephiM;
     end

     if II-ll >= 0
         B4(idr(II), idr(II-ll)) = B4(idr(II),idr(II-ll)) + sin(chi)*ephiP;
         B4(idi(II), idi(II-ll)) = B4(idi(II),idi(II-ll)) + sin(chi)*ephiP;
     else
         B4(idr(II), idr(-(II-ll))) = B4(idr(II), idr(-(II-ll))) + sin(chi)*ephiP;
         B4(idi(II), idi(-(II-ll))) = B4(idi(II), idi(-(II-ll))) - sin(chi)*ephiP;
     end

    end
%     for II = 0 : N
%         if II+1 <= N
%             idr = 2*(II) + 1;
%             idi = 2*(II+1);
%             B4(idr, idr+2) = B4(idr, idr+2) + 1;
%             B4(idi, idi+2) = B4(idi, idi+2) + 1;
%         end
% 
%         if II-1 >= 0
%             idr = 2*(II) + 1;
%             idi = 2*(II+1);
%             B4(idr, idr-2) = B4(idr, idr-2) + 1;
%             B4(idi, idi-2) = B4(idi, idi-2) + 1;
%         else
%             idr = 2*(II) + 1;
%             idi = 2*(II+1);
%             B4(idr,idr+2) = B4(idr, idr+2) + 1;
%             B4(idi, idi+2) = B4(idi, idi+2) - 1;
%         end
%     end


    B4 = B4*0.5*(rho2-rho1)*k^2;
    AB = A\B4;
    EV = eig(AB);
    amplitude = 1./EV;
    amplitude = amplitude(amplitude>0);
    amplitude = amplitude(amplitude<amax);
    amplitude = amplitude(imag(amplitude)==0);
    for ii = 1 : length(amplitude)
        harmonicData = [harmonicData; amplitude(ii)/g, k/1000];
    end
end

scatter(subharmonicData(:,2),subharmonicData(:,1),'filled','blue');
hold on;
scatter(harmonicData(:,2),harmonicData(:,1),'filled','red');
xlabel('k (mm^{-1})');
ylabel('a/g');
ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'Bold';
ax.LineWidth = 2;
box on;

function id = idr(ii)
    id = 2*ii + 1;
end

function id = idi(ii)
    id = 2*(ii+1);
end