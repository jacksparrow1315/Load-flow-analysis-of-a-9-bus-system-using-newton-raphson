function x = Newton_Raphson(nbus,Y,busd,baseMVA)
bus = busd(:,1);
type = busd(:,2);
V = busd(:,3);
del = busd(:,4);
Pg = busd(:,5)/baseMVA;
Qg = busd(:,6)/baseMVA;
Pl = busd(:,7)/baseMVA;
Ql = busd(:,8)/baseMVA;
Qmin = busd(:,9)/baseMVA;
Qmax = busd(:,10)/baseMVA;

P = Pg - Pl;
Q = Qg - Ql;
Psp = P;
Qsp = Q;
G = real(Y);
B = imag(Y);

pv = find(type == 2 | type == 1);
pq = find(type == 3);
npv = length(pv);
npq = length(pq);

Tol = 1;
Iter = 1;
while (Tol>0.0000001)
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)*V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)*V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    %%Q limit
    if Iter<=7 && Iter>2
        for n = 2:nbus
            if type(n)==2
                QG = Q(n) + Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
        end
    end

    %%Change from specified vaule
    dPa = Psp-P;
    dQa = Qsp-Q;

    k=1;
    dQ = zeros(npq,1);
    for i =1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k= k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];    %% Mismatch Vector


    %%Jacobians

    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k)+ V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k)-V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)*V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end


    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end


    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)*V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k)-V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end


    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1 :nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end


    J = [J1 J2; J3 J4];     %%Final jacobian matrix

    X = J\M;                 %%Correcting
    dTh = X(1:(nbus-1));  %%Corected Voltage angle
    dV = X(nbus:end);       %%Corrected Volatage Magnitude

    %%Correction
    del(2:nbus) = del(2:nbus) + dTh;

    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = V(i) + dV(k);
            k = k+1;
        end
    end

    Iter = Iter + 1;
    Tol = max (abs(M));

end

calculation(nbus,V,del,baseMVA);

fprintf("\nIteration needed::")
disp(Iter);
% fprintf('\n');
% fprintf('Jacobians:\n');
% fprintf("J1::\n")
% disp(J1);
% fprintf("J2::\n");
% disp(J2);
% fprintf("J3::\n");
% disp(J3);
% fprintf("J4::\n");
% disp(J4);
fprintf('\nBus admittance matrix:\n');
disp(Y)
fprintf('\nJacobian::\n');
disp(J)

subplot(211); bar(V);
title("Bus Voltage Magnitude"); xlabel("Bus Number"); ylabel('|V| pu');

subplot(212); bar(del*(180/pi));
title("Bus Voltage Angle in Degree"); xlabel("Bus Number"); ylabel('Voltage angle');


end