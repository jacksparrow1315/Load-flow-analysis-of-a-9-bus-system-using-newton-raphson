function [Pi, Qi, Pg, Qg, Pl, Ql] = calculation(Nbus,V,del,baseMVA)
Y = Y_bus(Nbus);
lined = Line_Data(Nbus);
busd = Bus_Data(Nbus);
Vm = pol2rec(V,del);
Del = del*(180/pi);
fb = lined(:,1);
tb = lined(:,2);
nl = length(fb);
Pl = busd(:,7);
Ql = busd(:,8);

Iij = zeros(Nbus,Nbus);
Sij = zeros(Nbus,Nbus);
Si  = zeros(Nbus,1);

%%Injecting Bus current
I = Y*Vm;
Im = abs(I);
Ia = angle(I);

%%Line Current
for m = 1:nl
    p = fb(m);
    q = tb(m);
    Iij(p,q) = -(Vm(p) - Vm(q))*Y(p,q);
    Iij(q,p) = -Iij(p,q);
end

Iijm = abs(Iij);
Iija = angle(Iij);


%%Line power
for m = 1:Nbus
    for n = 1:Nbus
        if m ~= n
            Sij(m,n) = Vm(m)*conj(Iij(m,n))*baseMVA;
        end
    end
end

Pij = real(Sij);
Qij = imag(Sij);

%%Line Loss
Lij = zeros(nl,1);
for m = 1:nl
    p = fb(m);
    q = tb(m);
    Lij(m) = Sij(p,q) + Sij(q,p);
end

Lpij = real(Lij);
Lqij = imag(Lij);

%%Bus power Injections
for i = 1:Nbus
    for j = 1 : Nbus
        Si(i) = Si(i) + conj(Vm(i))*Vm(j)*Y(i,j)*baseMVA;
    end
end

Pi = real(Si);
Qi = - imag(Si);
Pg = Pi + Pl;
Qg = Qi + Ql;


disp('                                  Newton Raphson Load Flow Analysis');
disp('--------------------------------------------------------------------------------------------------------------');
disp('|  Bus   |     V      |     Angle   |    Injected Power     |    Generated Power    |         Load        |');
disp('|  no.   |    pu      |     Degree  |   MW      |   MVar    |   MW      |    MVar   |     MW    |   MVar  |');

for m = 1:Nbus
    disp('----------------------------------------------------------------------------------------------------------');
    fprintf('%3g',m); fprintf('       %8.4f',V(m)); fprintf('     %8.4f',Del(m)); fprintf('     %8.4f',Pi(m)); 
    fprintf('     %8.4f',Qi(m)); fprintf('     %8.4f',Pg(m)); fprintf('     %8.4f',Qg(m)); 
    fprintf('    %8.4f',Pl(m)); fprintf('   %8.4f',Ql(m)); fprintf('\n');
end
disp('..............................................................................................................');
fprintf('Total                          '); fprintf('     %8.4f',sum(Pi)); fprintf('      %8.4f',sum(Qi)); 
fprintf('    %8.4f',sum(Pg)); fprintf('     %8.4f',sum(Qg)); fprintf('    %8.4f',sum(Pl)); fprintf('   %8.4f',sum(Ql));
fprintf('\n\n\n');

disp('                                         Line Flow And Losses');
disp('--------------------------------------------------------------------------------------------------------------');
disp('|  From  |  To    |    P     |    Q       | |  From  |  To    |    P     |    Q       |        Line Loss         |');
disp('|  Bus   |  Bus   |    MW    |    MVar    | |  Bus   |  Bus   |    MW    |    MVar    |     MW      |      MVar  |');
for m = 1:nl
    p = fb(m); q = tb(m);
    disp('----------------------------------------------------------------------------------------------------------');
    fprintf('  %4g',p); fprintf(  '%7g',q); fprintf('      %8.4f',Pij(p,q));  fprintf('   %8.4f',Qij(p,q));
    fprintf('  %6g',q); fprintf('  %7g',p); fprintf('      %8.4f',Pij(q,p));  fprintf('      %8.4f',Qij(q,p));
    fprintf('    %8.4f',Lpij(m)); fprintf('    %8.4f',Lqij(m)); 
    fprintf('\n');
end
disp('..............................................................................................................');
fprintf('Total Loss                                                                           ');
fprintf('  %8.4f',sum(Lpij)); fprintf('    %8.4f',sum(Lqij));
fprintf('\n\n\n');

LL = lined(:,7);

disp('                                   Overloaded Lines and Cables');
disp('------------------------------------------------------------------------------------------------------------');
disp('|  From  |  To    |    P     |    Q       |     S      |  Loading Limit |');
disp('|  Bus   |  Bus   |    MW    |    MVar    |    Mvar    |      Mvar      |');

for m = 1:nl
    p = fb(m); q = tb(m);
    temp = abs(Sij(p,q));
    if temp > LL(m)
        disp('----------------------------------------------------------------------------------------------------------');
        fprintf('  %4g',p); fprintf(  '%7g',q); fprintf('      %8.4f',Pij(p,q));  fprintf('    %8.4f',Qij(p,q));
        fprintf('      %8.4f',abs(Sij(p,q)));  fprintf('       %8.4f',LL(m));
        fprintf('\n');
    end
end


fprintf('\n\n\n');
disp('                                    Underloaded Lines and Cables');
disp('------------------------------------------------------------------------------------------------------------');
disp('|  From  |  To    |    P     |    Q       |      S      |  Loading Limit |');
disp('|  Bus   |  Bus   |    MW    |    MVar    |     Mvar    |      Mvar      |');

for m = 1:nl
    p = fb(m); q = tb(m);
    temp = abs(Sij(p,q));
    if temp < 0.5*LL(m)
        disp('----------------------------------------------------------------------------------------------------------');
        fprintf('  %4g',p); fprintf(  '%7g',q); fprintf('      %8.4f',Pij(p,q));  fprintf('    %8.4f',Qij(p,q));
        fprintf('      %8.4f',abs(Sij(p,q)));   fprintf('       %8.4f',LL(m));
        fprintf('\n');
    end
end

end