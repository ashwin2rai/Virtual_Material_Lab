function [E1_fib,E_mat] = localizer(eps_xx,eps_yy,eps_xy,E_mat,E1_fib,E2_fib,nu12_fib,vffib,nu_mat,G12_fib)

Sfib=[1/E1_fib -nu12_fib/E2_fib 0; -nu12_fib/E2_fib 1/E2_fib 0; 0 0 G12_fib];
Smat=[1/E_mat -nu_mat/E_mat 0; -nu_mat/E_mat 1/E_mat 0; 0 0 E_mat/(2*(1+nu_mat))];

Cfib=inv(Sfib);
Cmtrx=inv(Smat);

sgmfib=Cfib*[eps_xx; eps_yy; eps_xy];
sgmmtrx=Cmtrx*[eps_xx; eps_yy; eps_xy];

Plotfibmat(sgmmtrx(1),0,3e10,sgmfib(1),0,3e10,3)
pause(0.1)

% Very simple max stress failure criterion for fiber

if sgmfib(1)>=1e10;
    E1_fib=1 
end

% Very simple strain dependent damage criterion for matrix
E_mat=E_mat*(1-eps_yy);

end

