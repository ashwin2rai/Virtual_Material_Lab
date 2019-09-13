function [ Cmat ] = Stiffness(E_mat,E1_fib,E2_fib,nu12_fib,vffib,nu_mat,G12_fib)
% E_mat=3*10^9;
% E1_fib=110*10^9
% E2_fib=8*10^9
% nu12_fib=0.23
% vffib=0.57
% nu_mat=0.3
% 
% G12_fib=5*10^9


G12_mat=E_mat/(2*(1+nu_mat));
vfmat=1-vffib;

% Reuss and Voight models to homogenize and find material properties of the
% composite

E1=E1_fib*vffib+E_mat*vfmat;
E2=1/((vffib/E2_fib)+(vfmat/E_mat));
v12=vffib*nu12_fib+vfmat*nu_mat;
v21=1/(vffib/nu12_fib+vfmat/nu_mat);
G12=1/((vffib/G12_fib)+vfmat/G12_mat)

S=[(1/E1) (-v12/E2) 0; -v12/E2 1/E2 0; 0 0 (1/G12)];
C=S^-1
Cmat=C;

end

