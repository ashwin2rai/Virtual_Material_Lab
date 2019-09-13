clc; clear; close all;
%%%%%%%%%%%%%%%%%%%Displacement Field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ux = A*x + B*y
% uy = C*x + D*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'test2.gif';   %save the resulting gif with this filename
Choice = [ 0 0 0 1];
% Choice = [1 0 0 0] is uniaxial strain in the x direction
% Choice = [0 1 0 0] is uniaxial strain in the y direction
% Choice = [0 0 1 0] is pure shear strain in the xy direction
% Choice = [0 0 0 1] is any combination of strains, see below
%
% if you use choice = [0 0 0 1] please enter known_vec
% Enter any combination of strains and stress, please note only 3 values
% can be provided otherwise the system is overconstrained
% for example, if you want to apply strain in 11 and 12 (shear), then enter
% known_vec = [1 3 5] which means strain applied in 11, 12 and stress in
% 22 is 0 (unconstrained).
% another example, uniaxial strain in 11 would be known_vec = [1 5 6]
% strain applied in x direction (1) and stress in y direction (5) and shear
% direction (6) is free
%
% known_vec : eps_xx is 1 eps_yy is 2 eps_xy is 3 sigma_xx is 4 
% sigma_yy is 5  sigma_xy is 6
%
known_vec = [ 1  2 3 ]; 
gif=1;
time_step=25;
percentage_strain=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material properties: in this example, an isotropic material
% But feel free to change these properties to see how non isotropic 
% materials behave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod=100e9;
nu=0.3;

C11 = mod/(1-nu^2);
C12 = nu*C11;
C13 = 0;
C21 = C12;
C22 = C11;
C23 = 0;
C13 = 0;
C23 = 0;
C33 = C11*(1-nu)/2;

C_mat = [ C11  C12 C13 ; C12 C22 C23; C13 C23 C33 ];

hold on;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Node numbering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |
%        |
%        7------6------5
%        |             |
%        |             |
%        8      9      4
%        |             |
%        |             |
%        1------2------3------

Nodes.x=[0 1 2 2 2 1 0 0 1];
Nodes.y=[0 0 0 1 2 2 2 1 1];

hold on;
fill(Nodes.x(1:8),Nodes.y(1:8),'r');
axis equal;
if gif ==1
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if Choice(1) == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option for entering strains in x direction in per centage strain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ctr=1;
    eps_xx_p = percentage_strain;   %percentage strain
    total_time=1;    %pseudotime 1 sec default
    
    sigma_yy(ctr) = 0;    %boundary conditions for normal strains in x direction
    sigma_xy(ctr) = 0;
    
    delta_t=total_time/time_step;
    eps_xx = eps_xx_p/100;
    eps_xx_dot=eps_xx/total_time;     %strain rate
    
    A=zeros(time_step+1,1);
    B=zeros(time_step+1,1);
    C=zeros(time_step+1,1);
    A_dot=zeros(time_step,1);
    B_dot=zeros(time_step,1);
    C_dot=zeros(time_step,1);
    eps_xx_cur=zeros(time_step+1,1);     %initialization
    sigma_xx = zeros(time_step+1,1);
    sigma_xx_dot = zeros(time_step,1);
    
    for i=1:time_step+1
        ux(i).nodes=zeros(length(Nodes.x),1);
        uy(i).nodes=zeros(length(Nodes.x),1);
    end
    
    for i=1:time_step
        
        A_dot(ctr) = eps_xx_dot;
        newvals=-(C_mat(2:3,1)*A_dot(ctr))'/C_mat(2:3,2:3);
        C_dot(ctr)=newvals(1);
        B_dot(ctr)=newvals(2)/2;
        eps_xx_cur(ctr+1)=eps_xx_cur(ctr)+eps_xx_dot*delta_t;
        A(ctr+1)=A(ctr)+A_dot(ctr)*delta_t;
        B(ctr+1)=B(ctr)+B_dot(ctr)*delta_t;
        C(ctr+1)=C(ctr)+C_dot(ctr)*delta_t;
        sigma_xx_dot(ctr) = C11*A_dot(ctr) + C12*C_dot(ctr) + 2*B_dot(ctr)*C13;
        sigma_xx(ctr+1) = sigma_xx(ctr) + sigma_xx_dot(ctr)*delta_t;
        ux(ctr+1).nodes = A(ctr+1).*Nodes.x + B(ctr+1).*Nodes.y;
        uy(ctr+1).nodes=  B(ctr+1).*Nodes.x + C(ctr+1).*Nodes.y;
        NewNodes.x=Nodes.x+ux(ctr+1).nodes;
        NewNodes.y=Nodes.y+uy(ctr+1).nodes;
        figure(1);
        fill(NewNodes.x(1:8),NewNodes.y(1:8),'c'); hold on; axis equal;
        if gif==1
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1);
        end
        ctr=ctr+1;
    end
    figure(2);
    hold on;
    plot(eps_xx_cur,sigma_xx,'LineWidth',2);
end

if Choice(2) == 1
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option for entering strains in y direction in per centage strain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ctr=1;
    eps_yy_p = percentage_strain;   %percentage strain
    total_time=1;    %pseudotime 1 sec default
    
    sigma_xx(ctr) = 0;    %boundary conditions for normal strains in x direction
    sigma_xy(ctr) = 0;
    
    delta_t=total_time/time_step;
    eps_yy = eps_yy_p/100;
    eps_yy_dot=eps_yy/total_time;     %strain
    
    A=zeros(time_step+1,1);
    B=zeros(time_step+1,1);
    C=zeros(time_step+1,1);
    A_dot=zeros(time_step,1);
    B_dot=zeros(time_step,1);
    C_dot=zeros(time_step,1);
    eps_yy_cur=zeros(time_step+1,1);     %initialization
    sigma_yy = zeros(time_step+1,1);
    sigma_yy_dot = zeros(time_step,1);
    
    for i=1:time_step+1
        ux(i).nodes=zeros(length(Nodes.x),1);
        uy(i).nodes=zeros(length(Nodes.x),1);
    end
    
    for i=1:time_step
        
        C_dot(ctr) = eps_yy_dot;
        newvals=-(C_mat(1:2:3,2)*C_dot(ctr))'/C_mat(1:2:3,1:2:3);
        A_dot(ctr)=newvals(1);
        B_dot(ctr)=newvals(2)/2;
        eps_yy_cur(ctr+1)=eps_yy_cur(ctr)+eps_yy_dot*delta_t;
        A(ctr+1)=A(ctr)+A_dot(ctr)*delta_t;
        B(ctr+1)=B(ctr)+B_dot(ctr)*delta_t;
        C(ctr+1)=C(ctr)+C_dot(ctr)*delta_t;
        sigma_yy_dot(ctr) = C12*A_dot(ctr) + C22*C_dot(ctr) + 2*B_dot(ctr)*C23;
        sigma_yy(ctr+1) = sigma_yy(ctr) + sigma_yy_dot(ctr)*delta_t;
        ux(ctr+1).nodes = A(ctr+1).*Nodes.x + B(ctr+1).*Nodes.y;
        uy(ctr+1).nodes=  B(ctr+1).*Nodes.x + C(ctr+1).*Nodes.y;
        NewNodes.x=Nodes.x+ux(ctr+1).nodes;
        NewNodes.y=Nodes.y+uy(ctr+1).nodes;
        figure(1);
        fill(NewNodes.x(1:8),NewNodes.y(1:8),'g'); hold on; axis equal;
        if gif==1
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1);
        end
        ctr=ctr+1;
    end
    figure(2);
    hold on;
    plot(eps_yy_cur,sigma_yy,'--','LineWidth',2);
end

if Choice(3) == 1
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option for entering strains in xy direction in per centage strain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ctr=1;
    eps_xy_p = percentage_strain;   %percentage strain
    total_time=1;    %pseudotime 1 sec default
    
    sigma_xx(ctr) = 0;    %boundary conditions for normal strains in x direction
    sigma_yy(ctr) = 0;
    
    delta_t=total_time/time_step;
    eps_xy = eps_xy_p/100;
    eps_xy_dot=eps_xy/total_time;     %strain rate
    
    A=zeros(time_step+1,1);
    B=zeros(time_step+1,1);
    C=zeros(time_step+1,1);
    A_dot=zeros(time_step,1);
    B_dot=zeros(time_step,1);
    C_dot=zeros(time_step,1);
    eps_xy_cur=zeros(time_step+1,1);     %initialization
    sigma_xy = zeros(time_step+1,1);
    sigma_xy_dot = zeros(time_step,1);
    
    for i=1:time_step+1
        ux(i).nodes=zeros(length(Nodes.x),1);
        uy(i).nodes=zeros(length(Nodes.x),1);
    end
    
    for i=1:time_step
        
        B_dot(ctr) = eps_xy_dot;
        newvals=-(C_mat(1:2,3)*B_dot(ctr))'/C_mat(1:2,1:2);
        A_dot(ctr)=newvals(1);
        C_dot(ctr)=newvals(2);
        eps_xy_cur(ctr+1)=eps_xy_cur(ctr)+eps_xy_dot*delta_t;
        A(ctr+1)=A(ctr)+A_dot(ctr)*delta_t;
        B(ctr+1)=B(ctr)+B_dot(ctr)*delta_t;
        C(ctr+1)=C(ctr)+C_dot(ctr)*delta_t;
        sigma_xy_dot(ctr) = C13*A_dot(ctr) + C23*C_dot(ctr) + 2*B_dot(ctr)*C33;
        sigma_xy(ctr+1) = sigma_xy(ctr) + sigma_xy_dot(ctr)*delta_t;
        ux(ctr+1).nodes = A(ctr+1).*Nodes.x + B(ctr+1).*Nodes.y;
        uy(ctr+1).nodes=  B(ctr+1).*Nodes.x + C(ctr+1).*Nodes.y;
        NewNodes.x=Nodes.x+ux(ctr+1).nodes;
        NewNodes.y=Nodes.y+uy(ctr+1).nodes;
        figure(1);
        fill(NewNodes.x(1:8),NewNodes.y(1:8),'y'); hold on; axis equal;
        if gif==1
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1);
        end
        ctr=ctr+1;
    end
    figure(2);
    hold on;
    plot(eps_xy_cur,sigma_xy,'LineWidth',2);
    xlabel('Strain','Interpreter','tex');
    ylabel('Stress (Pa)','Interpreter','tex');
end

if Choice(4) == 1
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option for entering any combination of strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ctr=1;
    
    eps_xx_p = percentage_strain;
    eps_yy_p = percentage_strain/2;
    eps_xy_p = percentage_strain/5;
    sigma_xx(ctr) = 0;
    sigma_yy(ctr) = 0;
    sigma_xy(ctr) = 0;
    
    % Enter any combination strains and stress, please note only 3 values
    % can be provided otherwise the system is overconstrained
    % for example, if you want to apply strain in 11 and 12, then enter
    % known_vec = [1 3 5] which means strain applied in 11 12 and stress in
    % 22 is 0 (unconstrained)
    
    % known_vec = [eps_xx eps_yy eps_xy sigma_xx sigma_yy sigma_xy]
    
    unknown_vec = setdiff(1:6,known_vec);
    comp_C_mat = [C_mat [1 0 0; 0 1 0; 0 0 1]];
    
    
    total_time=1;    %pseudotime 1 sec default
    eps_xx_dot = 0;
    eps_yy_dot = 0;
    eps_xy_dot = 0;
    
    A=zeros(time_step+1,1);
    B=zeros(time_step+1,1);
    C=zeros(time_step+1,1);
    A_dot=zeros(time_step,1);
    B_dot=zeros(time_step,1);
    C_dot=zeros(time_step,1);
    
    eps_xx_cur=zeros(time_step+1,1);     %initialization
    sigma_xx = zeros(time_step+1,1);
    sigma_xx_dot = zeros(time_step,1);
    eps_yy_cur=zeros(time_step+1,1);     %initialization
    sigma_yy = zeros(time_step+1,1);
    sigma_yy_dot = zeros(time_step,1);
    eps_xy_cur=zeros(time_step+1,1);     %initialization
    sigma_xy = zeros(time_step+1,1);
    sigma_xy_dot = zeros(time_step,1);
    
    if ~isempty(find(known_vec==1))
        delta_t=total_time/time_step;
        eps_xx = eps_xx_p/100;
        eps_xx_dot=eps_xx/total_time;     %strain rate
    end
    
    if ~isempty(find(known_vec==2))
        delta_t=total_time/time_step;
        eps_yy = eps_yy_p/100;
        eps_yy_dot=eps_yy/total_time;     %strain
    end
    
    if ~isempty(find(known_vec==3))
        delta_t=total_time/time_step;
        eps_xy = eps_xy_p/100;
        eps_xy_dot=eps_xy/total_time;     %strain rate
    end
    con_vec = [eps_xx_dot eps_yy_dot eps_xy_dot 0 0 0 ]';
    
    for i=1:time_step+1
        ux(i).nodes=zeros(length(Nodes.x),1);
        uy(i).nodes=zeros(length(Nodes.x),1);
    end
    
    for i=1:time_step
        
        A_dot(ctr) = eps_xx_dot;
        C_dot(ctr) = eps_yy_dot;
        B_dot(ctr) = eps_xy_dot;
        
        newvals = comp_C_mat(:,unknown_vec)\(-comp_C_mat(:,known_vec)*con_vec(known_vec));
        newval_vec(unknown_vec) = newvals;
        newval_vec(known_vec) = con_vec(known_vec);
        
        A_dot(ctr)=newval_vec(1);
        B_dot(ctr)=newval_vec(3);
        C_dot(ctr)=newval_vec(2);
        
        eps_xy_cur(ctr+1)=eps_xy_cur(ctr)+eps_xy_dot*delta_t;
        eps_xx_cur(ctr+1)=eps_xx_cur(ctr)+eps_xx_dot*delta_t;
        eps_yy_cur(ctr+1)=eps_yy_cur(ctr)+eps_yy_dot*delta_t;
        A(ctr+1)=A(ctr)+A_dot(ctr)*delta_t;
        B(ctr+1)=B(ctr)+B_dot(ctr)*delta_t;
        C(ctr+1)=C(ctr)+C_dot(ctr)*delta_t;
        []
        sigma_xx_dot(ctr) = C11*A_dot(ctr) + C12*C_dot(ctr) + 2*B_dot(ctr)*C13;
        sigma_xx(ctr+1) = sigma_xx(ctr) + sigma_xx_dot(ctr)*delta_t;
        sigma_yy_dot(ctr) = C12*A_dot(ctr) + C22*C_dot(ctr) + 2*B_dot(ctr)*C23;
        sigma_yy(ctr+1) = sigma_yy(ctr) + sigma_yy_dot(ctr)*delta_t;
        sigma_xy_dot(ctr) = C13*A_dot(ctr) + C23*C_dot(ctr) + 2*B_dot(ctr)*C33;
        sigma_xy(ctr+1) = sigma_xy(ctr) + sigma_xy_dot(ctr)*delta_t;
        ux(ctr+1).nodes = A(ctr+1).*Nodes.x + B(ctr+1).*Nodes.y;
        uy(ctr+1).nodes=  B(ctr+1).*Nodes.x + C(ctr+1).*Nodes.y;
        NewNodes.x=Nodes.x+ux(ctr+1).nodes;
        NewNodes.y=Nodes.y+uy(ctr+1).nodes;
        figure(1);
        fill(NewNodes.x(1:8),NewNodes.y(1:8),'m'); hold on; axis equal;
        if gif==1
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1);
        end
        ctr=ctr+1;
    end
    figure(2);
    hold on;
    plot(eps_xx_cur,sigma_xx,'LineWidth',2);
    plot(eps_yy_cur,sigma_yy,'--','LineWidth',2);
    plot(eps_xy_cur,sigma_xy,'LineWidth',2);
    xlabel('Strain','Interpreter','tex');
    ylabel('Stress (Pa)','Interpreter','tex');
end

Nodes.x=[0 1 2 2 2 1 0 0 1];
Nodes.y=[0 0 0 1 2 2 2 1 1];
figure(1);
hold on;
fill(Nodes.x(1:8),Nodes.y(1:8),'r');
axis equal;
figure(1);
if gif==1
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
end
