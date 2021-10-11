
% Wasserstein metric based mixture fitting: 2-component mixture family
% 'theta_initial' is the initial parameter for the mixture family 

function [iteration_hist, density_figure ] = mixture_W2obj_minimization( T_start, T_end, dt, theta_initial, density_target,method )
%% time window
    tt=T_start:dt:T_end;

%% target distribution
    gg=density_target;  

%% initial distribution by theta_initial
    ff=mixture_density(theta_initial,tt);

%% optimization
    step=1;
    %stop criteria
        stepmax=300;
        norm_of_grad=1e-2;
        line_serach_step_min=1e-4;
    
    cost=[];  % objective value
    while(1)
        if step==1
            thetanow=theta_initial;
            theta=thetanow;
            [W2,T,varphi]=one_dim_OT(T_start,T_end,dt,ff,gg);
        end
        cost=[cost,W2];
    
        %Gradient computation
        density_grad=mixture_gradient(theta,tt);
        grad=zeros(1,5);
        for ii=1:5
            grad(ii)=sum(density_grad(ii,:).*varphi)*dt;
        end
    
    %-----------------------Computation of preconditioning matrix--------------
    switch method
        case {'WGD'}
            G=computeG_2mixture(T_start,T_end,dt,theta);  
            Hv=inv(G)*0.5;
        case {'FRGD'}
            G=computeF_2mixture(T_start,T_end,dt,theta);  
            Hv=inv(G);
        case {'GD'}
            G=eye(5);  
            Hv=inv(G);
        case {'pGD'}
            G=eye(5);
            G(1,1)=30;
            Hv=inv(G);
        case {'WGDm'} %modified WGD
            Gm=computeGm_2mixture(T_start,T_end,dt,theta,T);
            Hv=inv(Gm)*0.5;
    end

    %line search
    stepeps=line_serach_step_min; %lower bound of stepsize
    ss=1; %initial stepsize
    
    W2_old=W2;
       while(1)
           thetanow=theta-ss*grad*Hv;
            %parameter constraints
            while (max(isnan(thetanow))>0.5 ||thetanow(1)>=1 || thetanow(1) <=0 || thetanow(2)<-13 || thetanow(2)>13 || thetanow(3)<=0.01 || thetanow(3)>30 || thetanow(4)<=-13 || thetanow(4)>=13 || thetanow(5)<0.01 || thetanow(5)>30)
                ss=ss/2;
                if ss<stepeps    
                    fprintf('out of the region\n');
                    break
                end
                thetanow=theta-ss*grad*Hv;
            end
            if thetanow(1)>=1 
                thetanow(1)=1;
            elseif thetanow(1)<=0
                thetanow(1)=0;
            end
        
            if ss<stepeps
                break
            end
            [W2,T,varphi]=one_dim_OT(T_start,T_end,dt,mixture_density(thetanow,tt),gg);

            if W2_old<=W2
                ss=ss/2;
                if ss<stepeps
                    fprintf('line search failed\n');
                end
            else
                disp(['Stepsize ss =',num2str(ss)]);
                break
            end
       end
       
        %updating
        theta=thetanow;
        disp(['Step=',num2str(step),' and the result now is ',num2str(theta)]);
        step=step+1;
        if step>=stepmax || norm(grad,2)<norm_of_grad || ss<stepeps
            cost=[cost,W2];
            break
        end
    end


    density_figure=figure(1);
    plot(tt,mixture_density(theta_initial,tt));
    hold on
    plot(tt,gg);
    plot(tt,mixture_density(theta,tt));
    legend('initial','true','converged');
    title('figure of densities');

    iteration_hist=figure(2);
    semilogy(0:step-1,cost);
    title('iteration history');
    fprintf('objective value is %e \n',W2)
    ylabel('W^2 distance');xlabel('iteration step');

end

