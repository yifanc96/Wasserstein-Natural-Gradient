function Gm = computeGm_2mixture(T_start,T_end,dt,theta,T)

    tt=T_start:dt:T_end;
    Nt=length(tt);
    ff=mixture_density(theta,tt);
            eps=1e-5;
            F_grad=zeros(5,Nt);
            F_grad(1,:)=cdf('Normal',tt,theta(2),sqrt(theta(3)))-cdf('Normal',tt,theta(4),sqrt(theta(5)));
            F_grad(2,:)=(theta(1).*cdf('Normal',tt,theta(2)+eps,sqrt(theta(3)))-theta(1).*cdf('Normal',tt,theta(2)-eps,sqrt(theta(3))))./(2*eps);
            F_grad(3,:)=(theta(1).*cdf('Normal',tt,theta(2),sqrt(theta(3)+eps))-theta(1).*cdf('Normal',tt,theta(2),sqrt(theta(3)-eps)))./(2*eps);
            F_grad(4,:)=((1-theta(1)).*cdf('Normal',tt,theta(4)+eps,sqrt(theta(5)))-(1-theta(1)).*cdf('Normal',tt,theta(4)-eps,sqrt(theta(5))))./(2*eps);
            F_grad(5,:)=((1-theta(1)).*cdf('Normal',tt,theta(4),sqrt(theta(5)+eps))-(1-theta(1)).*cdf('Normal',tt,theta(4),sqrt(theta(5)-eps)))./(2*eps);
    
            temp=F_grad;
            %T'
            dT=zeros(1,length(T));
            dT(1)=1;dT(length(T))=1;
            ratio=1;
            for ii=2:length(T)-1
                dT(ii)=min([max([ratio*(T(ii+1)-T(ii-1))/(2*dt)+(1-ratio),0]),20]);   
            end
    
            for i=1:5
                temp(i,:)=dT.*F_grad(i,:)./ff;
                loc=find(abs(ff)<1e-4);
                temp(i,loc)=0;
            end
            Gm=temp*F_grad'*dt;

end

