function [ f_grad ] = mixture_gradient( theta,tt )
eps=1e-5;
        Nt=length(tt);
        f_grad=zeros(5,Nt);    
        f_grad(1,:)=pdf('Normal',tt,theta(2),sqrt(theta(3)))-pdf('Normal',tt,theta(4),sqrt(theta(5)));
        f_grad(2,:)=(theta(1).*pdf('Normal',tt,theta(2)+eps,sqrt(theta(3)))-theta(1).*pdf('Normal',tt,theta(2)-eps,sqrt(theta(3))))./(2*eps);
        f_grad(3,:)=(theta(1).*pdf('Normal',tt,theta(2),sqrt(theta(3)+eps))-theta(1).*pdf('Normal',tt,theta(2),sqrt(theta(3)-eps)))./(2*eps);
        f_grad(4,:)=((1-theta(1)).*pdf('Normal',tt,theta(4)+eps,sqrt(theta(5)))-(1-theta(1)).*pdf('Normal',tt,theta(4)-eps,sqrt(theta(5))))./(2*eps);
        f_grad(5,:)=((1-theta(1)).*pdf('Normal',tt,theta(4),sqrt(theta(5)+eps))-(1-theta(1)).*pdf('Normal',tt,theta(4),sqrt(theta(5)-eps)))./(2*eps);
end

