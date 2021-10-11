function [ output ] = mixture_density( theta,x )
output= theta(1).*pdf('Normal',x,theta(2),sqrt(theta(3)))+(1-theta(1)).*pdf('Normal',x,theta(4),sqrt(theta(5)));

end

