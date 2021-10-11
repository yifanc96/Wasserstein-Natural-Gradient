function [ W2,T,varphi ] = one_dim_OT(T_start,T_end,dt,ff,gg)
    tt=T_start:dt:T_end;
    Nt=length(tt);
    %target
    GG=zeros(1,Nt);
    for i=1:Nt
        GG(i)=sum(gg(1:i)*dt)-1/2*gg(1)*dt-1/2*gg(i)*dt;
    end
    gg=gg./max(GG);
    GG=GG./max(GG);
    down_line=min(GG);
    up_line=max(GG);

    %initial
    FF=zeros(1,Nt);
    for i=1:Nt
        FF(i)=sum(ff(1:i))*dt-1/2*ff(1)*dt-1/2*ff(i)*dt;
    end
    ff=ff./max(FF);
    FF=FF./max(FF);
    T=zeros(1,Nt); %optimal transport map G^{-1}(F)
    for it=1:Nt
        if FF(it) <= down_line
            T(it)=T_start;
        elseif FF(it) >= up_line
            T(it)=T_end;
        else 
            place=find(GG(:)>=FF(it),1,'first');
            if place==1
                T(it)=T_start;
            else 
                T(it)=T_start-dt+(place-1)*dt+(FF(it)-GG(place-1))*(dt)/(GG(place)-GG(place-1)); %linear interpolation
            end
        end 
    end

    W2=sum((tt-T).^2.*ff)*dt-1/2*(tt(1)-T(1)).^2*ff(1)*dt-1/2*(tt(Nt)-T(Nt)).^2*ff(Nt)*dt;
    varphi=zeros(1,Nt);
    temp=tt-T;
    for j=1:Nt
        varphi(j)=2*sum(temp(1:j))*dt-temp(1)*dt-temp(j)*dt;
    end

end

