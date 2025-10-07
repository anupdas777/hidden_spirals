
function eigen_value= Jacobian(inp_frame,c,d,k)
siz=size(inp_frame);
dia=siz(1,1);

hp=@(x) cos(x)-c*sin(x)-d*2*sin(2*x)+2*k*cos(2*x);



frm=inp_frame;
StaMa=zeros(dia^2,dia^2);

for i=1:dia^2
   for j=1:dia^2

    ri=min((i-mod(i,dia))/dia+1,dia);

    if mod(i,dia)==0
    ci=dia;
    else
        ci=mod(i,dia);
    end

    rj=min((j-mod(j,dia))/dia+1,dia);

     if mod(j,dia)==0
       cj=dia;
    else
        cj=mod(j,dia);
    end
    
    if abs(ri-rj)+abs(ci-cj)==1
     exi=1;
    else
     exi=0;
    end

    StaMa(i,j)=exi*hp((1)*(frm(rj,cj)-frm(ri,ci)));

    if i==j
        StaMa(i,j)=0;
    end
end

StaMa(i,i)=-sum(StaMa(i,:));


end
eigen_value=sort(real(eig(StaMa)));
end







