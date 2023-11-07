function[F_i_cortex]=Cortical_Dynein_Force(LengthFac,xRadius,nc,MT_state,c_attach,centers,centers_old,dt,vf,f0_dynein,v0_dynein,mt_vec,MT_length);

MT_bound_cortex=find(MT_state==4);
Centrosome_MT_bound_cortex=c_attach(MT_bound_cortex);
F_i_cortex=zeros(2,nc);
for i=1:nc
%    magvel(1,i)=sqrt((centers(1,i)-centers_old(1,i)).^2+(centers(2,i)-centers_old(2,i)).^2)/dt;
vel(:,i)=[(centers(1,i)-centers_old(1,i));(centers(2,i)-centers_old(2,i))]/dt;
end
%vel=[magvel; magvel];
if(isempty(MT_bound_cortex)==0)
for i=1:nc
    mindistcortex=xRadius-sqrt(centers(1,i).^2+centers(2,i).^2);
    MTBoundOni=find(Centrosome_MT_bound_cortex==i);
    if(isempty(MTBoundOni)==0)
        vd=sum(mt_vec(:,MT_bound_cortex(MTBoundOni)).*(vel(:,i)));
        fi_cortex=f0_dynein*(1-vd/v0_dynein); %EQ 2 does not appear to be multiplied by angle
        F_i_cortex(1,i)=F_i_cortex(1,i)+sum(-fi_cortex(1).*(mt_vec(1,MT_bound_cortex(MTBoundOni))).*exp(-MT_length(MT_bound_cortex(MTBoundOni))./(LengthFac*mindistcortex)));%MT_length(MTBoundOni))); %EQ 3
        F_i_cortex(2,i)=F_i_cortex(2,i)+sum(-fi_cortex(1).*(mt_vec(2,MT_bound_cortex(MTBoundOni))).*exp(-MT_length(MT_bound_cortex(MTBoundOni))./(LengthFac*mindistcortex)));%MT_length(MTBoundOni)));
        %F_i_cortex(2,i)=F_i_cortex(2,i)+sum(-fi_cortex(2).*(mt_vec(2,MT_bound_cortex(MTBoundOni))).*exp(-MT_length(MT_bound_cortex(MTBoundOni))./(LengthFac*mindistcortex)));%MT_length(MTBoundOni)));
        clear MTBoundOni
    end
end  
end

