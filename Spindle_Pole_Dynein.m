function[DCLState F_DCL]=Spindle_Pole_Dynein(c_rad,LengthFac,MinDBind,probD,gv_astral,mt_vec,dt,nc,c_attach,X_MT,centers,centers_old,vf,f0_dynein,v0_dynein,MT_length,MT_angle);

DCLState=[];
vec=[1:1:nc];
F_DCL=zeros(2,nc);
for i=1:nc
%    magvel(1,i)=sqrt((centers(1,i)-centers_old(1,i)).^2+(centers(2,i)-centers_old(2,i)).^2)/dt;
vel(:,i)=[(centers(1,i)-centers_old(1,i));(centers(2,i)-centers_old(2,i))]/dt;
end
%vel=[magvel; magvel];
for i=1:nc
     veci=vec;
     veci(i)=[];
     for jj=1:nc-1
         j=veci(jj);
         [Indexi Indexj CtoCVecItoj CtoCVecjtoI]=InterpolarMTs(i,j,mt_vec,centers,c_attach,MT_angle); 
         for k=1:length(Indexj)
         %vector fromf j to i
         mindistCtoC=sqrt((centers(1,i)-centers(1,j)).^2+(centers(2,i)-centers(2,j)).^2);
         if(MT_length(Indexj(k))>=(mindistCtoC-MinDBind)) 
         CtoCvec=[(centers(1,j)-centers(1,i))/mindistCtoC; (centers(2,j)-centers(2,i))/mindistCtoC];
         slope=(X_MT(2,Indexj(k))-centers(2,j))./(X_MT(1,Indexj(k))-centers(1,j));
         a=-slope;
         b=1;
         c=-(X_MT(2,Indexj(k))-slope.*X_MT(1,Indexj(k)));
         xclose=(b.*(b*centers(1,i)-a*centers(2,i))-a.*c)./(a.^2+b.^2);
         yclose=(a.*(-b*centers(1,i)+a*centers(2,i))-b.*c)./(a.^2+b.^2);
         distMTtoi=sqrt((xclose-centers(1,i)).^2+(yclose-centers(2,i)).^2);
         MT_DCL=find(distMTtoi<1);
          n=rand(1,k);
        if(isempty(MT_DCL)==0)%0 if not an empty vector
            if(n(k)<0.1);%probD)
                vd=sum((vel(:,i)-vel(:,j)-vf)).*mt_vec(:,Indexj(k));
                fi_DCL=f0_dynein*(1-vd/v0_dynein);
                DCLState=[DCLState Indexj(k)];
                F_DCL(1,j)=F_DCL(1,j)+sum(-fi_DCL(1).*(mt_vec(1,Indexj(k))).*exp(-mindistCtoC./(LengthFac*mindistCtoC)));
                F_DCL(2,j)=F_DCL(2,j)+sum(-fi_DCL(2).*(mt_vec(2,Indexj(k))).*exp(-mindistCtoC./(LengthFac*mindistCtoC)));
          
                F_DCL(1,i)=F_DCL(1,i)+sum(fi_DCL(1).*(mt_vec(1,Indexj(k))).*exp(-mindistCtoC./(LengthFac*mindistCtoC)));
               F_DCL(2,i)=F_DCL(2,i)+sum(fi_DCL(2).*(mt_vec(2,Indexj(k))).*exp(-mindistCtoC./(LengthFac*mindistCtoC)));
  end
        end
         end
         end
     end
end
DCLState=unique(DCLState);

