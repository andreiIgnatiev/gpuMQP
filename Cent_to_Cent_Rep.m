function[F_repulsive_centc]=Cent_to_Cent_Rep(nc,centers,C,repd);

F_repulsive_centc=zeros(2,nc); 
for i=1:nc
     for j=(i+1):nc
         distCtoC=sqrt((centers(1,j)-centers(1,i)).^2+(centers(2,j)-centers(2,i)).^2);
         if(distCtoC<repd)
            CtoCvec=[(centers(1,j)-centers(1,i))/distCtoC; (centers(2,j)-centers(2,i))/distCtoC];
            F_repulsive_centc(1,i)=F_repulsive_centc(1,i)+(CtoCvec(1)*C)/(1+distCtoC); %letting the constant "C" be 1 in the numerator
            F_repulsive_centc(2,i)=F_repulsive_centc(2,i)+(CtoCvec(2)*C)/(1+distCtoC);
            F_repulsive_centc(1,j)=F_repulsive_centc(1,j)-(CtoCvec(1)*C)/(1+distCtoC); %letting the constant "C" be 1 in the numerator
            F_repulsive_centc(2,j)=F_repulsive_centc(2,j)-(CtoCvec(2)*C)/(1+distCtoC);
         end
     end
end
