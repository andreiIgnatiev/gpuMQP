function [F_cent_chrB F_chr_chrB cross_bindchr]=ForceChrBoundCalculation(c_attach,ChrPrime,n_MT, mt_centers, ChrBlobPts, LengthFac, X_chr, X_chr_old, centers, centers_old, ChrCentMTInteractions,n_chr,L_chr,nc,MT_length,X_MT,mt_vec,dt,f0_chk,v0_chk);

F_cent_chrB=zeros(2,nc);
F_chr_chrB=zeros(2,n_chr);
cross_bindchr=zeros(1,n_chr);

for i=1:n_chr 
      vchk(:,i)=[(X_chr(1,i)-X_chr_old(1,i));(X_chr(2,i)-X_chr_old(2,i))]/dt;
end

 for j=1:nc 
     vc(:,j)=[(centers(1,j)-centers_old(1,j));(centers(2,j)-centers_old(2,j))]/dt;
 end

%ChrCentMTInteractions=2xn_MT, 1st row=chr #, 2nd row centrosome #. 0 if no interaction, # otherwise

%0.005 binding rate

for i=1:n_chr 
    MTsoni=and((ChrCentMTInteractions(2,:)==1),(mod(ChrCentMTInteractions(1,:),ChrPrime(i))==0));
    for j=(1:nc)
        distChrtoC=sqrt((centers(1,j)-X_chr(1,i)).^2+(centers(2,j)-X_chr(2,i)).^2);
        CtoChrVec=[(X_chr(1,i)-centers(1,j)); (X_chr(2,i)-centers(2,j))]/distChrtoC;
        vidot=real(acos(dot(CtoChrVec,vc(:,j))));
        MTsonj=and(MTsoni,c_attach==j);
%BINDING forces on each centrosome
        Fi_c=f0_chk*(1-(vidot/v0_chk)); 
        x_bind=Fi_c.*(mt_vec(1,MTsonj)).*(exp(-distChrtoC./(LengthFac.*distChrtoC)));
        y_bind=Fi_c.*(mt_vec(2,MTsonj)).*(exp(-distChrtoC./(LengthFac.*distChrtoC)));
        F_cent_chrB(1,j)=F_cent_chrB(1,j)+sum(-x_bind);
        F_cent_chrB(2,j)=F_cent_chrB(2,j)+sum(-y_bind);
        F_chr_chrB(1,i)=F_chr_chrB(1,i)+sum(x_bind); 
        F_chr_chrB(2,i)=F_chr_chrB(1,i)+sum(y_bind);

        ri_vec=[(X_MT(1,MTsonj)-X_chr(1,i)); (X_MT(2,MTsonj)-X_chr(2,i))];
        ri_vec_u=vecnorm(ri_vec,2,1);
        ri_vec=ri_vec./ri_vec_u;
        ri_vec(3,:)=0;
        cross_bindchr(i)=cross_bindchr(i)+sum(sum(cross(ri_vec,[x_bind;y_bind;zeros(1,length(x_bind))],1)));
    end


    
end


