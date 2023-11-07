function[F_repulsive_chr]=Chr_to_Chr_Rep(n_chr,X_chr,C,L_chr);

F_repulsive_chr=zeros(2,n_chr); 
for i=1:n_chr
     for j=(i+1):n_chr
         distChrtoChr=sqrt((X_chr(1,j)-X_chr(1,i)).^2+(X_chr(2,j)-X_chr(2,i)).^2);
         if(distChrtoChr<2*L_chr)
            ChrtoChrvec=[(X_chr(1,j)-X_chr(1,i))/distChrtoChr; (X_chr(2,j)-X_chr(2,i))/distChrtoChr];
            F_repulsive_chr(1,i)=F_repulsive_chr(1,i)+(ChrtoChrvec(1)*C)/(1+distChrtoChr); %letting the constant "C" be 1 in the numerator
            F_repulsive_chr(2,i)=F_repulsive_chr(2,i)+(ChrtoChrvec(2)*C)/(1+distChrtoChr);
            F_repulsive_chr(1,j)=F_repulsive_chr(1,j)-(ChrtoChrvec(1)*C)/(1+distChrtoChr); %letting the constant "C" be 1 in the numerator
            F_repulsive_chr(2,j)=F_repulsive_chr(2,j)-(ChrtoChrvec(2)*C)/(1+distChrtoChr);
         end
     end
end