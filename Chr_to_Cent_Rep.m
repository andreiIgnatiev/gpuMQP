%Chr to Cent Rep

function[F_repulsive_chrcent F_repulsive_centchr]=Chr_to_Cent_Rep(n_chr,nc, C,L_chr,centers, X_chr); 

F_repulsive_chrcent=zeros(2,n_chr);
F_repulsive_centchr=zeros(2,nc);
for i=1:n_chr 
    for j=(1:nc)
       % cent centers and chr centers?
        distChrtoC=sqrt((centers(1,j)-X_chr(1,i)).^2+(centers(2,j)-X_chr(2,i)).^2);
        if(distChrtoC<2*L_chr)
            ChrtoCvec=[(centers(1,j)-X_chr(1,i))/distChrtoC; (centers(2,j)-X_chr(2,i))/distChrtoC];
            F_repulsive_chrcent(1,i)=F_repulsive_chrcent(1,i)+(ChrtoCvec(1)*C)/(1+distChrtoC); %letting the constant "C" be 1 in the numerator
            F_repulsive_chrcent(2,i)=F_repulsive_chrcent(2,i)+(ChrtoCvec(2)*C)/(1+distChrtoC);
            F_repulsive_centchr(1,j)=F_repulsive_centchr(1,j)-(ChrtoCvec(1)*C)/(1+distChrtoC); %letting the constant "C" be 1 in the numerator
            F_repulsive_centchr(2,j)=F_repulsive_centchr(2,j)-(ChrtoCvec(2)*C)/(1+distChrtoC);
        end
        %L_chr & C 
    end
end