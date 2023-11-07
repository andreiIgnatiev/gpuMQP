
%Chromosome to cortex repulsion
function[F_repulsive_chrcortex,mindistchrcortex]=Chr_to_Cortex_Rep(xRadius,n_chr,C,L_chr,x,y,X_chr);


F_repulsive_chrcortex=zeros(2,n_chr);

for i=1:n_chr
        mindistchrcortex=xRadius-sqrt(X_chr(1,i).^2+X_chr(2,i).^2);
        if(mindistchrcortex<2*L_chr)
            [mindistChrtoCortex II]=min(sqrt((x-X_chr(1,i)).^2+(y-X_chr(2,i)).^2));
            ChrtoCortexvec=[(x(II(1))-X_chr(1,i))/mindistChrtoCortex; (y(II(1))-X_chr(2,i))/mindistChrtoCortex];
            F_repulsive_chrcortex(1,i)=F_repulsive_chrcortex(1,i)+(ChrtoCortexvec(1)*C)/(1+mindistchrcortex);%
            F_repulsive_chrcortex(2,i)=F_repulsive_chrcortex(2,i)+(ChrtoCortexvec(2)*C)/(1+mindistchrcortex);% Blows up if distance is really small (close to zero)
        end
end
