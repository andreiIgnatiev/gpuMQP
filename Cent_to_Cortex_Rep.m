function[F_repulsive_centcortex]=Cent_to_Cortex_Rep(xRadius,nc,C,repd,x,y,centers);

F_repulsive_centcortex=zeros(2,nc);
for i=1:nc
        mindistcortex=xRadius-sqrt(centers(1,i).^2+centers(2,i).^2);
        if(mindistcortex<repd)
            [mindistCtoCortex II]=min(sqrt((x-centers(1,i)).^2+(y-centers(2,i)).^2));
            CtoCortexvec=[(x(II(1))-centers(1,i))/mindistCtoCortex; (y(II(1))-centers(2,i))/mindistCtoCortex];
            F_repulsive_centcortex(1,i)=F_repulsive_centcortex(1,i)+(CtoCortexvec(1)*C)/(1+mindistcortex);%
            F_repulsive_centcortex(2,i)=F_repulsive_centcortex(2,i)+(CtoCortexvec(2)*C)/(1+mindistcortex);% Blows up if distance is really small (close to zero)
        end
end