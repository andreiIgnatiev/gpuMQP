function[mt_vec X_MT Fi Fk]=InterpolarForces(LengthFac,intAngle,overlap,CrossLink,nc,fi,fk,FMToni,FMTonk,i,k,mt_vec,X_MT,L_i,L_k,mindistCtoC,MT_length,f0);

Fi=zeros(2,nc);
Fk=zeros(2,nc);


Fx_i=zeros(1,nc);
Fy_i=zeros(1,nc);
Fx_k=zeros(1,nc);
Fy_k=zeros(1,nc);
       f1=find(intAngle>90*pi/180 & intAngle<=120*pi/180);
       f2=find(intAngle>120*pi/180 & intAngle<=180*pi/180);
       if(isempty(f1)==0)
            vdi=f0-(mt_vec(1,FMToni(f1)).*(fi(f1,1)'))+(mt_vec(2,FMToni(f1)).*(fi(f1,2)'));
            vdk=f0-(mt_vec(1,FMTonk(f1)).*(fk(f1,1)'))+(mt_vec(2,FMTonk(f1)).*(fi(f1,2)'));
            Fx_i(1,i)=sum(-vdi.*mt_vec(1,FMToni(f1)).*(1+overlap(f1)).*CrossLink.*exp(-L_i(f1)./((LengthFac)*mindistCtoC)));%*CrossLink);% all part of EQ 6,7
            Fy_i(1,i)=sum(-vdi.*mt_vec(2,FMToni(f1)).*(1+overlap(f1)).*CrossLink.*exp(-L_i(f1)./((LengthFac)*mindistCtoC)));%*CrossLink);
            Fx_k(1,k)=sum(-vdk.*mt_vec(1,FMTonk(f1)).*(1+overlap(f1)).*CrossLink.*exp(-L_k(f1)./((LengthFac)*mindistCtoC)));%*CrossLink);
            Fy_k(1,k)=sum(-vdk.*mt_vec(2,FMTonk(f1)).*(1+overlap(f1)).*CrossLink.*exp(-L_k(f1)./((LengthFac)*mindistCtoC)));%*CrossLink);
       end
       if(isempty(f2)==0)
            vdi=f0-(mt_vec(1,FMToni(f2)).*(fi(f2,1)'))+(mt_vec(2,FMToni(f2)).*(fi(f2,2)'));
            vdk=f0-(mt_vec(1,FMTonk(f2)).*(fk(f2,1)'))+(mt_vec(2,FMTonk(f2)).*(fi(f2,2)'));
            Fx_i(1,i)=sum(-vdi.*mt_vec(1,FMToni(f2)).*(2.*(1+overlap(f2))).*CrossLink.*exp(-L_i(f2)./((LengthFac)*mindistCtoC)));%*CrossLink);
            Fy_i(1,i)=sum(-vdi.*mt_vec(2,FMToni(f2)).*(2.*(1+overlap(f2))).*CrossLink.*exp(-L_i(f2)./((LengthFac)*mindistCtoC)));%*CrossLink);
            Fx_k(1,k)=sum(-vdk.*mt_vec(1,FMTonk(f2)).*(2.*(1+overlap(f2))).*CrossLink.*exp(-L_k(f2)./((LengthFac)*mindistCtoC)));%*CrossLink);
            Fy_k(1,k)=sum(-vdk.*mt_vec(2,FMTonk(f2)).*(2.*(1+overlap(f2))).*CrossLink.*exp(-L_k(f2)./((LengthFac)*mindistCtoC)));%*CrossLink);
       end
Fi=Fi+[Fx_i; Fy_i];
Fk=Fk+[Fx_k; Fy_k];

