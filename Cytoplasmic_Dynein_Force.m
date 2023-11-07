function[F_c]=Cytoplasmic_Dynein_Force(eta,X_MT,LengthFac,xRadius,nc,MT_state,c_attach,centers,centers_old,dt,vf,f0_dynein,v0_dynein,mt_vec,MT_length);

F_c=zeros(2,nc);

for i=1:nc
    vel(:,i)=[(centers(1,i)-centers_old(1,i));(centers(2,i)-centers_old(2,i))]/dt;
%      magvel(1,i)=sqrt((centers(1,i)-centers_old(1,i)).^2+(centers(2,i)-centers_old(2,i)).^2)/dt;
end
%  vel=[magvel; magvel];
MT_stateNew=MT_state;
MT_stateNew(find(MT_state==2))=[];
MT_stateNew=find(MT_state(MT_stateNew));
 for i=1:nc
     MTsoni=find(c_attach==i);
%      BoundMTsoni=intersect(MTsoni,MT_stateNew);
     for j=1:length(MTsoni)
     %for j=1:length(BoundMTs)
     vd=-vel(:,i)-vf;
     fc=f0_dynein*(1-vd/v0_dynein);
     F_c(1,i)=F_c(1,i)+sum(-eta.*MT_length(MTsoni(j)).*mt_vec(1,MTsoni(j)));
     F_c(2,i)=F_c(2,i)+sum(-eta.*MT_length(MTsoni(j)).*mt_vec(2,MTsoni(j)));
     end
 end
%  F_c
