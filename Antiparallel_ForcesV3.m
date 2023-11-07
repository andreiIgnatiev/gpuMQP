function[F_antiparallel_Eg5 F_antiparallel_HSET Eg5StateFinal HSETStateFinal mt_vec X_MT]=Antiparallel_ForcesV3(LengthFac,vf,clock,CrossLink,drag,X_MT,probHSET,probEg5,MT_state,incDist,MinDBind_Int,mt_vec,nc,MT_length,c_rad,mt_centers,MT_angle,c_attach,centers,centers_old,dt,f0_Eg5,v0_Eg5,f0_HSET,v0_HSET)

F_antiparallel_Eg5=zeros(2,nc);
Fi_Eg5=zeros(2,nc);
Fk_Eg5=zeros(2,nc);
F_antiparallel_HSET=zeros(2,nc);
Fi_HSET=zeros(2,nc);
Fk_HSET=zeros(2,nc);
Eg5StateFinal=[];
HSETStateFinal=[];
indicess=find(MT_state==2);
pbinde=0;
pbindh=0;

for i=1:nc
%      magvel(1,i)=sqrt((centers(1,i)-centers_old(1,i)).^2+(centers(2,i)-centers_old(2,i)).^2)/dt;
   vel(:,i)=[(centers(1,i)-centers_old(1,i));(centers(2,i)-centers_old(2,i))]/dt;
    
end

%  vel=[magvel; magvel];
%1 to 2
%1 to 2, 1 to 3, 2 to 3
% if(vel(:,:)<=v0_Eg5)
for i = 1:nc
    vec=[1:1:nc];
    vec(i)=[];
    MTsoni=find(c_attach==i);
   for kk=1:length(vec)
        k=vec(kk); 
        MTsonk=find(c_attach==k);
        maxMTsoni_length=max(MT_length(MTsoni));
        maxMTsonk_length=max(MT_length(MTsonk));
        mindistCtoC=sqrt((centers(1,k)-centers(1,i)).^2+(centers(2,k)-centers(2,i)).^2);
        
        if (maxMTsoni_length+maxMTsonk_length)>mindistCtoC
            [Indexi Indexk CtoCVecItoj CtoCVecjtoI]=InterpolarMTs(i,k,mt_vec,centers,c_attach,MT_angle);
            PtMTsoni=[];
            PtMTsonk=[];
            vector=[];
            for n=1:length(Indexi)
                increments=floor((MT_length(Indexi(n))+c_rad)/(incDist))+1;
                incMTsoni(n)=increments;
                for j=1:increments
                    ind=sum(incMTsoni(1:n-1))+j;
                    PtMTsoni(:,ind)=[mt_centers(1,Indexi(n))+(incDist*(j-1))*cos(MT_angle(Indexi(n)));mt_centers(2,Indexi(n))+(incDist*(j-1))*sin(MT_angle(Indexi(n)))];
                end
            end
            for n=1:length(Indexk)
                increments=floor((MT_length(Indexk(n))+c_rad)/(incDist))+1;
                incMTsonk(n)=increments;
                for j=1:increments
                    ind=sum(incMTsonk(1:n-1))+j;
                    PtMTsonk(:,ind)=[mt_centers(1,Indexk(n))+(incDist*(j-1))*cos(MT_angle(Indexk(n)));mt_centers(2,Indexk(n))+(incDist*(j-1))*sin(MT_angle(Indexk(n)))];
                    vector(ind)=Indexk(n); 
                end
            end
            L_iHSET=[];
            L_iEg5=[];
            L_kHSET=[];
            L_kEg5=[];
            ForceDistToCenterEg5=[];
            ForceDistToCenterHSET=[];
            Eg5State=[];
            HSETState=[];
            BoundMTsonk=[];
            for j=1:length(Indexi)
                MinDSave=10;%dumby variable that is updated later
                for qq=1:incMTsoni(j)
                    ind1=sum(incMTsoni(1:(j-1)))+qq;
                    distMTtoMT=sqrt((PtMTsonk(1,:)-PtMTsoni(1,ind1)).^2+(PtMTsonk(2,:)-PtMTsoni(2,ind1)).^2);
                    [MinD boundIndices]=min(distMTtoMT);
                    if(MinD<MinDSave)
                        MinDSave=MinD(1); %Distance between
                        boundIndex=boundIndices(1);
                        BoundMTs=vector(boundIndices);
                        BoundMTsonk=[BoundMTsonk BoundMTs];
                        BoundMTsonk=unique(BoundMTsonk);
                        MTFirstNumber=vector(boundIndex);%This is the MT on k that is bound to the closest MT on i
                        q=qq;
                    end
                end
                pbinde=0;
                pbindh=0;
%                 if(clock>=100)
%                     probEg5=0;
%                 end
                if(MinDSave<MinDBind_Int)
                        RVeg5=rand(1,1);
                        RVhset=rand(1,1);
                        pbinde(RVeg5<probEg5)=1;
                        pbindh(RVhset<probHSET)=1;
%                         if(pbindh==1)
%                             if(pbinde==1)
%                                 undo=rand(1,1);
%                                 if(undo<0.5)
%                                     pbinde=0;
%                                 elseif(undo>=0.5)
%                                     pbindh=0;
%                                 end
%                             end
%                         end
                        if(pbindh==1)
                            HSETState=[HSETState; Indexi(j) MTFirstNumber];
                            [L_i L_k]=motorPt(q,incDist,centers,PtMTsonk,boundIndex,i,k);
                            L_iHSET=[L_iHSET; L_i];
                            L_kHSET=[L_kHSET; L_k];
                        end
                        if(pbinde==1)
                            Eg5State=[Eg5State; Indexi(j) MTFirstNumber];
                            [L_i L_k]=motorPt(q,incDist,centers,PtMTsonk,boundIndex,i,k);
                            L_iEg5=[L_iEg5; L_i];
                            L_kEg5=[L_kEg5; L_k];
                        end
                end%if statement, if close and forces 
            end%loop on j MTs (Indexi(j))
            [rowEg5 colEg5]=size(Eg5State);
            [rowHSET colHSET]=size(HSETState);
            if(rowEg5>1) %Enter loop if there are two rows in Eg5State
                State=Eg5State;
                L_i=L_iEg5;
                L_k=L_kEg5;
                [Eg5State L_iEg5 L_kEg5]=RemoveRepeats(State,L_i,L_k);
            end
            if(rowHSET>1) %Enter loop if there are more than two rows in Eg5State
                State=HSETState;
                L_i=L_iHSET;
                L_k=L_kHSET; 
                [HSETState L_iHSET L_kHSET]=RemoveRepeats(State,L_i,L_k);
            end
            [rowEg5 colEg5]=size(Eg5State);
            [rowHSET colHSET]=size(HSETState);
            

            mag_vel_i(:,1)=vel(:,i);
            mag_vel_k(:,1)=vel(:,k);
            vi=mag_vel_i(:);
            vk=mag_vel_k(:);
            ViDot=acos(dot(CtoCVecItoj,vi));
            VkDot=acos(dot(CtoCVecjtoI,vk));
            if(pi/4<=ViDot<3*pi/4 && pi/4<=VkDot<3*pi/4)
                v=vi+vk;
            elseif((ViDot<=pi/4 || ViDot>=3*pi/4) && (VkDot<=pi/4 || VkDot>=3*pi/4))
                v=vi+vk;
            else
                v=vi-vk;
            end
%             v
%             vi_Eg5=(-vel(:,i)+vel(:,k)-vf);
%             vk_Eg5=(-vel(:,k)+vel(:,i)-vf);
%              vi_HSET=(vel(:,i)-vel(:,k)-vf);
%              vk_HSET=(vel(:,k)-vel(:,i)-vf);
            if(isempty(Eg5State)==0)
                FMToni=Eg5State(:,1);
                FMTonk=Eg5State(:,2);
                L_i=L_iEg5';
                    L_k=L_kEg5';
                for j=1:length(FMToni)
                    dotInt(j)=dot(mt_vec(:,FMToni(j)),mt_vec(:,FMTonk(j)));
                    intAngle(j)=acos(dotInt(j));%remove abs(dotInt)
                    a(j)=MT_length(FMToni(j))-L_i(j);
                    b(j)=MT_length(FMTonk(j))-L_k(j);
                    c(j)=sqrt(a(j)^2+b(j)^2-2*a(j)*b(j)*cos(intAngle(j)));
                    overlapMatrix(:,j)=[MT_length(FMTonk(j));MT_length(FMToni(j));c(j)/MT_length(FMTonk(j));c(j)/MT_length(FMToni(j))];%abs((MT_length(FMTonk)+MT_length(FMToni))-mindistCtoC)]%add abs of the final value of the matrix
                    overlap(j)=min(overlapMatrix(:,j));
                    
%                     v
%                     overlap
%                     overlap(j)
%                     (1-(v/v0_Eg5))
%                     CrossLink*overlap(j)
%                     (1-(vi)/v0_Eg5)
                    fi(j,:)=-f0_Eg5*(1-v/v0_Eg5);
                    fk(j,:)=-f0_Eg5*(1-v/v0_Eg5);
%                     fi(j,:)=-f0_Eg5*max(overlap(j)*CrossLink,(1-(v/v0_Eg5)));
%                     fk(j,:)=-f0_Eg5*max(overlap(j)*CrossLink,(1-(v/v0_Eg5)));
                end
%                 fi=-f0_Eg5*(1-(vi(:)/v0_Eg5));
%                 fk=-f0_Eg5*(1-(vk(:)/v0_Eg5));
                [mt_vec X_MT Fi_Eg5 Fk_Eg5]=InterpolarForces(LengthFac,intAngle,overlap,CrossLink,nc,fi,fk,FMToni,FMTonk,i,k,mt_vec,X_MT,L_i,L_k,mindistCtoC,MT_length,-f0_Eg5);
                clear intAngle overlap FMToni FMTonk c intAngle fi fk L_i L_k
            end

            if(isempty(HSETState)==0)
                FMToni=HSETState(:,1);
                FMTonk=HSETState(:,2);
                 L_i=L_iHSET';
                    L_k=L_kHSET';
                for j=1:length(FMToni)
                   
                    dotInt(j)=dot(mt_vec(:,FMToni(j)),mt_vec(:,FMTonk(j)));
                    intAngle(j)=acos(dotInt(j));%remove abs(dotInt)
                    a(j)=MT_length(FMToni(j))-L_i(j);
                    b(j)=MT_length(FMTonk(j))-L_k(j);
                    c(j)=sqrt(a(j)^2+b(j)^2-2*a(j)*b(j)*cos(intAngle(j)));
                    overlapMatrix(:,j)=[MT_length(FMTonk(j));MT_length(FMToni(j));c(j)/MT_length(FMTonk(j));c(j)/MT_length(FMToni(j))];%abs((MT_length(FMTonk)+MT_length(FMToni))-mindistCtoC)]%add abs of the final value of the matrix
                    overlap(j)=min(overlapMatrix(:,j));
%                     vi(:)/v0_HSET
%                     vk(:)/v0_HSET
%                     fi(j,:)=f0_HSET*max(overlap(j)*CrossLink,(1-(v/v0_HSET)));
%                     fk(j,:)=f0_HSET*max(overlap(j)*CrossLink,(1-(v/v0_HSET)));
                    fi(j,:)=f0_HSET*(1-v/v0_HSET);
                    fk(j,:)=f0_HSET*(1-v/v0_HSET);
%                     fi(j,:)=f0_HSET*max(CrossLink*overlap(j),(1-(abs(vi)/v0_HSET)));
%                     fk(j,:)=f0_HSET*max(CrossLink*overlap(j),(1-(abs(vk)/v0_HSET)));
                end
%                 fi=f0_HSET*(1-(vi(:)/v0_HSET));
%                 fk=f0_HSET*(1-(vk(:)/v0_HSET));
                [mt_vec X_MT Fi_HSET Fk_HSET]=InterpolarForces(LengthFac,intAngle,overlap,CrossLink,nc,fi,fk,FMToni,FMTonk,i,k,mt_vec,X_MT,L_i,L_k,mindistCtoC,MT_length,f0_HSET);
                clear intAngle overlap FMToni FMTonk c intAngle fi fk L_i L_k
            end
            Eg5StateFinal=[Eg5StateFinal; Eg5State];
            HSETStateFinal=[HSETStateFinal; HSETState];
        end%long enough MTs to potentially have forces
        clear Indexi Indexk 
        F_antiparallel_Eg5=F_antiparallel_Eg5+Fi_Eg5+Fk_Eg5;
        F_antiparallel_HSET=F_antiparallel_HSET+Fi_HSET+Fk_HSET;
    end%loop on k centrosome
end%loop on i centrosome
% end

