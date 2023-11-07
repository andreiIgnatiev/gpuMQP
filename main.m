    close all; clear all; clc;
%function MainCluster(param1,param2,param3,param4)
global nc dt
nc=2;% NUMbER OF CENTROSOMES
%nc=str2num(param1);
sims=1;%str2num(param2
%time step
dt=0.5;%str2num(param3);
%dt=str2num(param3);  
fac=0.01;%str2num(param4);
%fac=str2num(param4);
Tfinal=600;
%Tfinal=20*60;
ksims=1;
%while (ksims<=sims)
for(ksims=1)
    %Movie %%%%%%%%%%%%%%%%%%%%%%%%%
%     fid = figure(ksims); 
%     moviename=['nc' num2str(nc) 'sim' num2str(ksims)];
%     writerObj = VideoWriter(moviename,'MPEG-4');
%     writerObj.FrameRate = 4;
%     open(writerObj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %rng('shuffle')
    clear mt_vec X_MT MT_length MT_angle mt_centers MT_state centers
    prefix=['nc' num2str(nc) 'sim' num2str(ksims) 'T'];
    initialization; 
    parameters; %Microtubule_Parameters; 
    centers_old = centers;
    X_chr_old= X_chr;
    net_force=zeros(2,nc);
    clock=0;

    Component_F_stari_cortex_total=zeros(Tfinal,nc);
    Component_F_repulsive_centcortex_total=zeros(Tfinal,nc);
    Component_F_i_cortex_total=zeros(Tfinal,nc);
    Component_F_repulsive_centc_total=zeros(Tfinal,nc);
    Component_F_antiparallel_EG5_total=zeros(Tfinal,nc);
    Component_F_antiparallel_HSET_total=zeros(Tfinal,nc);
    Component_F_DCL_total=zeros(Tfinal,nc);

    for i=1:n_chr
        ChrPrime(i)=nthprime(i);
    end
   % while (clock<=Tfinal)
   for clock = 0:150%200
%         if(clock==0)       
%             filename = [prefix 'Initialization' '.mat'];
%             eval(['save ' filename ' mt_vec MT_angle MT_length c_attach']);
%         end
    clock=clock+1;
        t=dt*clock;
        %%%main time loop %%%%%%%%%%%%
        %CreatenewMTs
        %%
        if(clock>20)
            %%
        if(sum(MT_length)<6000) 
        [c_attachnew mt_centersnew MT_lengthnew MT_anglenew X_MTnew mt_vecnew, MT_statenew] = createnewMTs(mt_nuc, nc, c_rad, centers);
        n_MT=mt_nuc+n_MT; %update total # of MTs
        X_MT = [X_MT X_MTnew];  c_attach = [c_attach c_attachnew];  mt_centers = [mt_centers mt_centersnew];  
        MT_length = [MT_length MT_lengthnew]; MT_angle = [MT_angle MT_anglenew];  mt_vec = [mt_vec mt_vecnew];  
        MT_state = [MT_state MT_statenew]; 
        end
        end
        k2 =fac*gv_astral*dt*MT_length(1,:);
        %k2=0.0125*MT_length(1,:)/2;
        indicesg=find(MT_state==1);
        indicess=find(MT_state==2);
        indicesSlip=find(MT_state==3);
        indicesb=find(MT_state==4);
        Lg=length(indicesg);
        Ls=length(indicess);
        Lslip=length(indicesSlip);
        Lb=length(indicesb);
         
        %State_1
        [MT_state] = State_1(indicesg,Lg,MinDBind,MT_state,k2,dt,xRadius,X_MT,gv_astral,probD); 
 %State_2
        [MT_state] = State_2(indicess,Ls,k1,dt,MT_state);

        %State_3
        [MT_state MT_angle mt_vec X_MT] = State_3(mt_vec,mt_centers,c_rad,MT_length,indicesSlip,Lslip,MinDBind,n_MT,MT_angle,MT_state,k2,dt,xRadius,X_MT,gv_astral,probD,dtheta);

        %State_4    
        [MT_state] = State_4(indicesb,Lb,MinDBind,xRadius,X_MT,MT_state,net_force,mt_vec,f0_dynein,c_attach,gv_astral,dt);
  
        [ChrCentMTInteractions]=ChrCentMT_Calc(nc,ChrPrime,n_chr, n_MT, X_MT, ChrBlobPts, probChr, c_attach, X_chr, centers, r, MT_length, MT_angle);

        [ChrKTCentMTInteractions]=KTCentMT_Calc(nc,ChrPrime,probKT,n_chr, n_MT, X_MT, KTPts, c_attach, X_chr, centers, Krad, MT_length, MT_angle);
        
        %F_i_DCL
          [DCLState F_DCL]=Spindle_Pole_Dynein(c_rad,LengthFac,MinDBind,probD,gv_astral,mt_vec,dt,nc,c_attach,X_MT,centers,centers_old,vf,f0_dynein,v0_dynein,MT_length,MT_angle);
        %%% if its different..
          
          %CtoCortex_repulsive
          [F_repulsive_centcortex]=Cent_to_Cortex_Rep(xRadius,nc,C,repd,x,y,centers);
          
          %CtoC_repulsive
          [F_repulsive_centc]=Cent_to_Cent_Rep(nc,centers,C,repd);
          
          %ChrtoChr_repulsive
          [F_repulsive_chr]=Chr_to_Chr_Rep(n_chr,X_chr,C,L_chr);
          
          %ChrtoCortex_repulsive
          [F_repulsive_chrcortex mindistchrcortex]=Chr_to_Cortex_Rep(xRadius,n_chr,C,L_chr,x,y,X_chr);
          
          %ChrtoCent_repulsive
          [F_repulsive_chrcent F_repulsive_centchr]=Chr_to_Cent_Rep(n_chr,nc, C,L_chr,centers, X_chr);
          
          %F_stari_cortex
          [F_stari_cortex]=Slipping_Force(nc,MT_state,c_attach,f_stall,kappa,MT_length,mt_vec);
%         
%         %F_i_cortex
          [F_i_cortex]=Cortical_Dynein_Force(LengthFac,xRadius,nc,MT_state,c_attach,centers,centers_old,dt,vf,f0_dynein,v0_dynein,mt_vec,MT_length);
          
          %Antiparallel Forces
%           [F_antiparallel_Eg5 F_antiparallel_HSET Eg5State HSETState X_MT mt_vec]=Antiparallel_ForcesV3_old(X_MT,probHSET,probEg5,DCLState,LengthFac,MT_state,incDist,MinDBind,mt_vec,nc,n_MT,MT_length,c_rad,gv_astral,mt_centers,MT_angle,c_attach,centers,centers_old,dt,vf,f0_Eg5,v0_Eg5,f0_HSET,v0_HSET);
           [F_antiparallel_Eg5 F_antiparallel_HSET Eg5State HSETState mt_vec X_MT]=Antiparallel_ForcesV3(LengthFac,vf,clock,CrossLink,drag,X_MT,probHSET,probEg5,MT_state,incDist,MinDBind_Int,mt_vec,nc,MT_length,c_rad,mt_centers,MT_angle,c_attach,centers,centers_old,dt,f0_Eg5,v0_Eg5,f0_HSET,v0_HSET);
           
           [F_cent_chrS F_chr_chrS cross_slipchr]=ForceChrSlipCalculation(c_attach,ChrPrime,ChrCentMTInteractions,X_chr,n_chr,nc,MT_length,X_MT,mt_vec,kappa,f_stall);
        %ForceOnCentDueToMTChrSlip %ForceOnChrDueToMTChrSlip

            [F_cent_chrB F_chr_chrB cross_bindchr]=ForceChrBoundCalculation(c_attach,ChrPrime,n_MT, mt_centers, ChrBlobPts, LengthFac, X_chr, X_chr_old, centers, centers_old, ChrCentMTInteractions,n_chr,L_chr,nc,MT_length,X_MT,mt_vec,dt,f0_chk,v0_chk);
          %ForceOnCentDueToMTChrBind %ForceOnChrDueToMTChrBind
       %[F_cent_chrBKT F_chr_chrBKT cross_bindchrKT]=ForceChrBoundKTCalculation(c_attach,ChrPrime,n_MT, mt_centers, ChrBlobPts, LengthFac, X_chr, X_chr_old, centers, centers_old, ChrKTCentMTInteractions,n_chr,L_chr,nc,MT_length,X_MT,mt_vec,dt,f0_chk,v0_chk);

            centers_old=centers;     
       
              X_chr_old=X_chr;


       
        distCtoC=sqrt((centers(1,1)-centers(1,2)).^2+(centers(2,1)-centers(2,2)).^2);
        CtoCvec=[(centers(1,2)-centers(1,1))/distCtoC; (centers(2,2)-centers(2,1))/distCtoC];
        
        

        for i=1:nc
    AvgLength(i)=mean(MT_length(c_attach==i));
	AvgMTVol(i)=pi*0.0125^2*AvgLength(i);
	TotMTVol(i)=AvgMTVol(i)*length(find(c_attach==i));
	r_s(i)=radii(1)+AvgLength(i);
	SliceVol(i)=pi*r_s(i)^2*0.05;
	VolFrac(i)=TotMTVol(i)/SliceVol(i);
	j(i)=mean(X_MT(1,c_attach==i));
	k(i)=mean(X_MT(2,c_attach==i));
	CoM(:,i)=[j(i); k(i)];
	DistToCenter(i)=sqrt((-CoM(1,i)).^2+(-CoM(2,i)).^2);
	Drag(i)=(6*pi*2.1*r_s(i)*VolFrac(i))/(1-(DistToCenter(i)./xRadius).^2);
            net_force(:,i)=F_stari_cortex(:,i)+F_i_cortex(:,i)+F_repulsive_centcortex(:,i)+F_repulsive_centc(:,i)+F_antiparallel_Eg5(:,i)+F_antiparallel_HSET(:,i)+F_DCL(:,i)+F_repulsive_centchr(:,i)+F_cent_chrS(:,i)+F_cent_chrB(:,i);%+F_cent_chrBKT(:,i)
            centers(:,i) = centers(:,i)-(dt/drag)*net_force(:,i);
            
            Component_F_stari_cortex(i)=dot(CtoCvec,F_stari_cortex(:,i));
            Component_F_stari_cortex_total(clock,i)= Component_F_stari_cortex(i);
            Component_F_i_cortex(i)=dot(CtoCvec,F_i_cortex(:,i));
            Component_F_i_cortex_total(clock,i)= Component_F_i_cortex(i);
            Component_F_repulsive_centcortex(i)=dot(CtoCvec,F_repulsive_centcortex(:,i));
            Component_F_repulsive_centcortex_total(clock,i)= Component_F_repulsive_centcortex(i);
            Component_F_repulsive_centc(i)=dot(CtoCvec,F_repulsive_centc(:,i));
            Component_F_repulsive_centc_total(clock,i)= Component_F_repulsive_centc(i);
            Component_F_antiparallel_Eg5(i)=dot(CtoCvec,F_antiparallel_Eg5(:,i));
            Component_F_antiparallel_EG5_total(clock,i)= Component_F_antiparallel_Eg5(i);
            Component_F_antiparallel_HSET(i)=dot(CtoCvec,F_antiparallel_HSET(:,i));
            Component_F_antiparallel_HSET_total(clock,i)= Component_F_antiparallel_HSET(i);
            Component_F_DCL(i)=dot(CtoCvec,F_DCL(:,i));
            Component_F_DCL_total(clock,i)= Component_F_DCL(i);
        end
        
        
        for i=1:n_chr
            net_force_chr(:,i)= F_repulsive_chr(:,i)+ F_repulsive_chrcortex(:,i)+F_repulsive_chrcent(:,i)+F_chr_chrS(:,i)+F_chr_chrB(:,i);%++F_chr_chrBKT(:,i)
            X_chr(:,i) = X_chr(:,i)-(dt/drag)*net_force_chr(:,i);
            alpha(i)=alpha(i)-(dt/10)*(cross_slipchr(i)+cross_bindchr(i));
            ChrBlobPts(1,:,i)=r.*cos(theta+alpha(i))+X_chr(1,i);
            ChrBlobPts(2,:,i)=r.*sin(theta+alpha(i))+X_chr(2,i);
            KTPts(1,:,i)=Krad.*cos(delta+alpha(i))+X_chr(1,i);
            KTPts(2,:,i)=Krad.*sin(delta+alpha(i))+X_chr(2,i);
        end
        
        for i=1:n_chr
            [mindistChrtoCortex II]=min(sqrt((x-X_chr(1,i)).^2+(y-X_chr(2,i)).^2));
            DistChrtoCortex(i,clock)=mindistChrtoCortex;
        end
        
        for i=1:nc
            [mindistCtoCortex II]=min(sqrt((x-centers(1,i)).^2+(y-centers(2,i)).^2));
            DistCtoCortex(i,clock)=mindistCtoCortex;
        end

        %Center of centrosome MT is attached to
        for i = 1:n_MT
            mt_centers(:,i) = [centers(:,c_attach(i))];  
        end
        
      %%== plotting MTs (different color for each state; L1, L2...)  
      if(mod(t,1)==0)
%         figure(5)
%        % movegui(figure(5),[600 300])
%         for i=1:n_chr
%             plot(DistChrtoCortex(i,:),'LineWidth',2)
%             hold on
%         end
%             title(['Chromosome Distances from the Cortex'])
%         xlabel('Time')
%         ylabel('Distance')
%         hold off
% 
%         figure (6)
%        % movegui(figure(6),[600 0])
%         for i=1:nc
%             plot(DistCtoCortex(i,:),'LineWidth',3)
%             hold on 
%         end 
%         title(['Centrosome Distance from the Cortex'])
%         xlabel('Time')
%         ylabel('Distance')
%         hold off
%        filename = [prefix num2str(t) '.mat'];
%         eval(['save ' filename ' MT_state mt_centers centers c_attach n_MT MT_length X_MT mt_vec MT_angle Eg5State DCLState HSETState ComboState F_stari_cortex F_i_cortex F_antiparallel_Eg5 F_antiparallel_HSET F_DCL F_repulsive_centc F_repulsive_centcortex']);
        indices1=find(MT_state==1);
        indices2=find(MT_state==2);
        indices3=find(MT_state==3);
        indices4=find(MT_state==4);
        L1=length(indices1);
        L2=length(indices2);
        L3=length(indices3);
        L4=length(indices4);
        [L5 col5]=size(Eg5State);
        [L6 col6]=size(HSETState);
        L7=length(DCLState);
        
  %for i=1:nc       
   %      figure(i)
%MTVec=[Component_F_i_cortex(i); Component_F_stari_cortex(i); Component_F_antiparallel_Eg5(i); Component_F_antiparallel_HSET(i); Component_F_DCL(i)];
%MTVec=MTVec';
%if(i==2)
 %   MTVec=MTVec*(-1);
%end
%MTVecNeg=MTVec;
%MTVecNeg(MTVecNeg>0)=0;
%MTVecPos=MTVec;
%MTVecPos(MTVecPos<0)=0;
%H = bar(MTVecNeg,'stacked');
%hold on
 %G = bar(MTVecPos,'stacked');
 %ylim([-5 5])
%hold off
% colorSet = [];
% H(1).FaceColor='k';
% H(2).FaceColor='c';
% %H(3).FaceColor='b';
% H(3).FaceColor='r';
%   H(4).FaceColor='y';

 %     end


 figure(3)
for i=1:n_chr
    plot(KTPts(1,:,i),KTPts(2,:,i),'r','LineWidth',1)
    hold on
    plot(ChrBlobPts(1,:,i),ChrBlobPts(2,:,i),'r','LineWidth',1)
end


           if(L1>0)
           for i=1:L1
               plot([mt_centers(1,indices1(i))+c_rad*cos(MT_angle(indices1(i))) X_MT(1,indices1(i))], [mt_centers(2,indices1(i))+c_rad*sin(MT_angle(indices1(i))) X_MT(2,indices1(i))],'b','linewidth',0.025)
               hold on
           end 
           end
       if(L2>0)
           for i=1:L2
               plot([mt_centers(1,indices2(i))+c_rad*cos(MT_angle(indices2(i))) X_MT(1,indices2(i))], [mt_centers(2,indices2(i))+c_rad*sin(MT_angle(indices2(i))) X_MT(2,indices2(i))],'--b','linewidth',0.025)
               hold on
           end
       end
       if(L3>0)
            for i=1:L3
                plot([mt_centers(1,indices3(i))+c_rad*cos(MT_angle(indices3(i))) X_MT(1,indices3(i))], [mt_centers(2,indices3(i))+c_rad*sin(MT_angle(indices3(i))) X_MT(2,indices3(i))],'c','linewidth',0.025)
                hold on
            end
        end
        if(L4>0)
            for i=1:L4
                plot([mt_centers(1,indices4(i))+c_rad*cos(MT_angle(indices4(i))) X_MT(1,indices4(i))], [mt_centers(2,indices4(i))+c_rad*sin(MT_angle(indices4(i))) X_MT(2,indices4(i))],'k','linewidth',0.025)
                hold on
            end
        end
        if(L5>0)
            for i=1:L5
                plot([mt_centers(1,Eg5State(i,1))+c_rad*cos(MT_angle(Eg5State(i,1))) X_MT(1,Eg5State(i,1))], [mt_centers(2,Eg5State(i,1))+c_rad*sin(MT_angle(Eg5State(i,1))) X_MT(2,Eg5State(i,1))],'r','linewidth',0.025)
                plot([mt_centers(1,Eg5State(i,2))+c_rad*cos(MT_angle(Eg5State(i,2))) X_MT(1,Eg5State(i,2))], [mt_centers(2,Eg5State(i,2))+c_rad*sin(MT_angle(Eg5State(i,2))) X_MT(2,Eg5State(i,2))],'r','linewidth',0.025)
                hold on
            end
        end
        if(L6>0)
            for i=1:L6
                 plot([mt_centers(1,HSETState(i,1))+c_rad*cos(MT_angle(HSETState(i,1))) X_MT(1,HSETState(i,1))], [mt_centers(2,HSETState(i,1))+c_rad*sin(MT_angle(HSETState(i,1))) X_MT(2,HSETState(i,1))],'--m','linewidth',0.025)
                 plot([mt_centers(1,HSETState(i,2))+c_rad*cos(MT_angle(HSETState(i,2))) X_MT(1,HSETState(i,2))], [mt_centers(2,HSETState(i,2))+c_rad*sin(MT_angle(HSETState(i,2))) X_MT(2,HSETState(i,2))],'--m','linewidth',0.025)
                 hold on
            end
        end
         if(L7>0)
            for i=1:L7
                 plot([mt_centers(1,DCLState(i))+c_rad*cos(MT_angle(DCLState(i))) X_MT(1,DCLState(i))], [mt_centers(2,DCLState(i))+c_rad*sin(MT_angle(DCLState(i))) X_MT(2,DCLState(i))],'g','linewidth',0.025)
                hold on
            end
         end
      %  for i=1:length(ChrCentMTInteractions(5,:)==1)
      %     plot([mt_centers(1,(ChrCentMTInteractions(5,:)==1))+c_rad*cos(MT_angle(ChrCentMTInteractions(5,:)==1)) X_MT(1,(ChrCentMTInteractions(5,:)==1))], [mt_centers(2,(ChrCentMTInteractions(5,:)==1))+c_rad*sin(MT_angle(ChrCentMTInteractions(5,:)==1)) X_MT(2,(ChrCentMTInteractions(5,:)==1))],'y','linewidth',0.025)
      %          hold on
      %  end
         %Bound to Chr
         
        plot(X_MT(1,ChrCentMTInteractions(2,:)==1),X_MT(2,ChrCentMTInteractions(2,:)==1),'mo','MarkerSize',3)
        % Slipping on Chr
        plot(X_MT(1,ChrCentMTInteractions(3,:)==1),X_MT(2,ChrCentMTInteractions(3,:)==1),'go','MarkerSize',3)
        plot(X_MT(1,ChrKTCentMTInteractions(2,:)==1),X_MT(2,ChrKTCentMTInteractions(2,:)==1),'co','MarkerSize',5)
        for i = 1:nc
           h=viscircles(centers(:,i)',radii(i));
           hold on
           colors={'y','y','m','c'};
           h.Children(1).Color = colors{i};         
           h.Children(2).Color = 'k';
       end

        axis equal 
        axis([xmin xmax ymin ymax])
        axis equal 
        box off
        plot(x,y,'-k', 'linewidth',5) 
        hold on

        plottitle=['t= ' num2str(clock)]; 
        set(gcf,'color','white')
        title(plottitle,'FontWeight','bold','FontSize',18)
% frame = getframe(gcf);
%         writeVideo(writerObj, frame)

%3 is binding, 4 is slipping 
checkChrSlip=length(find(ChrCentMTInteractions(2,:)~=0));
checkChrB=length(find(ChrCentMTInteractions(3,:)~=0));

hold off

% figure(7) % plot of binding / slipping chr interactions, works
%    set(gcf,'color','white')
%    set(gcf,'defaultaxesfontsize',18)
%    chr_x1=categorical({'ChrSlip','ChrBound'})
%    chr_y1=[checkChrSlip;checkChrB];
%    b=bar(chr_x1,chr_y1)
%    b.FaceColor = 'flat';
%    b.CData(1,:) = [.5 0 .5 ];
%    b.CData(2,:) = [.5 0 .3 ];
%     
%    title(['MT Chr interactions at time =' num2str(clock)])
%    xlabel('MT Chr Interaction')
%    ylabel('# of Interactions')


  pause(.0001)
  hold off
%   pause
      end  
    clear k2
%     Eg5State=[];
%     HSETState=[];
%     ComboState=[];
    if (clock>15)
    [MT_length] = MT_update(MT_state,MT_length,dt,gv_astral,sv_astral,bsv_astral);
    end
    %deleteMTs  
    [c_attach mt_centers MT_length MT_angle X_MT mt_vec MT_state n_MT] = delete_MTs(n_MT, MT_length, c_attach, mt_centers, MT_angle, X_MT, mt_vec, MT_state, c_rad);
    DistToCenter=sqrt(centers(1,:).^2+centers(2,:).^2);
    CheckOut=find(DistToCenter>yRadius);
    if(isempty(CheckOut)==0)
        clock=Tfinal+1
        display('exiting')
        
    end
    end
%     if(isempty(CheckOut)==0)
%         display('deleting')
%         filemod=[prefix '*.mat']
%         delete(filemod)
%         
%     else
%         ksims=ksims+1;
%     end
%     close(writerObj); 
%figure (77) 


end