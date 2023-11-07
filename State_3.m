%State 3; slipping on cortex

function[MT_state MT_angle mt_vec X_MT] = State_3(mt_vec,mt_centers,c_rad,MT_length,indicesSlip,Lslip,MinDBind,n_MT,MT_angle,MT_state,k2,dt,xRadius,X_MT,gv_astral,probD,dtheta);

% indicesSlip = find(MT_state==3); %pulls out MTs that are slipping on cortex
% Lslip = length(indicesSlip);
betasign=randi([-1,1],n_MT);
betasign(betasign==0)=[];
%for MTs slipping on cortex
if (Lslip>0)
    for i=1:Lslip
        nb=rand(1,1);
         mindistcortex=xRadius-sqrt((X_MT(1,indicesSlip(i))).^2+(X_MT(2,indicesSlip(i))).^2);
         if mindistcortex<0
            MT_state(indicesSlip(i))=2;
         elseif((mindistcortex<MinDBind) & (nb<probD))
            %Might bind
            MT_state(indicesSlip(i)) = 4;  
         elseif(mindistcortex>=MinDBind)
             MT_state(indicesSlip(i))=1;
         else
            %beta is from pavin et al New Journal of Physics 2012
            beta=betasign(indicesSlip(i)); %~pi/100
            MT_angle(indicesSlip(i))=MT_angle(indicesSlip(i))-(beta*dtheta);
       	    X_MT(:,indicesSlip(i))=[(MT_length(indicesSlip(i))+c_rad).*cos(MT_angle(indicesSlip(i)))+mt_centers(1,indicesSlip(i)); (MT_length(indicesSlip(i))+c_rad).*sin(MT_angle(indicesSlip(i)))+mt_centers(2,indicesSlip(i))];
	    mt_vec(1,indicesSlip(i))=(1./(MT_length(indicesSlip(i)))).*(X_MT(1,indicesSlip(i))-(mt_centers(1,indicesSlip(i))+c_rad*cos(MT_angle(indicesSlip(i)))));
	    mt_vec(2,indicesSlip(i))=(1./(MT_length(indicesSlip(i)))).*(X_MT(2,indicesSlip(i))-(mt_centers(2,indicesSlip(i))+c_rad*cos(MT_angle(indicesSlip(i)))));
	 end
    end
end




