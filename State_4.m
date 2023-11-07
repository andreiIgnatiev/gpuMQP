%State 4; bound to cortex

function[MT_state] = State_4(indicesb,Lb,MinDBind,xRadius,X_MT,MT_state,net_force,mt_vec,f0_dynein,c_attach,gv_astral,dt);
%maybe pass in net_force and mt_vec instead of using diff_force to
%determine staying bound/unbinding; vector projection?
% indicesb = find(MT_state==4); %pulls out MTs that are bound to cortex
% Lb = length(indicesb);   
%for MTs bound to cortex
if(Lb>0)
    for i=1:Lb
        mindistcortex=xRadius-sqrt((X_MT(1,indicesb(i))).^2+(X_MT(2,indicesb(i))).^2);
        %dp=dot(net_force(:,c_attach(indicesb(i))),mt_vec(:,indicesb(i)));%*mt_vec(:,indicesb(i));
        %magproj=sign(dp)*sqrt((mt_vec(1,indicesb(i))*dp)^2+(mt_vec(2,indicesb(i))*dp)^2)
        if mindistcortex<0
            MT_state(indicesb(i))=3;
        elseif  mindistcortex>=MinDBind
            MT_state(indicesb(i))=2;
        end
    end
end

