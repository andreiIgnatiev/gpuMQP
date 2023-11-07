function [MT_length] = MT_update(MT_state,MT_length,dt,gv_astral,sv_astral,bsv_astral);

indicesg_new = find(MT_state==1); %pulls out MTs that are growing
Lg_new = length(indicesg_new);

if Lg_new>0
    MT_length(indicesg_new) = MT_length(indicesg_new)+dt*gv_astral;
end

indicess_new = find(MT_state==2); %pulls out MTs that are shrinking
Ls_new = length(indicess_new);

if Ls_new>0
    MT_length(indicess_new) = MT_length(indicess_new)-dt*sv_astral;
end

indicesSlip_new = find(MT_state==3); %pulls out MTs that are slipping on cortex
Lslip_new = length(indicesSlip_new);

if Lslip_new>0
    MT_length(indicesSlip_new) = MT_length(indicesSlip_new)+dt*gv_astral;
end

indicesb_new = find(MT_state==4); %pulls out MTs that are bound to cortex
Lb_new=length(indicesb_new);

if Lb_new>0
    MT_length(indicesb_new) = MT_length(indicesb_new)-dt*bsv_astral;
end

