function [c_attach mt_centers MT_length MT_angle X_MT mt_vec MT_state n_MT] = delete_MTs(n_MT, MT_length, c_attach, mt_centers, MT_angle, X_MT, mt_vec, MT_state, c_rad);

LocDeleteMT=find(MT_length<=0);
deleteMT=length(LocDeleteMT);
if(deleteMT>0)
    c_attach(LocDeleteMT)=[]; 
    mt_centers(:,LocDeleteMT) = [];
    MT_length(LocDeleteMT)=[];
    MT_angle(LocDeleteMT) = [];
    X_MT(:,LocDeleteMT) = [];
    mt_vec(:,LocDeleteMT) =[];
    MT_state(LocDeleteMT) = [];
end

n_MT=n_MT-deleteMT;
X_MT = [(MT_length+c_rad).*cos(MT_angle)+mt_centers(1,:); (MT_length+c_rad).*sin(MT_angle)+mt_centers(2,:)];
mt_vec(1,:)=(1./(MT_length)).*(X_MT(1,:)-(mt_centers(1,:)+c_rad*cos(MT_angle)));
        %unit vector -> scale or divide by length
        %x-location of mt minus x-location of peri centriolar material
mt_vec(2,:)=(1./(MT_length)).*(X_MT(2,:)-(mt_centers(2,:)+c_rad*sin(MT_angle)));
