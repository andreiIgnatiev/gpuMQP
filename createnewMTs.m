function [c_attachnew mt_centersnew MT_lengthnew MT_anglenew X_MTnew mt_vecnew MT_statenew] = createnewMTs(mt_nuc, nc, c_rad, centers)

newMT=mt_nuc; %nucleation rate of new MTs
    
if(newMT>0)
    c_attachnew(1:newMT)=randi([1 nc],1,newMT);
    mt_centersnew(:,1:newMT) = centers(:,c_attachnew);
    MT_lengthnew(1:newMT)=0.01*ones(1,newMT);
    MT_anglenew(1:newMT) = rand(1,newMT)*2*pi;
    X_MTnew(:,1:newMT) = [(MT_lengthnew(1:newMT)+c_rad).*cos(MT_anglenew(1:newMT))+mt_centersnew(1,(1:newMT));... 
        (MT_lengthnew(1:newMT)+c_rad).*sin(MT_anglenew(1:newMT))+mt_centersnew(2,(1:newMT))];
    mt_vecnew(1,1:newMT)=(1./(MT_lengthnew(1:newMT))).*(X_MTnew(1,1:newMT)-(mt_centersnew(1,1:newMT)+c_rad*cos(MT_anglenew(1:newMT))));
    mt_vecnew(2,1:newMT)=(1./(MT_lengthnew(1:newMT))).*(X_MTnew(2,1:newMT)-(mt_centersnew(2,1:newMT)+c_rad*sin(MT_anglenew(1:newMT))));
    MT_statenew(1:newMT) = ones(1,newMT);

end
     
    