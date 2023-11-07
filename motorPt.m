function[L_i L_k]=motorPt(q,incDist,centers,PtMTsonk,boundIndex,i,k);

L_i=q*incDist;
L_k=sqrt((centers(1,k)-PtMTsonk(1,boundIndex)).^2+(centers(2,k)-PtMTsonk(2,boundIndex)).^2);
% ForceDistToCenter=sqrt(PtMTsonk(1,boundIndex).^2+PtMTsonk(2,boundIndex).^2);
