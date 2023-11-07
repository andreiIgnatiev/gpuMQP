function[State L_i L_k]=RemoveRepeats(State,L_i,L_k);

A=sortrows(State,2);
C=find(diff(A(:,2))==0); %Find entries that are equal to zero (this is only for bb/MTFirstNumber..I haven't seen any repeats in aa/Indexi)
if(isempty(C)==0)
    A(C+1,:)=[];
    L_i(C+1)=[];
    L_k(C+1)=[];
end       
State=A;