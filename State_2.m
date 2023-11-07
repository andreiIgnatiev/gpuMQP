%State 2; shrinking

function[MT_state] = State_2(indicess,Ls,k1,dt,MT_state);

% indicess=find(MT_state==2);
% Ls=length(indicess);
if(Ls>0)
    probS=rand(1,Ls)-(1-exp(-k1*dt));
    MT_state(indicess(probS<=0))=1;
end
%Low probability of leaving this state (switching to growing.. about 8%
%probability)
