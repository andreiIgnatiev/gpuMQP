%state 1; growing

function[MT_state] = State_1(indicesg,Lg,MinDBind,MT_state,k2,dt,xRadius,X_MT,gv_astral,probD); 

%MTs that are in growing state and switch to shrinking or slipping
if (Lg>0)
    for i=1:Lg
        n2=rand(1,1);
        k2bar=1-exp(-k2(indicesg(i))*dt);
        %nb=rand(1,1);
        mindistcortex=xRadius-sqrt((X_MT(1,indicesg(i))).^2+(X_MT(2,indicesg(i))).^2);
        if mindistcortex<0
            MT_state(indicesg(i))=3;
        elseif(mindistcortex<MinDBind)
            if(n2<=probD)
                MT_state(indicesg(i))= 4; %binding on cortex
            else
               MT_state(indicesg(i))=3; %slipping on cortex
            end
        elseif((n2-k2bar)<=0)
            %Randomly choosing a subset of growing to move to shrinking  
            MT_state(indicesg(i))=2;
        end     
    end
end
