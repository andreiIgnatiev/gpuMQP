
function [Indexi Indexj CtoCVecItoj CtoCVecjtoI]=InterpolarMTs(i,j,mt_vec,centers,c_attach,MT_angle)
MTsoni=find(c_attach==i);
MTsonj=find(c_attach==j);
mindistCtoC=sqrt((centers(1,j)-centers(1,i)).^2+(centers(2,j)-centers(2,i)).^2);
CtoCVecItoj=[(centers(1,j)-centers(1,i)) (centers(2,j)-centers(2,i))]/mindistCtoC;
CtoCVecjtoI=-CtoCVecItoj;
MTRangei=[];
MTRangej=[];

Remove=[];
for k=1:length(MTsoni)
    MTRangei(k)=real(acos(dot(mt_vec(:,MTsoni(k)),CtoCVecItoj)));
     if(MTRangei(k)>pi/2)
        Remove=[Remove; k];
     end
end
  MTsoni(Remove)=[];
    Indexi=MTsoni;

Remove=[];
for l=1:length(MTsonj)
    MTRangej(l)=real(acos(dot(mt_vec(:,MTsonj(l)),CtoCVecjtoI)));
    if(MTRangej(l)>pi/2)
        Remove=[Remove; l];
    end
end
  MTsonj(Remove)=[];
    Indexj=MTsonj;
