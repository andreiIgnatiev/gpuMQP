
function [Indexi Indexj]=Slopes_2(c_attach,centers,MT_angle,i,j)

MToni=find(c_attach==i);
MTonj=find(c_attach==j);
SlopeCtoC=(centers(2,j)-centers(2,i))/(centers(1,j)-centers(1,i));
gamma1=atan2(centers(2,j)-centers(2,i),(centers(1,j)-centers(1,i)));
gamma2=gamma1+pi;
gamma3=2*pi+gamma1;


MTRangei=find(MT_angle(MToni)<gamma1+pi/2 & MT_angle(MToni)>gamma1-pi/2);
if gamma1-pi/2<0
    MTRangei=[MTRangei find(MT_angle(MToni)<=2*pi & MT_angle(MToni)>2*pi+gamma1-pi/2)];
end
if gamma1-pi/2>2*pi
    MTRangei=[MTRangei find(MT_angle(MToni)>=0 & MT_angle(MToni)<mod(gamma1-pi/2,2*pi))];
end
if gamma1+pi/2>2*pi
    MTRangei=[MTRangei find(MT_angle(MToni)>0 & MT_angle(MToni)<mod(gamma1+pi/2,2*pi))];
end
if gamma1+pi/2<0
    MTRangei=[MTRangei find(MT_angle(MToni)<=2*pi & MT_angle(MToni)>2*pi+gamma1+pi/2)];
end
Indexi=MToni(MTRangei);

MTRangej=find(MT_angle(MTonj)<gamma2+pi/2  & MT_angle(MTonj)>gamma2-pi/2);
if gamma2-pi/2>2*pi
    MTRangej=[MTRangej find(MT_angle(MTonj)>=0 & MT_angle(MTonj)<mod(gamma2+pi/2,2*pi))];
end
if gamma2-pi/2<0
    MTRangej=[MTRangej find(MT_angle(MTonj)<=2*pi & MT_angle(MTonj)>2*pi+gamma2-pi/2)];
end
if gamma2+pi/2>2*pi
    MTRangej=[MTRangej find(MT_angle(MTonj)>=0 & MT_angle(MTonj)<mod(gamma2+pi/2,2*pi))];
end
if gamma2+pi/2<0
    MTRangej=[MTRangej find(MT_angle(MTonj)<=2*pi & MT_angle(MTonj)>2*pi+gamma2+pi/2)];
end
Indexj=MTonj(MTRangej);

