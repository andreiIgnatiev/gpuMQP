centers=[2 1.5; 2 3];
X_MT=[1;1];
i=1;
j=2;
k=1;
mindistCtoC=sqrt((centers(1,i)-centers(1,j)).^2+(centers(2,i)-centers(2,j)).^2);
         CtoCvec=[(centers(1,j)-centers(1,i))/mindistCtoC; (centers(2,j)-centers(2,i))/mindistCtoC];
         slope=(X_MT(2,k)-centers(2,j))./(X_MT(1,k)-centers(1,j));
         a=-slope;
         b=1;
         c=-(X_MT(2,k)-slope.*X_MT(1,k));
         xclose=(b.*(b*centers(1,i)-a*centers(2,i))-a.*c)./(a.^2+b.^2);
         yclose=(a.*(-b*centers(1,i)+a*centers(2,i))-b.*c)./(a.^2+b.^2);
         %distMTtoi=sqrt((X_MT(1,MTonj)-centers(1,i)).^2+(X_MT(2,MTonj)-centers(2,i)).^2);
         distMTtoi=sqrt((xclose-centers(1,i)).^2+(yclose-centers(2,i)).^2);
         %MTtoCAngle=acos(dot(CtoCvec,mt_vec(:,k)));
         plot(