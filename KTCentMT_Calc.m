function [ChrKTCentMTInteractions]=KTCentMT_Calc(nc,ChrPrime,probKT,n_chr, n_MT, X_MT, KTPts, c_attach, X_chr, centers, Krad, MT_length, MT_angle);
slope=Inf(1,n_MT);
ChrKTCentMTInteractions=zeros(2,n_MT); %1st row=chr #, 2nd row if KT interaction is happening.
for i=1:n_chr
                    %query points are X_MT to check if end of MTs are in
                    %chromosome blob region

%SHOULD think about how to deal with a single MT interacting with more than
    [in_KT] = inpolygon(X_MT(1,:),X_MT(2,:),KTPts(1,:,i),KTPts(2,:,i));
    randvec=rand(1,n_MT);
    indices_chr_BKT=and(in_KT,(randvec<probKT));
    
    inSoloChrKT=and(in_KT,ChrKTCentMTInteractions(1,:)==0);
    inMultChrKT=and(in_KT,ChrKTCentMTInteractions(1,:)>0);
    ChrKTCentMTInteractions(1,inSoloChrKT)=ChrPrime(i); 
    ChrKTCentMTInteractions(1,inMultChrKT)=ChrPrime(i)*ChrKTCentMTInteractions(1,inMultChrKT); 
    ChrKTCentMTInteractions(2,indices_chr_BKT)=1;

    
    
    % ~~~~ Used for Chr Arms & For KTs 
    xstar=X_chr(1,i)-centers(1,:);
    ystar=X_chr(2,i)-centers(2,:);
    chr_angle=mod(atan2(ystar,xstar),2*pi);
    
    minL_KT=sqrt(xstar.^2+ystar.^2)-mean(Krad)/2;
    for j=1:nc
        %Finding the closest point on the MT to the chromosome
         CtoChrvec=[xstar(j)/sqrt(xstar(j)^2+ystar(j)^2); ystar(j)/sqrt(xstar(j)^2+ystar(j)^2)];
         slope=(X_MT(2,:)-centers(2,j))./(X_MT(1,:)-centers(1,j));
         a=-slope;
         b=ones(1,n_MT);
         c=-(X_MT(2,:)-slope.*X_MT(1,:));
         xclose=(b.*(b.*centers(1,j)-a.*centers(2,j))-a.*c)./(a.^2+b.^2);
         yclose=(a.*(-b*centers(1,j)+a*centers(2,j))-b.*c)./(a.^2+b.^2);
         distMTtoKT=sqrt((xclose-X_chr(1,i)).^2+(yclose-X_chr(2,i)).^2);
         MT_KT=(distMTtoKT<mean(Krad));
        %closest point meets distance criteria and MT is actually on j-th
        %centrosome
        indices_mt_pass=and(MT_KT,c_attach==j);
        %of the ones that meet above criteria, also those that are not
        %already accounted for in in_KT
        indices_mt_intersect=and(indices_mt_pass,~in_KT);

        mtinSoloChrKT=and(indices_mt_intersect,ChrKTCentMTInteractions(1,:)==0);
        mtinMultChrKT=and(indices_mt_intersect,ChrKTCentMTInteractions(1,:)>0);
        ChrKTCentMTInteractions(1,mtinSoloChrKT)=ChrPrime(i); 
        ChrKTCentMTInteractions(1,mtinMultChrKT)=ChrPrime(i)*ChrKTCentMTInteractions(1,mtinMultChrKT); 
        indices_chr_BKTmtin=and(indices_mt_intersect,(randvec<probKT));
       
        ChrKTCentMTInteractions(2,indices_chr_BKTmtin)=1;
     
    end
end
