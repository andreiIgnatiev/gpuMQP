%Could add a row 4 for binding to the kinetochore
function [ChrCentMTInteractions]=ChrCentMT_Calc(nc,ChrPrime,n_chr, n_MT, X_MT, ChrBlobPts, probChr, c_attach, X_chr, centers, r, MT_length, MT_angle);

ChrCentMTInteractions=zeros(3,n_MT);
% ChrCentMTInteractions=zeros(4,n_MT);
%3xn_MT, 1st row=chr #, 2nd row is 1 if binding, 3rd row is 1 if slipping. 0 otherwise
%for chr#, since a MT can interact with more than 1 chromosome, we utilize
%facChr=nthprime(i)
%for cent#, since we know which cent each MT is attached to, just knowing
%which MT is generating chr related forces will allow us to put the correct
%force on the chromosome and centrosome


for i=1:n_chr
    
    %query points are X_MT to check if end of MTs are in or on
                    %chromosome blob region
    [in] = inpolygon(X_MT(1,:),X_MT(2,:),ChrBlobPts(1,:,i),ChrBlobPts(2,:,i));
%SHOULD think about how to deal with a single MT interacting with more than
    
    randvec=rand(1,n_MT);
    indices_chr_B=and(in,(randvec<probChr));
    B=and(in,(randvec>=probChr));
    indices_chr_S=and(in,B);
        %     slope=(X_MT(2,Indexj(k))-centers(2,j))./(X_MT(1,Indexj(k))-centers(1,j));
%          a=-slope;
%          b=1;
%          c=-(X_MT(2,Indexj(k))-slope.*X_MT(1,Indexj(k)));
%          xclose=(b.*(b*centers(1,i)-a*centers(2,i))-a.*c)./(a.^2+b.^2);
%          yclose=(a.*(-b*centers(1,i)+a*centers(2,i))-b.*c)./(a.^2+b.^2);
%          distMTtoi=sqrt((xclose-centers(1,i)).^2+(yclose-centers(2,i)).^2);
%          MT_K=find(distMTtoi<1);
    %indices_chr_K=and(in,MT_K);
    %The above is incomplete, just an idea
    
    inSoloChr=and(in,ChrCentMTInteractions(1,:)==0);
    inMultChr=and(in,ChrCentMTInteractions(1,:)>0);
    ChrCentMTInteractions(1,inSoloChr)=ChrPrime(i); 
    ChrCentMTInteractions(1,inMultChr)=ChrPrime(i)*ChrCentMTInteractions(1,inMultChr); 
    ChrCentMTInteractions(2,indices_chr_B)=1;
    ChrCentMTInteractions(3,indices_chr_S)=1;
%     ChrCentMTInteractions(4,indices_chr_K)=1;
    % ~~~~ Used for Chr Arms & For KTs 
    xstar=X_chr(1,i)-centers(1,:);
    ystar=X_chr(2,i)-centers(2,:);
    chr_angle=mod(atan2(ystar,xstar),2*pi);
%     chr_angle=real(atan2(ystar,xstar));
%     chr_angle(chr_angle<0)=abs(chr_angle(chr_angle<0))+pi;
    minL=sqrt(xstar.^2+ystar.^2)+mean(r)/2;

    for j=1:nc
        indices_mt_pass=and(MT_length>minL(j),c_attach==j);
        %any length mt that is longer than the minimum length
        indices_mt_angle=(MT_angle>(chr_angle(j)-pi/32)&MT_angle<(chr_angle(j)+pi/32));
        %any angle within area of the angle of the chr to get any mts
        %near.
        indices_mt_intersect=and(indices_mt_pass,indices_mt_angle);
        indices_mt_intersect=and(indices_mt_intersect,~in);
%         indices_mt_intersect=setdiff(indices_mt_intersect,in);
        mtinSoloChr=and(indices_mt_intersect,ChrCentMTInteractions(1,:)==0);
        mtinMultChr=and(indices_mt_intersect,ChrCentMTInteractions(1,:)>0);
        ChrCentMTInteractions(1,mtinSoloChr)=ChrPrime(i); 
        ChrCentMTInteractions(1,mtinMultChr)=ChrPrime(i)*ChrCentMTInteractions(1,mtinMultChr); 
        indices_chr_Bmtin=and(indices_mt_intersect,(randvec<probChr));
        indices_chr_Smtin=and(indices_mt_intersect,~and(indices_mt_intersect,indices_chr_Bmtin));
%         indices_chr_Kmtin=and(indices_mt_intersect,indices_chr_K);
        ChrCentMTInteractions(2,indices_chr_Bmtin)=1;
        ChrCentMTInteractions(3,indices_chr_Smtin)=1;
%         ChrCentMTInteractions(4,indices_chr_Kmtin)=1;
    end
%     clear indices_chr_Smtin indices_chr_Bmtin indices_mt_intersect indices_mt_pass indices_mt_angle
%     clear indices_chr_B indices_chr_S
end
% find(ChrCentMTInteractions(:,2)==1)
% find(ChrCentMTInteractions(:,3)==1)

%%%
%Debugging
% find(indices_chr_S>0)=[ 19    50    67    78    88 125   137   168
% 260];
%find(indices_chr_B>0)=[6    22    31   138   245];
%find(in>0)=[6    9    22    31    50  67    78    88   125   137   138
%168   245   260]
%%5+9=14
