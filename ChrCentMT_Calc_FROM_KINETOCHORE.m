function [ChrCentMTInteractions indices_mt_intersect]=ChrCentMT_Calc(n_chr, n_MT, X_MT, ChrBlobPts, probChr, c_attach, X_chr, centers, r, MT_length, MT_angle);

ChrCentMTInteractions=zeros(5,n_MT);%2xn_MT, 1st row=chr #, 2nd row centrosome #. 0 if no interaction, # otherwise

for i=1:n_chr
                    %query points are X_MT to check if end of MTs are in
                    %chromosome blob region
    [in] = inpolygon(X_MT(1,:),X_MT(2,:),ChrBlobPts(1,:,i),ChrBlobPts(2,:,i));
%SHOULD think about how to deal with a single MT interacting with more than
    
    randvec=rand(1,n_MT);
    indices_chr_B=and(in,(randvec<probChr));
    indices_chr_S=and(in,~and(in,indices_chr_B));
    
    ChrCentMTInteractions(1,in)=i;
    ChrCentMTInteractions(2,in)=c_attach(in);
    ChrCentMTInteractions(3,indices_chr_B)=1;
    ChrCentMTInteractions(4,indices_chr_S)=1;
    
    % ~~~~ Used for Chr Arms & For KTs 
    xstar=X_chr(1,i)-centers(1,:);
    ystar=X_chr(2,i)-centers(2,:);
    chr_angle=mod(atan2(ystar,xstar),2*pi);
    
    minL=sqrt(xstar.^2+ystar.^2)+mean(r)/2;
    indices_mt_pass=MT_length>minL(c_attach);
    %any length mt that is longer than the minimum length
    indices_mt_angle=(MT_angle>(chr_angle(c_attach)-pi/32)&MT_angle<(chr_angle(c_attach)+pi/32));
    %any angle within area of the angle of the chr to get any mts
    %near.
    indices_mt_intersect=and(indices_mt_pass,indices_mt_angle);
    
     ChrCentMTInteractions(1,indices_mt_intersect)=i;
     ChrCentMTInteractions(2,indices_mt_intersect)=c_attach(indices_mt_intersect);
     ChrCentMTInteractions(3,indices_mt_intersect)=1;
     ChrCentMTInteractions(5,indices_mt_intersect)=1;
    %intersection between those that both pass the min length and are in
    %that specific angle. 
end