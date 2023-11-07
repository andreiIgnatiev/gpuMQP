function [F_cent_chrS F_chr_chrS cross_slipchr]=ForceChrSlipCalculation(c_attach,ChrPrime,ChrCentMTInteractions,X_chr,n_chr,nc,MT_length,X_MT,mt_vec,kappa,f_stall);

F_cent_chrS=zeros(2,nc);
F_chr_chrS=zeros(2,n_chr);
cross_slipchr=zeros(1,n_chr);
%3xn_MT, 1st row=chr #, 2nd row is 1 if binding, 3rd row is 1 if slipping. 0 otherwise
%for chr#, since a MT can interact with more than 1 chromosome, we utilize
%facChr=nthprime(i)
%for cent#, since we know which cent each MT is attached to, just knowing
%which MT is generating chr related forces will allow us to put the correct
%force on the chromosome and centrosome
%fstar=min(f_stall,pi^2*kappa./(MT_length).^2); %Eqn 1
%fstar
%SLIPPING forces on each centrosome


%SLIPPING forces on each chromosome
for i=1:n_chr
    MTsoni=and((ChrCentMTInteractions(3,:)==1),(mod(ChrCentMTInteractions(1,:),ChrPrime(i))==0));
    for j=1:nc
        MTsonj=and(MTsoni,c_attach==j);
        x_slip=min(f_stall,pi^2*kappa./(MT_length(MTsonj)).^2).*(-mt_vec(1,MTsonj));
        y_slip=min(f_stall,pi^2*kappa./(MT_length(MTsonj)).^2).*(-mt_vec(2,MTsonj));
        F_cent_chrS(1,j)=F_cent_chrS(1,j)+sum(-x_slip);
        F_cent_chrS(2,j)=F_cent_chrS(1,j)+sum(-y_slip);
        F_chr_chrS(1,i)=F_chr_chrS(1,i)+sum(x_slip); 
        F_chr_chrS(2,i)=F_chr_chrS(2,i)+sum(y_slip);
        ri_vec=[(X_MT(1,MTsonj)-X_chr(1,i)); (X_MT(2,MTsonj)-X_chr(2,i))];
        ri_vec_u=vecnorm(ri_vec,2,1);
        ri_vec=ri_vec./ri_vec_u;
        ri_vec(3,:)=0;
        cross_slipchr(i)=cross_slipchr(i)+sum(sum(cross(ri_vec,[x_slip;y_slip;zeros(1,length(x_slip))],1)));
    end
end
    
    
    
    
    %cross(ri_vec,[x_slip;y_slip;zeros(1,length(x_slip))],1)
end