function[F_stari_cortex]=Slipping_Force(nc,MT_state,c_attach,f_stall,kappa,MT_length,mt_vec);

MT_slip_cortex=find(MT_state==3); 
Centrosome_MT_slip_cortex=c_attach(MT_slip_cortex); 
F_stari_cortex=zeros(2,nc);
if(isempty(MT_slip_cortex)==0)
for i=1:nc
    MTSlipOni=find(Centrosome_MT_slip_cortex==i); 
    if(isempty(MTSlipOni)==0)
        fstari_cortex=min(f_stall,pi^2*kappa./(MT_length(MT_slip_cortex(MTSlipOni)).^2));
        F_stari_cortex(1,i)=sum(fstari_cortex.*(mt_vec(1,MT_slip_cortex(MTSlipOni))));
        F_stari_cortex(2,i)=sum(fstari_cortex.*(mt_vec(2,MT_slip_cortex(MTSlipOni))));
        clear MTSlipOni
    end
end
end
%we triple checked that we are indexing on all of the correct things. 12/19
%TO DO: double check sign/direction of force.. 12/19: we removed - sign in
%front of mt_vec (comparison with F_i_cortex)
