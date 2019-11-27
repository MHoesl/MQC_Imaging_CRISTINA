function [SQ_Fl,TQ_Fl,...
            SQ_TE_Fl,TQ_TE_Fl,ZQ_TE_Fl,...
            SQXi90,TQXi90,SQXi0,TQXi0]=CRISTINA_getFleysherMQC(myimage_Xi0,myimage_Xi90,TEs_ms,B0,EvoT0)

l = length(TEs_ms);

S_Xi0 =fftc(myimage_Xi0,4); 
S_Xi90 =fftc(myimage_Xi90,4);
Splus  = (1/2).*(S_Xi0 + 1i.*S_Xi90);
Sminus = (1/2).*(S_Xi0 - 1i.*S_Xi90); 

% for completeness und comparsion of the SQ and TQ results of the individual cycles:
SQXi0  = abs(S_Xi0(:,:,:,6)) + abs(S_Xi0(:,:,:,8)); 
TQXi0  = abs(S_Xi0(:,:,:,4)) + abs(S_Xi0(:,:,:,10)); 
%DQXi0 = abs(S_Xi0(:,:,:,5)) + abs(S_Xi0(:,:,:,9))
SQXi90 = abs(S_Xi0(:,:,:,6)) + abs(S_Xi0(:,:,:,8)); 
TQXi90 = abs(S_Xi0(:,:,:,4)) + abs(S_Xi90(:,:,:,10)); 
%DQXi0 = abs(S_Xi0(:,:,:,5)) + abs(S_Xi0(:,:,:,9))


%% Fleysher reco
OMEGA = repmat(B0,[1 1 l size(Splus,4)]);
Stot = (Splus.*exp(1i*OMEGA*EvoT0*1e-3) - Sminus.*exp(-1i*OMEGA*EvoT0*1e-3));
for n = 1:l
    TE_s = TEs_ms(n)*1e-3;
    Stot(:,:,n,:) = Stot(:,:,n,:).*exp(1i*OMEGA(:,:,1,:)*TE_s);
end

% Along echo time
SQ_TE_Fl=abs(Stot(:,:,:,6)) + abs(Stot(:,:,:,8)); 
TQ_TE_Fl=abs(Stot(:,:,:,4)) + abs(Stot(:,:,:,10)); 
ZQ_TE_Fl=abs(Stot(:,:,:,7));

% SQ image at the first echo time and averaged TQ image over echo 3 to 7
SQ_Fl = SQ_TE_Fl(:,:,1);
TQ_Fl = mean(TQ_TE_Fl(:,:,3:7),3);

end