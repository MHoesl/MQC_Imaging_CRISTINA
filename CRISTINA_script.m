% CRISTINA Multi echo TQ and SQ Imaging Reco
% 
% data from Siemens magnetom VB17 (/T,3T)
%
% last update: 2019/11
% Michaela Hösl


mykspace0_Xi90 = rawdata_Xi90;
mykspace0_Xi0 = rawdata_Xi0;


%% POCS for asymmetric k-space reco instead of 0 fill:
iter=50;
watch = false; 
pst = 360/ph_inc;

ksz=size(mykspace0_Xi0);
l=length(NTEs_ms);
mykspace0_Xi0_rs =mean(reshape(mykspace0_Xi0,[ksz(1:2)  l  pst ksz(4)/pst]),5);
mykspace0_Xi90_rs=mean(reshape(mykspace0_Xi90,[ksz(1:2) l pst ksz(4)/pst]),5);
pocs_TE=11;  % pocs fails for later echoes due to low SNR
for j = 1:pocs_TE
    for i=1:pst
    [imPocs_Xi90(:,:,j,i), mykspace0_Xi90_pocs(:,:,j,i)]  = pocs(mykspace0_Xi90_rs(:,:,j,i),iter,watch);
    [imPocs_Xi0(:,:,j,i), mykspace0_Xi0_pocs(:,:,j,i)]  = pocs(mykspace0_Xi0_rs(:,:,j,i),iter,watch);
    end
end

%% Zero Filling
NCol=size(mykspace0_Xi0_pocs,1);
mykspace0_Xi0_pocs(:,:,pocs_TE:l,:) = mykspace0_Xi0_rs(:,:,pocs_TE:l,:);
mykspace0_Xi90_pocs(:,:,pocs_TE:l,:) = mykspace0_Xi0_rs(:,:,pocs_TE:l,:);

% 2. Zero fill k-space and apply 2D Hann window to get Image:
k_0fill = 2; 
[mykspace_Xi90, NCol, NLin] = zerofillkspace(mykspace0_Xi90_pocs,NCol0,NLin0, k_0fill);
[mykspace_Xi0, NCol, NLin]  = zerofillkspace(mykspace0_Xi0_pocs,NCol0,NLin0, k_0fill);

myimage_Xi90 = fft2c(mykspace_Xi90); 
myimage_Xi0  = fft2c(mykspace_Xi0); 


%% GetMQ
[SQ,TQ,B0,imask,SQ_TE,TQ_TE,ZQ_TE] = CRISTINA_getMQC(myimage_Xi0,myimage_Xi90,NTEs_ms);

as(cat(3,TQ,SQ,(TQ./SQ).*imask), 'ColorMap', 'Parula')

as(B0)
%% Fleysher Reco using the B0 offset value of the B0 map
[SQ_Fl,TQ_Fl, SQ_TE_Fl,TQ_TE_Fl,ZQ_TE_Fl,...
  SQXi90,TQXi90,SQXi0,TQXi0] = CRISTINA_getFleysherMQC(myimage_Xi0,myimage_Xi90,NTEs_ms,B0,EvoTimeInit);

as(cat(3,TQ_Fl,SQ_Fl,(TQ_Fl./SQ_Fl).*imask), 'ColorMap', 'Parula');
as(cat(3, mean(TQXi90(:,:,3:7),3),mean(TQXi0(:,:,3:7),3)), 'ColorMap', 'Parula')


%% Fit Results
[SQFitImage_0TEms,TQFitImage,...
SQFitImage_subsampled,TQFitImage_subsampled,...
SQ_fitresult_maps,TQ_fitresult_maps,...
SQ_NormvalImage,TQ_NormvalImage,NTEs_ex,timeVec_interpol] = fit_CRISTINA(SQ_TE_Fl,TQ_TE_Fl, NTEs_ms,EvoT0,MixTime,NCol,NLin,imask);


