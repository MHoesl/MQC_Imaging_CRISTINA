% CRISTINA Multi echo TQ and SQ Imaging Reco
% 
% data from Siemens magnetom VB17 (/T,3T)
%
% last update: 2019/10
% Michaela Hösl

%% Zero fill k-space and apply 2D Hann window to get Image:

mykspace0_Xi90 = rawdata_Xi90;
mykspace0_Xi0 = rawdata_Xi0;

k_0fill = 2; symmetry = 0;
[mykspace_zf_Xi90, NCol, NLin] = zerofillkspace(mykspace0_Xi90,NCol0,NLin0, k_0fill, symmetry);
[mykspace_zf_Xi0, NCol, NLin] = zerofillkspace(mykspace0_Xi0,NCol0,NLin0, k_0fill, symmetry);

myimage_zf_Xi90 = fft2c(mykspace_zf_Xi90); 
myimage_zf_Xi0  = fft2c(mykspace_zf_Xi0); 

%as(mykspace_zf_Xi90)
%as(myimage_zf_Xi90)


%% 3. RECO to obtain spectrum per voxel --> SQ, DQ, TQ 

myimage1 = myimage_zf_Xi90(:,:,:,1:NRep);
myimage2 = myimage_zf_Xi0(:,:,:,1:NRep);

Imagesignal_Xi90 = (myimage1 - mean(myimage1,4));
Imagesignal_Xi0  = (myimage2 - mean(myimage2,4));

Imagespec_Xi90 = fftshift(fft(Imagesignal_Xi90,NRep,4),4)./NRep;
Imagespec_Xi0  = fftshift(fft(Imagesignal_Xi0,NRep,4),4)./NRep;


%View4D(abs(Imagespec_Xi90))
%View4D(abs(Imagespec_Xi0))


%% frequencies vector:
nyquist = (360/ph_inc)/2;
spectral_resolution_points = size(Imagespec_Xi90,4);
f_vec = linspace(-nyquist,+nyquist,floor(spectral_resolution_points));

spectral_distance = (spectral_resolution_points/nyquist/2);
factor_nyquist=nyquist-3; 

p_TQl = spectral_distance*factor_nyquist+1      ; p_TQr = spectral_distance*(factor_nyquist+6) +1 ;
p_DQl = spectral_distance*(factor_nyquist+1)+1  ; p_DQr = spectral_distance*(factor_nyquist+5) +1 ; 
p_SQl = spectral_distance*(factor_nyquist+2)+1  ; p_SQr = spectral_distance*(factor_nyquist+4) +1 ; %p_SQ0 = spectral_distance*(factor_nyquist+3)+1  ;

SQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_SQl)) + abs(Imagespec_Xi90(:,:,:,p_SQr));
DQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_DQl)) + abs(Imagespec_Xi90(:,:,:,p_DQr));

SQ_image_Xi0 = abs(Imagespec_Xi0(:,:,:,p_SQl)) + abs(Imagespec_Xi0(:,:,:,p_SQr));
DQ_image_Xi0 = abs(Imagespec_Xi0(:,:,:,p_DQl)) + abs(Imagespec_Xi0(:,:,:,p_DQr));

if p_TQr > spectral_resolution_points % e.g. for ph_inc of 60% the TQ signal is at the Nyquist limit hence only at one side
   TQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_TQl));
   TQ_image_Xi0  = abs(Imagespec_Xi0(:,:,:,p_TQl));
else
    TQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_TQl)) + abs(Imagespec_Xi90(:,:,:,p_TQr));
    TQ_image_Xi0  = abs(Imagespec_Xi0(:,:,:,p_TQl)) + abs(Imagespec_Xi0(:,:,:,p_TQr));
end

figure;imagesc(SQ_image_Xi90(:,:,1));axis equal; axis off; colormap 'parula'
figure;imagesc(SQ_image_Xi0(:,:,1));axis equal; axis off; colormap 'parula'

figure;imagesc(TQ_image_Xi90(:,:,5));axis equal; axis off; colormap 'parula'
figure;imagesc(TQ_image_Xi0(:,:,5));axis equal; axis off; colormap 'parula'

%as(TQ_image_Xi90,'ColorMap','parula'); as(TQ_image_Xi0,'ColorMap','parula')
% as(DQ_image_Xi90,'ColorMap','parula'); as(DQ_image_Xi0,'ColorMap','parula')
% as(SQ_image_Xi90,'ColorMap','parula'); as(SQ_image_Xi0,'ColorMap','parula')
 

%% Make Body Mask

level = multithresh(SQ_image_Xi90(:,:,1),4);
BodyMask = imquantize(SQ_image_Xi90(:,:,1),level(1));
BodyMask = imbinarize(BodyMask,1);
pixelclosing = strel('disk',2);
BodyMask = imclose(BodyMask,pixelclosing);


%as(BodyMask)

%% B0 map from multiple echo time phase data
% needs the folder pdff-master
% (PRESCO) = Phase Regularized Estimation using Smoothing and Constrained Optimization 
% Works best with phase unwrapping (e.g. unwrap2.m & unwrap3.m from https://github.com/marcsous/unwrap).
% Outputs:
%  params.B0 is B0 (Hz)
%  params.R2 is R2* (1/s)
%  params.FF is PDFF (%)
%  params.PH is PH (rad)
%  sse is sum of squares error

figure;
gammaRatio = 11.262 / 42.57747892;
mydata = myimage_zf_Xi90(:,:,:,2);
myfieldstrength = 6.89;
[params, sse] = presco(NTEs_ms*1e-3,mydata,gammaRatio*myfieldstrength); 
B0map = params.B0;
B0map_masked = params.B0.*BodyMask;

figure;imagesc(B0map);axis equal; axis off; colormap 'gray'
hold on; visboundaries(BodyMask)


%% Fleysher Reco using the B0 offset value of the B0 map

Spec_plus = 1/2*(Imagespec_Xi0 + 1i*Imagespec_Xi90);
Spec_minus= 1/2*(Imagespec_Xi0 - 1i*Imagespec_Xi90); 

t1_s = EvoTimeInit*1e-3;

Spec_image_total = zeros(NCol,NLin,length(NTEs_ms),NRep);
for n = 1:1:length(NTEs_ms)
    TE_s = NTEs_ms(n) * 1e-3;  

    for x_VOX = 1:1:NCol  
        for y_VOX = 1:1:NLin
            
                omega =abs(B0map_masked(x_VOX,y_VOX)); 
                Spec_image_total(x_VOX,y_VOX,n,:) = (Spec_plus(x_VOX,y_VOX,n,:).*exp(+1i*omega*t1_s) - Spec_minus(x_VOX,y_VOX,n,:).*exp(-1i*omega*t1_s)).*exp(-1i*omega*TE_s);
            
        end
    end
end

%View4D(abs(Spec_image_total),'ColorMap',parula)

%% TOTAl SQ, DQ and TQ images
SQ_image_total = abs(Spec_image_total(:,:,:,p_SQl)) + abs(Spec_image_total(:,:,:,p_SQr));
DQ_image_total = abs(Spec_image_total(:,:,:,p_DQl)) + abs(Spec_image_total(:,:,:,p_DQr));

if p_TQr > spectral_resolution_points
   TQ_image_total = abs(Spec_image_total(:,:,:,p_TQl));
else
    TQ_image_total = abs(Spec_image_total(:,:,:,p_TQl)) + abs(Spec_image_total(:,:,:,p_TQr));
end


% as(SQ_image_total,'ColorMap','parula')
% as(DQ_image_total,'ColorMap','parula')
% as(TQ_image_total,'ColorMap','parula')


%save('mydata_forfit_onserver.mat', 'SQ_image_total','TQ_image_total', 'NTEs_ms','EvoTimeInit','MixTime','NCol','NLin','BodyMask')




%% Look at the Spectra for N specific voxels
N = 8;
VOIs = readPoints(abs(squeeze(myimage_zf_Xi90(:,:,1))), N);
pnts_cycle = 360/ph_inc;
pnts = size(Imagespec_Xi90,4);
NCycles = pnts/pnts_cycle;

nyq_cycle = pnts/2; 
f= linspace(-nyq_cycle,+nyq_cycle,floor(pnts));


figure; hold on;
for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);

    vec = squeeze(Imagespec_Xi0(xi,yi,4,:));
    
    plot((f/NCycles),abs(vec)./max(abs(vec)));
    %plot((f/NCycles),abs(vec));
end
ylabel({'amplitude [a.u]'}); xlabel({'n \cdot \omega'});
set(gca,'FontName','Arial','FontSize',14,'XTick',[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]);
hold on;legend('1','2','3','4','5','6','7','8');




%% Look at spectra for specific voxels along TE
for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);
    
    %as(squeeze(abs(Imagespec_Xi90(xi,yi,:,:))),'ColorMap','parula');
    %as(squeeze(abs(Imagespec_Xi0(xi,yi,:,:))),'ColorMap','parula');
    
    %as(squeeze(abs(Spec_plus(xi,yi,:,:))),'ColorMap','parula');
    %as(squeeze(abs(Spec_minus(xi,yi,:,:))),'ColorMap','parula');
    
    as(squeeze(abs(Spec_image_total(xi,yi,:,:))),'ColorMap','parula');
    
end



%% Multiparametric Fit

% alternatively load fit_result_phantom.mat

tic
[SQFitImage,TQFitImage,...
    SQFitImage_0TEms,...
    SQFitImage_subsampled,TQFitImage_subsampled,...
    SQFitImage_neg, ...
    SQ_fitresult_maps,TQ_fitresult_maps,...
    SQ_NormvalImage,TQ_NormvalImage,NTEs_extended,...
    NTEs_negative,NTEs_interpolated]  = fit_CRISTINA(SQ_image_total,TQ_image_total, NTEs_ms,EvoTimeInit,MixTime,NCol,NLin,BodyMask);
toc

%myfilename = 'myfitresult_CRISTINA.mat';
%save(myfilename)

%% TQ-images add zero signal at TE = 0;
TQ_images_0TEms = zeros(NCol,NLin,length(NTEs_ms)+1);
TQ_images_0TEms(:,:,2:end) = TQ_image_total;

% Maps:
%as(SQ_fitresult_maps,'ColorMap','Parula')
%as(TQ_fitresult_maps,'ColorMap','Parula')

% fitted MQC images:
as(SQFitImage_0TEms,'ColorMap','Parula')
as(SQFitImage_subsampled,'ColorMap','Parula')
as(TQFitImage,'ColorMap','Parula')


%% Maps:
ASQl = medfilt2((SQ_fitresult_maps(:,:,1)));
ASQs = medfilt2((SQ_fitresult_maps(:,:,2)));
ATQ = medfilt2(medfilt2(TQ_fitresult_maps(:,:,1)));
T2star_l_TQ = medfilt2(TQ_fitresult_maps(:,:,2));
T2star_s_TQ = medfilt2(TQ_fitresult_maps(:,:,3));
T2star_l_SQ = medfilt2(SQ_fitresult_maps(:,:,3));
T2star_s_SQ = medfilt2(SQ_fitresult_maps(:,:,4));
T2star_l= 1/2.*(T2star_l_SQ + T2star_l_TQ);
T2star_s= 1/2.*(T2star_s_SQ + T2star_s_TQ);


as(ASQl,'ColorMap','Parula')
as(ASQs,'ColorMap','Parula')


as(T2star_l,'ColorMap','Parula')
as(T2star_s,'ColorMap','Parula')

