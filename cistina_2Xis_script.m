% Multi echo cistina data
% 
% data from Siemens magnetom VB17
% double resonant head coil (NCha = 2,first one for 1H, the second one for 23Na)
%
% resolution in frequency direction is automatically doubled to evade aliasing at no extra time cost
%
% 2019/07
% Michaela


%% 1. Import k-space data
clear; close all;

% Import Data Xi90
[file_Xi90,folder_Xi90] = uigetfile({'*.dat'},'FileSelector');
myfile_Xi90 = fullfile(folder_Xi90,file_Xi90); 
[rawdata_Xi90,~,twix_Xi90] = mapVBVD_centered(myfile_Xi90); 


% Import Data Xi0
[file_Xi0,folder_Xi0] = uigetfile({'*.dat'},'FileSelector');
myfile_Xi0 = fullfile(folder_Xi0,file_Xi0); 
[rawdata_Xi0,~,twix_Xi0] = mapVBVD_centered(myfile_Xi0); 



% Get Twix Infos
[NCol0, NLin0,NCha,NAve, NRep,...
    EvoTimeInit,EvoTimeStep,MixTime,ph_inc, ...
    phase_axis, phase_axis_pi,...
    InitialPhase,NTEs,Xi] = getmytwixinfos(twix_Xi90, 3);



if (NAve == 2 && NCha == 2) %first channel is 1H, second is 23Na
    mykspace0_Xi90 = permute(rawdata_Xi90,[1 3 4 2 5]);  % x y NAve NCha NRep
    mykspace0_Xi0 = permute(rawdata_Xi0,[1 3 4 2 5]);  
    
else
    mykspace0_Xi90 = rawdata_Xi90;
    mykspace0_Xi0 = rawdata_Xi0;
end  
myimage0_Xi90 = fft2c(mykspace0_Xi90);
myimage0_Xi0 = fft2c(mykspace0_Xi0);


% 7T Marseille double tuned head coil:
if length(size(mykspace0_Xi90)) > 4
    mykspace0_Xi90 = mykspace0_Xi90(:,:,:,2,:); mykspace0_Xi90 = permute(mykspace0_Xi90,[1 2 3 5 4]);
    mykspace0_Xi0 = mykspace0_Xi0(:,:,:,2,:);   mykspace0_Xi0 = permute(mykspace0_Xi0,[1 2 3 5 4]);
    
    myimage0_Xi90 = myimage0_Xi90(:,:,:,2,:);   myimage0_Xi90  = permute(myimage0_Xi90,[1 2 3 5 4]);
    myimage0_Xi0 = myimage0_Xi0(:,:,:,2,:);     myimage0_Xi0  = permute(myimage0_Xi0,[1 2 3 5 4]);
end


%as(mykspace0_Xi90)
%as(myimage0_Xi90)
%as(myimage0_Xi0)


% 2. Zero fill k-space, apply 2D Hann window and get Image:
k_0fill = 2; symmetry = 0;
[mykspace_zf_Xi90, NCol, NLin] = zerofillkspace(mykspace0_Xi90,NCol0,NLin0, k_0fill, symmetry);

[mykspace_zf_Xi0, NCol, NLin] = zerofillkspace(mykspace0_Xi0,NCol0,NLin0, k_0fill, symmetry);

myimage_zf_Xi90 = fft2c(mykspace_zf_Xi90); 
myimage_zf_Xi0 = fft2c(mykspace_zf_Xi0); 

as(mykspace_zf_Xi90)
as(myimage_zf_Xi90)
as(myimage_zf_Xi0)


%% voxel of interest signal plot
VOI = readPoints(abs(squeeze(myimage_zf_Xi90(:,:,1))), 1);

figure;plot(real(squeeze(myimage_zf_Xi90(VOI(1,1),VOI(2,1),1,:))),'.-b')
hold on;plot(imag(squeeze(myimage_zf_Xi90(VOI(1,1),VOI(2,1),1,:))),'.-r')

plot(real(squeeze(myimage_zf_Xi0(VOI(1,1),VOI(2,1),1,:))),'.--b')
hold on;plot(imag(squeeze(myimage_zf_Xi0(VOI(1,1),VOI(2,1),1,:))),'.--r')



%% 3. Reconstruction to obtain SQ, DQ, TQ images
[Imagesignal_Xi90,Imagespec_Xi90] = reco_cistina(myimage_zf_Xi90, NLin,NCol, NRep,length(NTEs));
[Imagesignal_Xi0 ,Imagespec_Xi0]  = reco_cistina(myimage_zf_Xi0, NLin,NCol, NRep,length(NTEs));


%as(spectroscopic_im)
View4D(abs(Imagespec_Xi90(:,:,:,:)))
%View4D(abs(Imagespec_Xi0(:,:,:,:)))


%% 3. cistina SQ DQ TQ images:

%frequencies vector:
nyquist = (360/ph_inc)/2;
points = size(Imagespec_Xi90,4);
f_vec = linspace(-nyquist,+nyquist,floor(points));

dist = (points/nyquist/2);
a=nyquist-3; % works for nyquist 3 and nyquist 6; ... in general?!

p_TQl = dist*a+1      ; p_TQr = dist*(a+6) +1 ;
p_DQl = dist*(a+1)+1  ; p_DQr = dist*(a+5) +1 ; 
p_SQl = dist*(a+2)+1  ; p_SQ0 = dist*(a+3) +1  ; p_SQr = dist*(a+4) +1  ;

SQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_SQl)) + abs(Imagespec_Xi90(:,:,:,p_SQr));
DQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_DQl)) + abs(Imagespec_Xi90(:,:,:,p_DQr));
TQ_image_Xi90 = abs(Imagespec_Xi90(:,:,:,p_TQl)) + abs(Imagespec_Xi90(:,:,:,p_TQr));


SQ_image_Xi0 = abs(Imagespec_Xi0(:,:,:,p_SQl)) + abs(Imagespec_Xi0(:,:,:,p_SQr));
DQ_image_Xi0 = abs(Imagespec_Xi0(:,:,:,p_DQl)) + abs(Imagespec_Xi0(:,:,:,p_DQr));
TQ_image_Xi0 = abs(Imagespec_Xi0(:,:,:,p_TQl)) + abs(Imagespec_Xi0(:,:,:,p_TQr));


as(TQ_image_Xi90,'ColorMap','parula')
as(DQ_image_Xi90,'ColorMap','parula')
as(SQ_image_Xi90,'ColorMap','parula')

as(TQ_image_Xi0,'ColorMap','parula')
as(DQ_image_Xi0,'ColorMap','parula')
as(SQ_image_Xi0,'ColorMap','parula')


%% 1. Make Body Mask

level = multithresh(SQ_image_Xi90(:,:,1),4);
BodyMask = imquantize(SQ_image_Xi90(:,:,1),level(1));
BodyMask = imbinarize(BodyMask,1);
pixelclosing = strel('disk',2);
BodyMask = imclose(BodyMask,pixelclosing);

figure();imagesc(BodyMask);axis equal; axis tight; axis off;colormap gray


BodyMask3D = zeros(size(SQ_image_Xi90));
for l = 1:1:length(NTEs)
    BodyMask3D(:,:,l)= BodyMask;
    SQ_images_Xi90m(:,:,l) = (SQ_image_Xi90(:,:,l)) .* BodyMask;
    DQ_images_Xi90m(:,:,l) = (DQ_image_Xi90(:,:,l)) .* BodyMask;
    TQ_images_Xi90m(:,:,l) = (TQ_image_Xi90(:,:,l)) .* BodyMask;
    
    SQ_images_Xi0m(:,:,l) = (SQ_image_Xi0(:,:,l)) .* BodyMask;
    DQ_images_Xi0m(:,:,l) = (DQ_image_Xi0(:,:,l)) .* BodyMask;
    TQ_images_Xi0m(:,:,l) = (TQ_image_Xi0(:,:,l)) .* BodyMask;
end

as(SQ_images_Xi90m,'ColorMap','parula')
as(DQ_images_Xi90m,'ColorMap','parula')
as(TQ_images_Xi90m,'ColorMap','parula')

as(SQ_images_Xi0m,'ColorMap','parula')
as(DQ_images_Xi0m,'ColorMap','parula')
as(TQ_images_Xi0m,'ColorMap','parula')


%Difference
%delTQ = (TQ_images_Xi90m - TQ_images_Xi0m);
%delDQ = (DQ_images_Xi90m - DQ_images_Xi0m);
%delSQ = (SQ_images_Xi90m - SQ_images_Xi0m);


%% B0 map from multiple echo time phase data

% Phase Images
PhaseImage_TE1 = atan(imag(myimage_zf_Xi90(:,:,1,2))./real(myimage_zf_Xi90(:,:,1,2)));
PhaseImage_TE3 = atan(imag(myimage_zf_Xi90(:,:,3,2))./real(myimage_zf_Xi90(:,:,3,2)));

% conventionally 
B0map = (PhaseImage_TE1-PhaseImage_TE3) / (2*pi*(NTEs(3)-NTEs(1))*1e-3);


B0map_m = B0map.* BodyMask;
as(B0map_m)

%%
%Spec_sum_of_squares_reco  = sqrt(Imagespec_Xi0.^2 + Imagespec_Xi90.^2);
%View4D(abs(Spec_sum_of_squares_reco),'ColorMap',parula)


%% Fleysher Reco using the B0 offset value of the B0 map

Spec_plus = 1/2*(Imagespec_Xi0 + 1i*Imagespec_Xi90);
Spec_minus= 1/2*(Imagespec_Xi0 - 1i*Imagespec_Xi90); 

t1 = EvoTimeInit*1e-3;

Spec_image_total = zeros(NCol,NLin,length(NTEs),NRep);
for n = 1:1:length(NTEs)
    TE_s = NTEs(n) * 1e-3;  

    for x_VOX = 1:1:NCol  
        for y_VOX = 1:1:NLin 
            
            if ( squeeze(BodyMask(x_VOX,y_VOX,:)) == 1 )
            
                omega =abs(B0map_m(x_VOX,y_VOX)); 
                Spec_image_total(x_VOX,y_VOX,n,:) = (Spec_plus(x_VOX,y_VOX,n,:).*exp(+1i*omega*t1) - Spec_minus(x_VOX,y_VOX,n,:).*exp(-1i*omega*t1)).*exp(-1i*omega*TE_s);

            end
            
        end
    end
end

View4D(abs(Spec_image_total),'ColorMap',parula)


SQ_image_total = abs(Spec_image_total(:,:,:,p_SQl)) + abs(Spec_image_total(:,:,:,p_SQr));
DQ_image_total = abs(Spec_image_total(:,:,:,p_DQl)) + abs(Spec_image_total(:,:,:,p_DQr));
TQ_image_total = abs(Spec_image_total(:,:,:,p_TQl)) + abs(Spec_image_total(:,:,:,p_TQr));


as(SQ_image_total,'ColorMap','parula')
as(DQ_image_total,'ColorMap','parula')
as(TQ_image_total,'ColorMap','parula')



%% Multiparametric Fit

%about 70min

tic
[SQFittedData_zf,TQFittedData_zf,...
    SQFittedData_0TEms_zf,...
    SQ_fitresult_vector_im,...
    TQ_fitresult_vector_im,...
    SQ_normval_im_zf,...
    TQ_normval_im_zf, NTEs_ex] = fit_cistina(SQ_image_total,TQ_image_total, NTEs,EvoTimeInit,MixTime,NCol,NLin,BodyMask);

toc


TQ_images_0TEms = zeros(NCol,NLin,length(NTEs)+1);
TQ_images_0TEms(:,:,2:end) = TQ_image_total;

% T2 Maps:
as(SQ_fitresult_vector_im,'ColorMap','Parula')
as(TQ_fitresult_vector_im,'ColorMap','Parula')

% fitted MQC images:
as(SQFittedData_0TEms_zf,'ColorMap','Parula')
as(TQ_images_0TEms,'ColorMap','Parula')


ratio = zeros(x_VOX,y_VOX);
for x_VOX = 1:1:NCol  
    for y_VOX = 1:1:NLin 
        if (mean(squeeze(SQ_images_Xi90m(x_VOX,y_VOX,:))) ~= 0 && std(squeeze(SQ_images_Xi90m(x_VOX,y_VOX,:))) ~= 0 )
            ASQ_all(x_VOX,y_VOX,:) = (SQ_fitresult_vector_im(x_VOX,y_VOX,1)+SQ_fitresult_vector_im(x_VOX,y_VOX,2));
            ratio(x_VOX,y_VOX,:) = TQ_fitresult_vector_im(x_VOX,y_VOX,1)/ASQ_all(x_VOX,y_VOX,:);
        end 
    end
end


%% select Voxels of Interest:

N = 8;
VOIs = readPoints(abs(squeeze(myimage_zf_Xi90(:,:,1))), N);


figure;hold on;title('TQC signal')
ylabel({'amplitude [a.u]'}); xlabel({'TE [ms]'});
set(gca,'FontName','Arial','FontSize',14);set(gcf,'PaperSize',[10 6]);

for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);
    
    plot(NTEs_ex,squeeze(abs(TQ_images_0TEms(xi,yi,:))),'.-','LineWidth',1) 
    plot(NTEs_ex,squeeze(abs(TQFittedData_zf(xi,yi,:))),'r-','LineWidth',1)
    hold on
    
    maxpoint = find(abs(TQFittedData_zf(xi,yi,:))==max(abs(TQFittedData_zf(xi,yi,:))));
    T2long  = round(SQ_fitresult_vector_im(xi,yi,3),1); 
    T2short = round(SQ_fitresult_vector_im(xi,yi,4),1);
    
    if size(maxpoint,1)>1
        maxpoint = (maxpoint(1));
    end

    text(13,max(TQFittedData_zf(xi,yi,:)), num2str(i))
    text(15,max(TQFittedData_zf(xi,yi,:)),{['maxTQ=' num2str(NTEs(maxpoint)) 'ms']})
    text(30,max(TQFittedData_zf(xi,yi,:)),{['T2s=' num2str(T2short) 'ms ' ', T2l=' num2str(T2long) 'ms']})
end



figure;hold on;title('SQC signal')
ylabel({'amplitude [a.u]'}); xlabel({'TE [ms]'});
set(gca,'FontName','Arial','FontSize',14);set(gcf,'PaperSize',[10 6]);

for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);
    
    plot(NTEs,squeeze(abs(SQ_image_total(xi,yi,:))),'.-','LineWidth',1)
    plot(NTEs_ex,squeeze(abs(SQFittedData_0TEms_zf(xi,yi,:))),'r-','LineWidth',1) 
    hold on
    
    T2long  = round(SQ_fitresult_vector_im(xi,yi,3),1); 
    T2short = round(SQ_fitresult_vector_im(xi,yi,4),1);


    text(8,max(SQFittedData_0TEms_zf(xi,yi,:)), num2str(i))
    text(10,max(SQFittedData_0TEms_zf(xi,yi,:)),{['T2s=' num2str(T2short) 'ms ' ', T2l=' num2str(T2long) 'ms']})
end



%% Look at spectra for specific voxels
figure; hold on;
for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);
    
    vec = squeeze(Imagespec_Xi90(xi,yi,5,:));
    
    plot(abs(vec)./max(abs(vec)));
end

ylabel({'amplitude [a.u]'}); xlabel({'n \cdot \omega'});
set(gca,'FontName','Arial','FontSize',14,'XTick',[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]);
hold on;legend('1','2','3','4','5','6','7','8');

%% Look at spectra for specific voxels along TE
for i = 1:1:N
    yi =  VOIs(1,i);
    xi =  VOIs(2,i);
    as(squeeze(Imagespec_Xi90(xi,yi,:,:)),'ColorMap','parula');
end






%%
% vec_0P_154mM = Imagespec_Xi90(29,19,5,:);
% 
% vec_2P_100mM = Imagespec_Xi90(32,34,5,:);
% vec_2P_130mM = Imagespec_Xi90(42,34,5,:);
% 
% vec_4P_100mM = Imagespec_Xi90(23,34,5,:);
% vec_4P_154mM = Imagespec_Xi90(37,27,5,:);
% 
% vec_5P_100mM = Imagespec_Xi90(28,27,5,:);
% vec_5P_154mM = Imagespec_Xi90(37,42,5,:);
% m = [vec_0P_154mM, vec_2P_100mM,vec_2P_130mM, vec_4P_100mM, vec_4P_154mM, vec_5P_154mM, vec_5P_154mM]; 
% 
% m = permute(m,[4 2 3 1]);
% 
% figure; plot(f_vec,abs(m),'LineWidth',1) 
% ylabel({'amplitude [a.u]'}); xlabel({'n \cdot \omega'});
% set(gca,'FontName','Arial','FontSize',14,'XTick',[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]);
% hold on;legend('0% agar 154mM','2% agar 100mM','2% agar 130 mM','4% agar 100mM','4% agar 154mM','5% agar 100mM','5% agar 154mM');
% 
% 
% figure;plot(f_vec,abs(m./max(m)),'LineWidth',1)
% ylabel({'amplitude normalized[a.u]'}); xlabel({'n \cdot \omega'});
% set(gca,'FontName','Arial','FontSize',14,'XTick',[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]);
% hold on;
% legend('0% agar 154mM','2% agar 100mM','2% agar 130 mM','4% agar 100mM','4% agar 154mM','5% agar 100mM','5% agar 154mM');



