function [SQ,TQ,B0,imask,...
            SQ_TE,TQ_TE,ZQ_TE]=CRISTINA_getMQC(myim_Xi0,myim_Xi90,NTEs_ms)

imsz=size(myim_Xi0);

%assumes 12 steps per cycle
myimgs=reshape(myim_Xi0,[imsz(1:3) 12 imsz(4)/12]);
myimgs2=reshape(myim_Xi90,[imsz(1:3) 12 imsz(4)/12]);

%1.Quadrature combination
myic=mean(myimgs+ 1i*myimgs2,5);

%2.second FFT to obtain MQ spectrum along the phasecycle dimension: 4:
myit2=fftc(myic,4);


%3.Spectral resolution equals Number of phasesteps.
% spectral points: 1  2  3  4  5  6  7  8   9  10  11  12
% corresponding Q:-6 -5 -4 -3 -2 -1  0  1   2   3   4  5
%                  -  -  - TQ DQ SQ  -  SQ  DQ TQ -
% FFT positive and negative sides are combined except Nyquist position and 0 position:
myit3=cat(4,myit2(:,:,:,7),squeeze(my_coilcombine2(cat(5,myit2(:,:,:,8:12),myit2(:,:,:,6:-1:2)),5)),myit2(:,:,:,1));


%4.BodyMask, optional: dilate mask with strel to see a little more boundaries
imask=abs(myit3(:,:,1,2))>10*mean(abs(myit3(:)));
se = strel('disk',2);
imask2=logical(imdilate(single(imask),se));


%5. Select SQ, TQ, spectral position and TE times
SQ = abs(myit3(:,:,1,2));
TQ = mean(abs(myit3(:,:,3:7,4)),3);

SQ_TE = abs(myit3(:,:,:,2));
TQ_TE = abs(myit3(:,:,:,4));
ZQ_TE = abs(myit3(:,:,:,1));


% 6. B0 Mask
Dphis=(diff(angle(myic(:,:,1:2,:)),1,3));
B0off=median(unwrap3(squeeze(Dphis),repmat(imask,[1 1 12])),3)/diff(2*pi*1e-3*NTEs_ms(1:2));
% center back due to unwrapping
B0=B0off-mean(B0off(find(imask(:))));

% display
%as(repmat(imask,[1 1 3]).*cat(3,TQ,SQ,TQ./SQ), 'ColorMap', 'Parula');