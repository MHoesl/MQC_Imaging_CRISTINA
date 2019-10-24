%-----------------------------------------------------------------------
% This function imports the raw data of the TQTPPI sequence 
% tested for Trio VB19
% 
% Michaela Hoesl 08/2018/08/22
%
%-----------------------------------------------------------------------

function [tqSpectra, tqFID,  tqFID_combined, freqVec, ...
    myspectra_1stDim, fid_data_av,...
    Sequence_length_min, Sequence_start_hours, Sequence_stop_hours, ...
    filename, NPhaseCycles,Rep,NAcq, fid_data] = reco_TQTPPI_spectro(twix,zeroFill, filterFID, filterFacPost,NPhaseSteps)


%% 1. Import the Data
%twix = mapVBVD();
timestamp = twix.image.timestamp;
%pmutime   = twix.image.pmutime;
filename = twix.image.filename;
%dataSize_all = twix.image.dataSize;

Sequence_start_hours   = min(timestamp)*2.5/1000/60/60; 
Sequence_stop_hours    = max(timestamp)*2.5/1000/60/60;
Sequence_length_min    = (max(timestamp) - min(timestamp))*2.5/1000/60; 


NCha = twix.image.NCha;
NCol = twix.image.NCol;
NAcq = twix.image.NAcq;
fid_data = twix.image(:,:,:,:,:,:,:,:,:,:);


  
Rep = size(fid_data,10);
fid_data = squeeze(fid_data); %fid_datasize = size(fid_data);



%for 7T two channel coil:
if NCha == 2
    fid_data = fid_data(:,2,:);
    NCha = 1;
    fid_data = squeeze(fid_data);
end



NPhaseCycles = NAcq/(NPhaseSteps*2);

%recover the time step from the data
EvoTimeStep = twix.hdr.Phoenix.sWiPMemBlock.alFree{3}/1000;
fSample = 1/(EvoTimeStep * 1e-3);  %evotime step is in micro seconds
nfft = NPhaseCycles * NPhaseSteps * zeroFill;
freqVec = (fSample / 2 * [-1:2/nfft:1-2/nfft]) / 1e3;
tqFID = zeros(NCha,NAcq/2*zeroFill,Rep);


%%
if NCha == 1
   
    for R = 1:1:Rep
        
    % Averaging over succeeding two Measurements
    fid_data_av = zeros(size(fid_data,1), ceil(NAcq/2),R);
    for i = 1:2:NAcq
        k = i - floor(i/2);
        fid_data_av(:,k,R) = 1/2 * (fid_data(:,i,R) + fid_data(:,i+1,R));
        fid_dataA = fid_data(:,i,R);
        fid_dataB = fid_data(:,i+1,R);
    end 
    
    fid_data_av(1,:,R) = fid_data_av(1,:,R) * .5;
    
    % 1. Phase correction
    firstDim_Spec = fftshift(fft(fid_data_av(:,1,R)));
    peakvalue = max(abs(firstDim_Spec));
    peakpos = find(~(abs(firstDim_Spec)-peakvalue));
    banana = @(x) -real(firstDim_Spec(peakpos,1)*exp(1i*x));
    [x(1),fval(1)] = fminsearch(banana,[0]);
    fid_data_av_phasecorr(:,:,R) = fid_data_av(:,:,R) * exp(1i*x(1));

    % 2. Construct T_evo_FID:
    n = 2^nextpow2(size(fid_data_av_phasecorr,1));
    myspectra_1stDim  = ifftshift(ifft(fid_data_av_phasecorr(:,:,R),n,1),1);

    %test if taking the maximum of indivial point is better than taking the
    %same position (of the first) for all fid max
    myspectra_maxima = max(abs(myspectra_1stDim));       
    maximum_position = zeros(size(myspectra_maxima));
    for i = 1:1:size(myspectra_1stDim,2)
        maximum_position(1,i) = find(  abs(myspectra_1stDim(:,i)) == myspectra_maxima(1,i)  ); 
    end
    
    tqFID(1,1:NAcq/2,R) = myspectra_1stDim(maximum_position(1,1),:);        
    
    
    
    % get zero mean  
    tqFID(1,1:NAcq/2,R) = tqFID(1,1:NAcq/2,R) - mean(tqFID(1,1:NAcq/2,R));
    tqFID(1,NAcq/2+1:end,R) = 0;                                             %figure(), plot(real(tqFID(1,:,R)))
    
   
    
    % Optional cos^2 filter of final evolution-time FID:
    if filterFID
        timeVecCos = linspace(0,1/2 * pi,fix(size(tqFID,2)/filterFacPost))';
        for idx = 1:size(tqFID,1)
            filterVec = cos(timeVecCos).^2;
            filterVec(end+1:size(tqFID,2)) = 0;
            tqFID(idx,:,R) = tqFID(idx,:,R) .* filterVec';
            tqFID(idx,length(timeVecCos)+1:end,R) = 0;
        end
    end
    
    
    % 2.Dim FFT to obtain Spectrum
    %tqSpectra = zeros(size(tqFID));
    tqFID(:,1,R) = .5 * tqFID(:,1,R);                                         
    tqSpectra(1,:,R) = fftshift(fft(squeeze(real(tqFID(1,:,R)))));
    tqFID(:,1,R) = tqFID(:,1,R) * 2;                                          %figure(), plot(freqVec,real(tqSpectra(1,:,R)))
    
   
     
    end

    
    
    tqFID_combined=tqFID;
   
    
    
elseif NCha >=1

    for R = 1:1:Rep
        % Averaging over succeeding two Measurements 
        fid_data_av = zeros(size(fid_data,1), NCha, NAcq/2,R);
        for i = 1:2:NAcq
            k = i - floor(i/2);
            %fid_data_av(:,:,k,R) = 1/2 * (fid_data(:,:,i,R) + fid_data(:,:,i+1,R));
            fid_data_av(:,1,k,R) = 1/2 * (fid_data(:,2,i,R) + fid_data(:,2,i+1,R));
            
            fid_dataA(:,k) = fid_data(:,2,i,R);
            fid_dataB(:,k) = fid_data(:,2,i+1,R);
            %fid_data_av(:,:,k,R) = (fid_data(:,:,2,R)); %for not averaging
        end
        fid_data_av(1,:,:,:) = fid_data_av(1,:,:,:) * 0.5; % multiplication of 1st point of all FIDs by 0.5 --> for FFT


        for chan = 1:1:NCha

                % 1. Phase correction of the single fids before combining:    %test = fftshift(fft(fid_data_av(:,1,:))); plot(real(squeeze(test)))
                firstDim_Spec = fftshift(fft(fid_data_av(:,1,1,R)));  
                peakvalue = max(abs(firstDim_Spec)); 
                peakpos = find(~(abs(firstDim_Spec)-peakvalue));

                banana = @(x) -real(firstDim_Spec(peakpos,1,1)*exp(1i*x));
                [x(1),fval(1)] = fminsearch(banana,[0]);
                fid_data_av_phasecorr = fid_data_av * exp(1i*x(1));


                % 2. Construct T_evo_FID:
                n = 2^nextpow2(size(fid_data_av_phasecorr,1)); 
                
                myspectra_1stDim  = ifftshift(ifft(fid_data_av_phasecorr(:,chan,:,R),n,1),1);  % plot(abs(squeeze(myspectra)))
                myspectra_1stDimA  = ifftshift(ifft(fid_dataA,n,1)); myspectra_1stDimB  = ifftshift(ifft(fid_dataB,n,1));
                
                
                myspectra_1stDim = squeeze(myspectra_1stDim(:,1,:));
                myspectra_maxima = max(abs(myspectra_1stDim));


%                 maximum_position = zeros(size(myspectra_maxima));
%                 for i = 1:1:size(myspectra_1stDim,2)
%                      maximum_position(1,i) = find(  abs(myspectra_1stDim(:,i)) == myspectra_maxima(1,i)  ); 
%                 end
                
                maximum_position = find(  abs(myspectra_1stDim(:,1)) == myspectra_maxima(1,1)  );

                % individual maxima of spec, gives not a better result, further effects ...
                % take maximum of the first spec.
                tqFID(chan,1:NAcq/2,R) = myspectra_1stDim(maximum_position(1,1),:);
                tqFID(chan,:,R) = tqFID(chan,:,R) .* exp(-1i*angle(tqFID(chan,1,R)));

                
                tqFIDA = myspectra_1stDimA(maximum_position(1,1),:);
                tqFIDB = myspectra_1stDimB(maximum_position(1,1),:);
                
                tqFID(chan,:,R) = tqFID(chan,:,R) .* exp(-1i*angle(tqFID(chan,1,R)));


                % Get zero mean       
                tqFID(chan,1:NAcq/2,R) = tqFID(chan,1:NAcq/2,R) - mean(tqFID(chan,1:NAcq/2,R));
                tqFID(chan,NAcq/2+1:end,R) = 0;

                tqFIDA = tqFIDA(chan,1:NAcq/2,R) - mean(tqFIDA(chan,1:NAcq/2,R));
                tqFIDB = tqFIDB(chan,1:NAcq/2,R) - mean(tqFIDB(chan,1:NAcq/2,R));
                tqFIDA(chan,NAcq/2+1:end,R) = 0;
                tqFIDB(chan,NAcq/2+1:end,R) = 0;
                

                % Optional cos^2 filter of final evolution-time FID:
                if filterFID
                    timeVecCos = linspace(0,1/2 * pi,fix(size(tqFID,2)/filterFacPost))';
                    for idx = 1:size(tqFID,1)
                        filterVec = cos(timeVecCos).^2;
                        filterVec(end+1:size(tqFID,2)) = 0;
                        tqFID(idx,:,R) = tqFID(idx,:,R) .* filterVec';
                        tqFID(idx,length(timeVecCos)+1:end,R) = 0;
                    end
                end





        end %end loop over channels
    
%           figure(3+R)
%           tqFID4plot=permute(tqFID,[2 1 3]);
%           subplot(3,1,1), plot(real(tqFID4plot(:,:,R)));
%           title('FID (real) for all channels'), xlabel('time [ms]'), ylabel('amplitude [a.u]');
%           subplot(3,1,2), plot(imag(tqFID4plot(:,:,R)));
%           title('FID (imag) for all channels'), xlabel('time [ms]'), ylabel('amplitude [a.u]')
%           subplot(3,1,3), plot(abs(tqFID4plot(:,:,R)));
%           title('FID (abs) for all channels' ), xlabel('time [ms]'), ylabel('amplitude [a.u]')
          

    % Combine the FIDs from the channels:
    tqFID_combined   = mean((tqFID),1);
    %figure, plot(real(tqFID_comb))

    % 2.Dim FFT to obtain Spectrum
    tqSpectra = zeros(size(tqFID_combined));
    tqFID_combined(1,1,R) = .5 * tqFID_combined(:,1,R);
    tqSpectra(1,:,R) = fftshift(fft(squeeze(real(tqFID_combined(1,:,R)))));
    tqFID_combined(1,1,R) = tqFID_combined(1,1,R) * 2;
    
    tqSpectraA(1,:,R) = fftshift(fft(squeeze(real(tqFIDA(1,:,R)))));
    tqSpectraB(1,:,R) = fftshift(fft(squeeze(real(tqFIDB(1,:,R)))));
    %figure(), plot(freqVec,real(tqSpectra))
    
    end
end


end
