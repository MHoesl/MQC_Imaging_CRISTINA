% Script for TQTPPI VB19 tested
% Michaela Hoesl 08/2018/08/22
%
%-----------------------------------------------------------------------

%clear;
%close all;
zeroFill = 2; % value 1 is no zero filling

%optional cos^2 filter: 0 = OFF, 1 = ON, 
filterFID = 0; filterFID_strength= 2;

twix = mapVBVD();

%everything in ms
EvoTimeInit = twix.hdr.Phoenix.sWiPMemBlock.alFree{2}/1000;
EvoTimeStep = twix.hdr.Phoenix.sWiPMemBlock.alFree{3}/1000;
PhaseIncrement = twix.hdr.Phoenix.sWiPMemBlock.adFree{1};
StartPhase = twix.hdr.Phoenix.sWiPMemBlock.adFree{2};
NPhaseSteps = 360/PhaseIncrement; % 8 for 45°


%%
[tqSpectra, tqFID           ,  tqFID_combined, ...
 freqVec  , myspectra_1stDim, fid_data_av    , ...
 Sequence_length_min, Sequence_start_hours   , ... 
 Sequence_stop_hours, filename, NPhaseCycles , ...
 Rep , NAcq , raw_data] = reco_TQTPPI_spectro(twix,zeroFill, filterFID, ...
                                              filterFID_strength,NPhaseSteps);

fprintf('EvoTime_t0= %0.2f ms EvoTimeStep==%0.2f ms and Nacq=%0.0f Max_tevo = %0.2f ms\n', ...
        EvoTimeInit,EvoTimeStep,NAcq/2, (NAcq/2*EvoTimeStep-EvoTimeInit));
 
timeVec = linspace(EvoTimeInit, EvoTimeInit + EvoTimeStep * NPhaseCycles * NPhaseSteps * zeroFill * Rep,    NAcq/2 * zeroFill); 
%timeVec = linspace(EvoTimeInit, EvoTimeStep * NPhaseCycles * NPhaseSteps * Rep,    size(tqFID,2) * Rep); 

for R = 1:1:Rep
% 1st Dimenstion FID data:
% figure(),    
% subplot(3,1,1), plot(abs(squeeze(fid_data_av(:,:,R))))
% subplot(3,1,2), plot(real(squeeze(fid_data_av(:,:,R))))
% subplot(3,1,3), plot(imag(squeeze(fid_data_av(:,:,R))))

% 1st Dimenstion spectra:
%figure(),plot(real(myspectra_1stDim(:,:)))

figure, hold on; plot(timeVec, real(tqFID(1,:,R))) ,title('FID in evolution time direction'), xlabel('evolution time [ms]'), ylabel('amplitude [a.u]')
figure, hold on; plot(freqVec,real(tqSpectra(1,:,R))),title('TQTPPI spectrum real values'), xlabel('frequency [kHz]'), ylabel('amplitude [a.u]')
%figure(3), hold on; plot(real(tqSpectra(1,:,R)))
%figure(3), hold on; plot(freqVec,abs(tqSpectra(1,:,R))), title('TQTPPI spectrum absolute values'), xlabel('frequency [kHz]'), ylabel('amplitude [a.u]')

end

%% Fit timesignal
close all
my_function_triple = @(x)x(:,3).*  sin(-2*pi*3.*x(:,4).*timeVec +  x(:,8)).* (exp(-timeVec./x(:,5)) - exp(-timeVec./x(:,6)));

if size(tqFID,3) > 1   
    timeVec = linspace(EvoTimeInit, EvoTimeStep * NPhaseCycles * NPhaseSteps * zeroFill,    NAcq/2 * zeroFill); % don't include Repetitions
    
   for i = 1:1:size(tqFID,3)
      tqFID_tr = tqFID(:,:,i)';
      [Fit_FID(i,:), x(i,:), ...
            fmin(i), FitSpec(i,:), ...
          FitSpec_norm(i,:)] = FitTQTPPI(timeVec, tqFID_tr, freqVec,EvoTimeStep, zeroFill,NPhaseSteps);
      index = find( my_function_triple(x(i,:)) == max(my_function_triple(x(i,:))));
      max_TQ_time(i,:) = timeVec(index);
   end
   
   x(:,1:6)
   mean(x(:,3)./(x(:,1)+x(:,2))), std(x(:,3)./(x(:,1)+x(:,2)))
   x(:,5:6)
   mean(max_TQ_time)
   timeVec = linspace(EvoTimeInit, EvoTimeStep * NPhaseCycles * NPhaseSteps * zeroFill*Rep, NAcq/2 * zeroFill); % include Repetitions
   plot_Trio_tqtppi((permute(tqSpectra,[3,2,1]))',freqVec,timeVec,(permute(tqFID,[3,2,1]))',zeroFill, FitSpec_norm,Fit_FID);
   
else
    tqFID_tr = tqFID';
    [Fit_FID, x, fmin, FitSpectra, FitSpectra_normalized] = FitTQTPPI(timeVec, tqFID_tr, freqVec,EvoTimeStep, zeroFill,NPhaseSteps);
    %figure(); plot(timeVec,my_function_triple(x))
    index = find( my_function_triple(x) == max(my_function_triple(x)));
    max_TQ_time = timeVec(index);
    fprintf('TQTPPI parameters: A_SQF==%0.2f, A_SQS=%0.2f, A_TQ= %0.2f, T2_1= %0.2f ms, T2_2= %0.2f ms, max_TQ_time = %0.2f ms\n',x(1),x(2), x(3)/(x(1)+x(2)), x(5),x(6),max_TQ_time);
    %plot_tqtppi(tqSpectra',freqVec,timeVec,tqFID_tr,zeroFill, FitSpectra_bc_norm,Fit_FID); 
end

%% Plot results
close all
[tqSpectra_mean,TQP_height, TQ_index, TQFreq, SQ_index] = plot_Trio_tqtppi(tqSpectra',freqVec,timeVec,tqFID',zeroFill, FitSpectra_normalized,Fit_FID, my_function_triple, x); 

