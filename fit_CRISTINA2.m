%% *****************************************************************************
%
%
%                       Do the fit on the whole image
%
%
% ******************************************************************************

function[SQFitImage_0TEms,TQFitImage,...
    SQFitImage_subsampled,TQFitImage_subsampled,...
    Fitresult_maps,SQ_NormvalImage,...
    TQ_NormvalImage,timeVec_0,...
    timeVec_interpol,fmin]  = fit_CRISTINA2(SQ_Images,TQ_Images,timeVec,t1,t2,BodyMask)

NCol= size(SQ_Images,1);
NLin= size(SQ_Images,2);
Fitresult_maps = zeros(NCol,NLin,7); 


%TQ-images add zero signal at TE = 0;
TQ_images_ext = zeros(NCol,NLin,length(timeVec)+1);
TQ_images_ext(:,:,1) = 0; 
TQ_images_ext(:,:,2:end) = TQ_Images;


timeVec_0 = zeros(1,length(timeVec)+1);
timeVec_0(1,2:end) = timeVec;
timeVec_interpol = linspace(min(timeVec_0),max(timeVec_0),10*length(timeVec_0));



% Start values, upper and lower bounds
%        A_SQslow, A_SQfast, A_TQ   , T2slow, T2fast,  SQoffset, TQoffset
x0    = [0.6      ,0.4      , 500,      20   ,   3     ,  0           ,0];
lb    = [0        ,0        ,   0,       8   ,   2     , -1e-7     ,-1e-7];
ub    = [500      ,500      ,5000,      55   ,   10    ,  1            ,1]; 
   
% fit functions 
% **************************************************************************** 
%fitfunction_TQ = @(x)  x(:,3).*((exp(-timeVec_0./x(:,4))-exp(-timeVec_0./x(:,5))).*(exp(-t1/x(:,4))-exp(-t1/x(:,5))).*exp(-t2/x(:,5))) + x(:,7); 

fitfunction_TQ_1 = @(x)  x(:,3).*((exp(-timeVec_0(1:11)./x(:,4))-exp(-timeVec_0(1:11)./x(:,5))).*(exp(-t1/x(:,4))-exp(-t1/x(:,5))).*exp(-t2/x(:,4))) + x(:,7); 
fitfunction_TQ_2 = @(x)  x(:,3).*((exp(-timeVec_0(12:end)./x(:,4))-exp(-timeVec_0(12:end)./x(:,5))).*(exp(-t1/x(:,4))-exp(-t1/x(:,5))).*exp(-t2/x(:,4))) + x(:,7); 

fitfunction_SQ = @(x) (x(:,1).*exp(-(timeVec+t1+t2)./x(:,4))   +  (x(:,2)).*exp(-(timeVec+t1+t2)./x(:,5))) + x(:,6);


% Options
%****************************************************************************
opts = optimoptions(@fmincon, 'Algorithm',' interior-point', 'Display','final', ...
    'TolFun',1e-12, 'TolX',1e-12, 'StepTolerance',1e-12, 'MaxFunEvals',5000,...%'ObjectiveLimit',1e-12,...
    'FiniteDifferenceStepSize',1e-10,'MaxIter',inf);


% III. define objective function, this captures Y = SQ_images(xi,yi,TEs) (measurement data)
%****************************************************************************
%fileID_SQ = fopen('SQresult_cistina_fit.txt','w');
%fileID_TQ = fopen('TQresult_cistina_fit.txt','w');

for x_VOX = 1:1:NCol  
     for y_VOX = 1:1:NLin
        
        if ( squeeze(BodyMask(x_VOX,y_VOX,:)) == 1 )
           
            voxeldataSQ   = permute(SQ_Images(x_VOX,y_VOX,:),[3 2 1]);
            normvalueSQ   = max(voxeldataSQ);
            voxeldataSQ_n = voxeldataSQ./normvalueSQ;

            voxeldataTQ   = permute(TQ_images_ext(x_VOX,y_VOX,:),[3 2 1]);
            normvalueTQ   = max(voxeldataTQ);
            voxeldataTQ_n = voxeldataTQ./normvalueTQ;
            
 

            % III. objective function 
            %****************************************************************************
            % myfit_L2_norm = @(x) sum(0.8*(voxeldataTQ_n - fitfunction_TQ(x)').^2)+ 0.2*sum((voxeldataSQ_n - fitfunction_SQ(x)').^2); %has to be one value!
            myfit_L2_norm = @(x) sum(0.8*(voxeldataTQ_n(1:11) - fitfunction_TQ_1(x)').^2)+ sum(0.2*(voxeldataTQ_n(12:end) - fitfunction_TQ_2(x)').^2) + 0.2*sum((voxeldataSQ_n - fitfunction_SQ(x)').^2); %has to be one value!
            

            
            % IV: do the fit
            %****************************************************************************
            MQ_problem = createOptimProblem('fmincon','objective',myfit_L2_norm,'x0',x0,'lb',lb,'ub',ub,'options',opts);
            gs = GlobalSearch('MaxTime',7,'NumTrialPoints',400, 'NumStageOnePoints',200);
            [x,fmin_mq] = run(gs,MQ_problem);

            
            x(1) = x(1).*normvalueSQ;%ASQslow
            x(2) = x(2).*normvalueSQ;%ASQfast
            x(3) = x(3).*normvalueTQ; %ATQ
            x(6) = x(6).*normvalueSQ;%SQoffset
            x(7) = x(7).*normvalueTQ; %TQoffset

            
            % store the fit of the voxel
            %****************************************************************************
            
            fmin(x_VOX,y_VOX,:)  = fmin_mq;
            Fitresult_maps(x_VOX,y_VOX,:)  = x;
            
            SQ_NormvalImage(x_VOX,y_VOX,:)  = normvalueSQ;
            TQ_NormvalImage(x_VOX,y_VOX,:)  = normvalueTQ;
            
            %fprintf(fileID_TQ,'%12.8f', x); fprintf(fileID_TQ,'%d',x_VOX);fprintf(fileID_TQ,' ');fprintf(fileID_TQ,'%d',y_VOX);fprintf(fileID_TQ,'%12.8f',normvalueTQ);fprintf(fileID_TQ,'\n');
                        
            SQFitImage(x_VOX,y_VOX,:) = fitfunction_SQ(x);
            TQFitImage(x_VOX,y_VOX,1:11) = fitfunction_TQ_1(x);
            TQFitImage(x_VOX,y_VOX,12:length(timeVec_0)) = fitfunction_TQ_2(x);

            % extrapolate the TE = 0ms value for the SQ image:    
            function_SQ = (x(1).*exp(-(timeVec_0+t1+t2)./x(4)) + (x(2)).*exp(-(timeVec_0+t1+t2)./x(5))) + x(6);
            SQFitImage_0TEms(x_VOX,y_VOX,:) = function_SQ; %.*normvalueSQ;
            %interpolate
            function_SQ_int = (x(1).*exp(-(timeVec_interpol+t1+t2)./x(4))  +  (x(2)).*exp(-(timeVec_interpol+t1+t2)./x(5))) + x(6);
            function_TQ_int =  x(3).*((exp(-timeVec_interpol./x(4)) - exp(-timeVec_interpol./x(5))).*(exp(-t1/x(4)) -  exp(-t1/x(5))).*exp(-t2/x(4))) + x(7); 
            SQFitImage_subsampled(x_VOX,y_VOX,:) = function_SQ_int; %.*normvalueSQ;
            TQFitImage_subsampled(x_VOX,y_VOX,:) = function_TQ_int; %.*normvalueTQ;
            
            x0 = x; % update Startvalue of fit
        end
    end
end


TQFitImage(end:NCol,end:NLin,:)=0;
%SQFitImage(end:NCol,end:NLin,:)=0;

SQFitImage_subsampled(end:NCol,end:NLin,:)=0;
TQFitImage_subsampled(end:NCol,end:NLin,:)=0;

SQFitImage_0TEms((size(SQFitImage_0TEms,1)):NCol,(size(SQFitImage_0TEms,1)):NLin,:) = 0;

SQ_NormvalImage((size(SQ_NormvalImage,1)):NCol,(size(SQ_NormvalImage,1)):NLin,:) = 0;
TQ_NormvalImage((size(TQ_NormvalImage,1)):NCol,(size(TQ_NormvalImage,1)):NLin,:) = 0;

filename_all = 'fit_CRISTINA2_result.mat';
save(filename_all)

end

