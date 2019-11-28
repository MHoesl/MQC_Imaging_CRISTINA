%% *****************************************************************************
%
%
%                       MQC CRISTINA Fit for whole image
%
%
% ******************************************************************************

function[SQFitImage_0TEms,TQFitImage,...
SQFitImage_subsampled,TQFitImage_subsampled,...
SQ_fitresult_maps,TQ_fitresult_maps,...
SQ_NormvalImage,TQ_NormvalImage,timeVec_0,timeVec_interpol]  = fit_CRISTINA(SQ_Images,TQ_Images,timeVec,t1,t2,NCol,NLin,BodyMask)


SQ_fitresult_maps = zeros(NCol,NLin,5); 
TQ_fitresult_maps = zeros(NCol,NLin,4);


%TQ-images add zero signal at TE = 0;
TQ_images_ext = zeros(NCol,NLin,length(timeVec)+1);
TQ_images_ext(:,:,2:end) = TQ_Images;


timeVec_0 = zeros(1,length(timeVec)+1);
timeVec_0(1,2:end) = timeVec;
timeVec_interpol = linspace(min(timeVec_0),max(timeVec_0),10*length(timeVec_0));



% Start values, upper and lower bounds
%            A_TQ ,                 T2long   , T2short,  offset, 
x0_tq    = [ 1,                        20  ,   4    ,  0.0 ];
lb_tq    = [ 0,                        10  ,   0    , -1 ];
ub_tq    = [ 20.0,                     30  ,   10   ,  1 ]; 
%              A_SQslow, A_SQfast,   T2long  , T2short
y0_sq      = [ 0.6      , 0.4 ,       20   ,   4      ,  0.0  ];
y_lb_sq    = [ 0.0      , 0.0 ,       10    ,  0       , -1 ];
y_ub_sq    = [ 100       , 100  ,     30   ,   10      ,  1 ];


% fit functions 
% **************************************************************************** 
fitfunction_TQ = @(x)  x(:,1).*((exp(-timeVec_0./x(:,2)) - exp(-timeVec_0./x(:,3))).*(exp(-t1/x(:,2)) -  exp(-t1/x(:,3))).*exp(-t2/x(:,2))) + x(:,4); 
fitfunction_SQ = @(y) (y(:,1).*exp(-(timeVec+t1+t2)./y(:,3))   +  y(:,2).*exp(-(timeVec+t1+t2)./y(:,4))) + y(:,5);


% Options
%****************************************************************************
opts = optimoptions(@fmincon, 'Algorithm',' interior-point', 'Display','iter', ...
    'TolFun',1e-12, 'TolX',1e-12, 'StepTolerance',1e-12, 'MaxFunEvals',2500,...
    'ObjectiveLimit',1e-12,...
    'FiniteDifferenceStepSize',1e-10,'MaxIter',inf);


% III. define objective function, this captures Y = SQ_images(xi,yi,TEs) (measurement data)
%****************************************************************************

fileID_SQ = fopen('SQresult_cistina_fit.txt','w');
fileID_TQ = fopen('TQresult_cistina_fit.txt','w');


for x_VOX =1:1:NCol  
     for y_VOX = 1:1:NLin
        

        if ( squeeze(BodyMask(x_VOX,y_VOX,:)) == 1 )
           
            voxeldataSQ  = permute(SQ_Images(x_VOX,y_VOX,:),[3 2 1]);
            normvalueSQ   = max(voxeldataSQ);
            voxeldataSQ_n = voxeldataSQ./normvalueSQ;

            voxeldataTQ  = permute(TQ_images_ext(x_VOX,y_VOX,:),[3 2 1]);
            normvalueTQ   = max(voxeldataTQ);
            voxeldataTQ_n = voxeldataTQ./normvalueTQ;
            
            

            % III. objective function 
            %****************************************************************************
            myfit_L2_norm_TQ = @(x) sum((voxeldataTQ_n - fitfunction_TQ(x)').^2); %has to be one value!
            myfit_L2_norm_SQ = @(y) sum((voxeldataSQ_n - fitfunction_SQ(y)').^2); %has to be one value!


            
            % IV: do the TQ fit
            %****************************************************************************
            TQ_problem = createOptimProblem('fmincon','objective',myfit_L2_norm_TQ,'x0',x0_tq,'lb',lb_tq,'ub',ub_tq,'options',opts);
            gs = GlobalSearch('MaxTime',10,'NumTrialPoints',800, 'NumStageOnePoints',200);
            [xtq,fmin_tq] = run(gs,TQ_problem);

            
            % V: do the SQ fit
            %****************************************************************************
            %the fitted T2 values from the TQ signal curve are used as
            %start values for the fit of the SQ signal curve
            y0_sq(3)   =  xtq(2)          ;y0_sq(4)   = xtq(3);
            
            
            SQ_problem = createOptimProblem('fmincon','objective', myfit_L2_norm_SQ,'x0',y0_sq,'lb',y_lb_sq,'ub',y_ub_sq,'options',opts);

            [xsq,fmin_sq] = run(gs,SQ_problem);


            
            % store the fit of the voxel
            %****************************************************************************
            SQ_fitresult_maps(x_VOX,y_VOX,:)  = xsq;
            TQ_fitresult_maps(x_VOX,y_VOX,:)  = xtq;
            
            SQ_NormvalImage(x_VOX,y_VOX,:)  = normvalueSQ;
            TQ_NormvalImage(x_VOX,y_VOX,:)  = normvalueTQ;
            
            fprintf(fileID_TQ,'%12.8f', xtq); fprintf(fileID_TQ,'%d',x_VOX);fprintf(fileID_TQ,' ');fprintf(fileID_TQ,'%d',y_VOX);fprintf(fileID_TQ,'%12.8f',normvalueTQ);fprintf(fileID_TQ,'\n');
            fprintf(fileID_SQ,'%12.8f', xsq); fprintf(fileID_SQ,'%d',x_VOX);fprintf(fileID_SQ,' ');fprintf(fileID_SQ,'%12.8f',normvalueSQ);fprintf(fileID_SQ,'\n');
             
            SQFitImage(x_VOX,y_VOX,:) = fitfunction_SQ(xsq).*normvalueSQ;
            TQFitImage(x_VOX,y_VOX,:) = fitfunction_TQ(xtq).*normvalueTQ;

            % extrapolate the TE = 0ms value for the SQ image:
            function_SQ = (xsq(1).*exp(-(timeVec_0+t1+t2)./xsq(3))  +  xsq(2).*exp(-(timeVec_0+t1+t2)./xsq(4))) + xsq(5);
            SQFitImage_0TEms(x_VOX,y_VOX,:) = function_SQ.*normvalueSQ;
            
            %fit data at intermediate time points
            function_SQ_int = (xsq(1).*exp(-(timeVec_interpol+t1+t2)./xsq(3))  +  xsq(2).*exp(-(timeVec_interpol+t1+t2)./xsq(4))) + xsq(5);
            function_TQ_int =  xtq(1).*((exp(-timeVec_interpol./xtq(2)) - exp(-timeVec_interpol./xtq(3))).*(exp(-t1/xtq(2)) -  exp(-t1/xtq(3))).*exp(-t2/xtq(2))) + xtq(4); 

            SQFitImage_subsampled(x_VOX,y_VOX,:) = function_SQ_int.*normvalueSQ;
            TQFitImage_subsampled(x_VOX,y_VOX,:) = function_TQ_int.*normvalueTQ;
            
            x0 = xsq; % update Startvalue of fit
        end
    end
end


fclose(fileID_TQ);
fclose(fileID_SQ);


TQFitImage(end:NCol,end:NLin,:)=0;
SQFitImage(end:NCol,end:NLin,:)=0;

SQFitImage_subsampled(end:NCol,end:NLin,:)=0;
TQFitImage_subsampled(end:NCol,end:NLin,:)=0;


SQFitImage_neg(end:NCol,end:NLin,:)=0;

SQFitImage_0TEms((size(SQFitImage_0TEms,1)):NCol,(size(SQFitImage_0TEms,1)):NLin,:) = 0;


SQ_NormvalImage((size(SQ_NormvalImage,1)):NCol,(size(SQ_NormvalImage,1)):NLin,:) = 0;
TQ_NormvalImage((size(TQ_NormvalImage,1)):NCol,(size(TQ_NormvalImage,1)):NLin,:) = 0;


%filename_all = 'myfitresult_CRISTINA.mat';
%save(filename_all)

end

