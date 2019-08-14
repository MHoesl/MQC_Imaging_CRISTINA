%% *****************************************************************************
%
%
%                       Do the fit on the whole image
%
%
% ******************************************************************************

function [SQFittedData_zf,TQFittedData_zf,SQFittedData_0TEms_zf,SQ_fitresult_vector_im,TQ_fitresult_vector_im,SQ_normval_im_zf,TQ_normval_im_zf,NTEs_ex] = fit_cistina(SQ_Images_masked,TQ_Images_masked, NTEs,EvoTimeInit,MixTime,NCol,NLin,BodyMask)



SQ_fitresult_vector_im = zeros(NCol,NLin,5); 
TQ_fitresult_vector_im = zeros(NCol,NLin,4);


t1= EvoTimeInit; t2= MixTime; 
timeVec = NTEs; 


%TQ-images add zero signal at TE = 0;
TQ_images_ext = zeros(NCol,NLin,length(NTEs)+1);
TQ_images_ext(:,:,2:end) = TQ_Images_masked;

NTEs_ex = zeros(1,length(NTEs)+1);
NTEs_ex(1,2:end) = NTEs;
timeVec_0 = NTEs_ex;


%         Start values, upper and lower bounds
%            A_TQ ,                 T2slow   , T2fast,  offset, 
x0_tq    = [ 1,                        20  ,   4    ,  0.0 ];
lb_tq    = [ 0,                        10  ,   0    , -1 ];
ub_tq    = [ 20.0,                     40  ,   10   ,  1 ]; 
%              A_SQslow, A_SQfast,   T2slow  , T2fast
y0_sq      = [ 0.6      , 0.4 ,       20   ,   4      ,  0.0  ];
y_lb_sq    = [ 0.6      , 0.0 ,       10   ,   0       , -1 ];
y_ub_sq    = [ 60       , 40  ,       40   ,   40      ,  1 ];

% fit functions 
% **************************************************************************** 
fitfunction_TQ = @(x)  x(:,1).*((exp(-timeVec_0./x(:,2)) - exp(-timeVec_0./x(:,3))).*(exp(-t1/x(:,2)) -  exp(-t1/x(:,3))).*exp(-t2/x(:,2))) + x(:,4); 
fitfunction_SQ = @(y) (y(:,1).*exp(-(timeVec+t1+t2)./y(:,3))   +  y(:,2).*exp(-(timeVec+t1+t2)./y(:,4))) + y(:,5);


% Options
%****************************************************************************
opts = optimoptions(@fmincon, 'Algorithm',' interior-point', 'Display','iter', 'TolFun',1e-12, 'TolX',1e-12, 'StepTolerance',1e-12, 'MaxFunEvals',2500,'ObjectiveLimit',1e-12,...
    'FiniteDifferenceStepSize',1e-10,'MaxIter',inf);


% III. define objective function, this captures Y = SQ_images(xi,yi,TEs) (measurement data)
%****************************************************************************

fileID_SQ = fopen('SQresult_cistina_fit.txt','w');
fileID_TQ = fopen('TQresult_cistina_fit.txt','w');


for x_VOX = 1:1:NCol  
    for y_VOX = 1:1:NLin 
        

        if ( squeeze(BodyMask(x_VOX,y_VOX,:)) == 1 )
           
            voxeldataSQ  = permute(SQ_Images_masked(x_VOX,y_VOX,:),[3 2 1]);
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
            gs = GlobalSearch('MaxTime',10,'NumTrialPoints',500, 'NumStageOnePoints',50);
            [xtq,fmin_tq] = run(gs,TQ_problem);

            
            % V: do the SQ fit
            %****************************************************************************
            y0_sq(3)=xtq(2); y0_sq(4)=xtq(3);
            SQ_problem = createOptimProblem('fmincon','objective', myfit_L2_norm_SQ,'x0',y0_sq,'lb',y_lb_sq,'ub',y_ub_sq,'options',opts);

            [xsq,fmin_sq] = run(gs,SQ_problem);


            
            % store the fit of the voxel
            %****************************************************************************
            SQ_fitresult_vector_im(x_VOX,y_VOX,:)  = xsq;
            TQ_fitresult_vector_im(x_VOX,y_VOX,:)  = xtq;
            
            SQ_normval_im(x_VOX,y_VOX,:)  = normvalueSQ;
            TQ_normval_im(x_VOX,y_VOX,:)  = normvalueTQ;
            
            fprintf(fileID_TQ,'%12.8f', xtq); fprintf(fileID_TQ,'%d',x_VOX);fprintf(fileID_TQ,'%d',y_VOX);fprintf(fileID_TQ,'%12.8f',normvalueTQ);fprintf(fileID_TQ,'\n');
            fprintf(fileID_SQ,'%12.8f', xsq); fprintf(fileID_TQ,'%d',x_VOX);fprintf(fileID_TQ,'%d',y_VOX);fprintf(fileID_SQ,'%12.8f',normvalueSQ);fprintf(fileID_SQ,'\n');
             
            SQFittedData(x_VOX,y_VOX,:) = fitfunction_SQ(xsq).*normvalueSQ;
            TQFittedData(x_VOX,y_VOX,:) = fitfunction_TQ(xtq).*normvalueTQ;

            % extrapolate the TE = 0ms value for the SQ image:
            function_SQ = (xsq(1).*exp(-(timeVec_0+t1+t2)./xsq(3))  +  xsq(2).*exp(-(timeVec_0+t1+t2)./xsq(4))) + xsq(5);
            SQFittedData_0TEms(x_VOX,y_VOX,:) = function_SQ.*normvalueSQ;
            
            x0 = xsq; % update Startvalue of fit
        end
    end
end


fclose(fileID_TQ);
fclose(fileID_SQ);

filename_all = 'myresult_cistina.mat';
save(filename_all)

%%

TQFittedData_zf = TQFittedData;
TQFittedData_zf(end:NCol,end:NLin,:)=0;

SQFittedData_zf = SQFittedData;
SQFittedData_zf(end:NCol,end:NLin,:)=0;


SQFittedData_0TEms_zf= SQFittedData_0TEms;
SQFittedData_0TEms_zf((size(SQFittedData_0TEms,1)):NCol,(size(SQFittedData_0TEms,1)):NLin,:)=0;

SQ_normval_im_zf = SQ_normval_im;
SQ_normval_im_zf((size(SQ_normval_im,1)):NCol,(size(SQ_normval_im,1)):NLin,:)=0;


TQ_normval_im_zf = TQ_normval_im;
TQ_normval_im_zf((size(TQ_normval_im,1)):NCol,(size(TQ_normval_im,1)):NLin,:)=0;


end

