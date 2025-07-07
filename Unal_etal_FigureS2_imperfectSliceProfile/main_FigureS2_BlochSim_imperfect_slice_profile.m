
%% Hazar Benan Unal, PhD
%% Figure S2 of Unal et al. MRM 2025

%%-- This code will generate the result 
%%-- Dec-2016

clear all
close all 
disp(' '); disp(' '); disp(' '); 


%% Investigation of how slice profile affects estimated MBV

FA_2nds = [20];
TBW_sinc = 4;

for FA_2nd_idx = 1:length(FA_2nds)
    FA_2nd = FA_2nds(FA_2nd_idx);
    FA_array_degree = [3, FA_2nd];
    RFphase_increment = deg2rad( 50 ); 
    % RFphase_increment = deg2rad( 118.1 ); % this is an unstable but unbiased choice

    % sliceThickness = 8*1e-3; %m
    sliceThickness = 10; %mm    
    TR = 4.5; % ms %orig 4 
    TE = 1.7; % ms %orig 1.63
    df = 0; % Hz
    Nf = 200; % number of spin isochromats along the z axis 
    
    RFduration_sec = 1000*1e-6; % sec (600 microseconds orig)
    rf_samplingPeriod = 4; %microsec
    Gspoiler_moment = 20; % mT * ms/ m; corresponds to Deichmann MRM 2009, orig 16
    Gspoiler_effectiveAmp = 20; % mT / m %orig 20 %Hazar: change it back to 20 to make it divisible by 4

    %This calculation below is confirmed with an example in Berkeley EE225E Lecture 13 notes
    gamma = 4258;     % Hz/Gauss; note: 1 Tesla is 10^4 Gauss 
    RF_BW = TBW_sinc/RFduration_sec; %Hz    
    G_sliceSelect_amp = (RF_BW*100/sliceThickness)/gamma; %mT/m.
      
    T1_bases = [250:250:1500]%,1000];
    for T1_idx = 1:length(T1_bases)
        
        %Create the figures that will be filled in each index in for loop
      
        T1_base = T1_bases(T1_idx);
        %T1_base =         300; % ms %1500 orig
        %No change for slice profile question: Hazar
        peak_T1_change =   0;%round(T1_base/10);       % maximum +/- T1 change from T1_base in [ms] ;50 orig
        T2_base =         45; % ms 
        peak_T2_change =  0;         % maximum +/- T2 change from T1_base in [ms] 
        T2_LUT_constant = 45; % ms  (T2 is assumed constant for the LUT)

        Nex_LUT = 2000;  % the LUT is for steady-state imaging 

        LUT_T1range_true = [ (T1_base - (2.5*peak_T1_change)) : 1 : (T1_base + (2.5*peak_T1_change))  ]; 

        num_TRs_to_simulate = 1000; % # of TRs for all of the simulation period 
        num_TRs_init_fixed = 1;  % # of TRs for the initial part where T1 is assumed fixed (will be zero in real scenario) 
        R2R_period = 750;           % duration of R-R period in [ms]
        Nex_block = 1;  % # of TRs in each T1-constant time block; the T1 changes in step-woise fashion every Nex_block*TR milliseconds


        %% Exact Bloch Eq. sol. w/ realistic gradient spoiler parameters (ignores diffusion)

        TR_sec = TR*1e-3; % sec 
        TE_sec = TE*1e-3; % sec
        flips_radian = deg2rad( FA_array_degree ); 

        Gspoiler_duration = 1e-3 * (Gspoiler_moment / Gspoiler_effectiveAmp); % sec 

        x = [0]; % cm
        y = [0]; % cm
        z = linspace(-0.5,0.5,Nf); % cm
        [xx,yy,zz] = ndgrid(x,y,z); % array of spatial positions in (cm)

        off_resonance_freq = df * ones(prod(size(xx)),1); % Hz
        %% Hazar: There needs to be some array of frequecies for slice
        %excitation. This df = 0 assumes perfect excitation
%         off_resonance_freq = linspace(-RF_BW/2,RF_BW/2,Nf)'; %Hz
        %Edit: No, actually it's supposed to be 0. It's off-resonance by
        %hardware, not the RF pulse, so it's not related to bandwidth at
        %all.
        %% Hazar ends

        mode = 100; % 000 would result in printouts for every run; 100 would only print warnings; for time-evolution (time-resolved magnetization) see blochCim.m for mode definition
        sens = 1; % default coil sensitivities 

        % The following is the time points (in seconds) for a single TR: 
        % ** ** note that this ends at TE so that the Bloch Eq simul output corresponds to the ADC event (echo time during readout) ** ** 
        dt = [(TR_sec - Gspoiler_duration - RFduration_sec - TE_sec + (RFduration_sec/2)), Gspoiler_duration, RFduration_sec, TE_sec - (RFduration_sec/2)]; % all in seconds 
        
            %% Hazar: This assume at each sampling point the B1 is
            %FA/constant. Instead we need to change it dt*sinc(t)
            % B1_mag = [0, 0, flips_radian(1)/RFduration_sec/gamma/2/pi, 0]; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
            
            B1_mag_sinc = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), msinc(dt(3)*1e6/rf_samplingPeriod,TBW_sinc/4), zeros(1,dt(4)*1e6/rf_samplingPeriod)];        
            B1_mag_sinc = B1_mag_sinc * (flips_radian(2)*180/pi) / (gamma*2*pi*sum(B1_mag_sinc*rf_samplingPeriod/1e6) * 180/pi);
    
            %Now sanity check the scaling
            B1_mag = [0, 0, flips_radian(2)/RFduration_sec/gamma/2/pi, 0]; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
%             B1_mag_rect = max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod);            
            B1_mag_rect = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod), zeros(1,dt(4)*1e6/rf_samplingPeriod)];
                    
            FA_fromB1mag_rect = gamma*2*pi*sum(B1_mag_rect*rf_samplingPeriod/1e6) * 180/pi;
            if abs(FA_fromB1mag_rect - flips_radian(2)*180/pi) > 1e-4; ; disp('B1 scaling incorrect!'); hbu; end;
    
            FA_fromB1mag_sinc = gamma*2*pi*sum(B1_mag_sinc*rf_samplingPeriod/1e6) * 180/pi
            if abs(FA_fromB1mag_sinc - flips_radian(2)*180/pi) > 1e-4; disp('B1 scaling incorrect!'); hbu; end;
          
            t_allSamples = [1:length(B1_mag_sinc)]*rf_samplingPeriod/1e6;
    
            gx_allSamples = zeros(1,length(t_allSamples))';
            gy_allSamples = zeros(1,length(t_allSamples))';
            gz_allSamples = zeros(1,length(t_allSamples))';
            
            gz_allSamples(int32(dt(1)*1e6/rf_samplingPeriod)+1:int32((dt(1)+dt(2))*1e6/rf_samplingPeriod)) = 0.1*Gspoiler_effectiveAmp;
            gz_allSamples(int32((dt(1)+dt(2))*1e6/rf_samplingPeriod)+1:int32((dt(1)+dt(2)+dt(3))*1e6/rf_samplingPeriod)) = 0.1*G_sliceSelect_amp;
            gz_allSamples(int32((dt(1)+dt(2)+dt(3))*1e6/rf_samplingPeriod)+1:int32((dt(1)+dt(2)+dt(3))*1e6/rf_samplingPeriod) + dt(3)*1e6/rf_samplingPeriod/2) = -0.1*G_sliceSelect_amp;
            
            %%  Hazar: end

        disp(' '); disp(' dt (timing vector) [ms] :'); dt*1e3, disp(' '); disp(' ');
        if min(dt)<=0, disp(' '); disp([' !!!! error in dt vector !!!!']); pause; disp(' ');disp(' ');end

        % gradients in gauss units (vector of size dt): 
        gx = [0, 0, 0, 0]'; 
        gy = [0, 0, 0, 0]'; 
        gz = [0, 0.1*Gspoiler_effectiveAmp, 0, 0]'; % spoiler gradient in Gauss/cm units: 1 mT/m = 0.1 Gauss/cm --> Gspoiler_effectiveAmp is multiplied by 0.1

        %Hazar: add Gz during RF
%         gz = [0, 0.1*Gspoiler_effectiveAmp, 0.1*G_sliceSelect_amp, -0.1*G_sliceSelect_amp/2]'; % spoiler gradient in Gauss/cm units: 1 mT/m = 0.1 Gauss/cm --> Gspoiler_effectiveAmp is multiplied by 0.1

        % -- checks: 
        disp(' '); 
        disp(['TR  = ' num2str( (sum(dt)*1e3) ) ' ms ']); 
        disp(['gradient moment along spoiling direction  = ' num2str( (Gspoiler_effectiveAmp*Gspoiler_duration*1e3) ) ' mT * ms / m']); 
        disp(' ');


        T1_vec = zeros(num_TRs_to_simulate, 1); 
        T1_vec(1:num_TRs_init_fixed) = T1_base;
        T2_vec = zeros(num_TRs_to_simulate, 1); 
        T2_vec(1:num_TRs_init_fixed) = T2_base;
        parfor indexTR = (num_TRs_init_fixed+1):num_TRs_to_simulate
            time_point = TR * (indexTR - num_TRs_init_fixed); 
            time_point_R2R = mod( round(time_point), round(R2R_period)); 
            T1_vec(indexTR) = T1_base - ( peak_T1_change * sin( 2*pi*time_point_R2R/R2R_period ) ); 
            T2_vec(indexTR) = T2_base - ( peak_T2_change * sin( 2*pi*time_point_R2R/R2R_period ) ); 
        end


        %% Building look-up table (dictionary) for correcting estimated T1s (mapping apparent T1 to actual T1); note that T2 is assumed constant 

        %LUT_T1range_true = unique(T1_vec,'sorted'); % unique values of T1 in the variable-T1 simulation setting 
        LUT_T1apparent_bloch = zeros(length(LUT_T1range_true), 1); 

        tic;
        parfor indx_T1 = 1:length(LUT_T1range_true)

            T1_sec = 1e-3 * LUT_T1range_true(indx_T1); 
            T2_sec = T2_LUT_constant * 1e-3; % sec                <--- T2 is assumed constant for the LUT 

            % -- bloch equation (Nex_LUT repetitions) for flip angle #1 
            mxsol = 0; mysol = 0; mzsol = 1; % initial magnetization is normalized thermal equilibrium
            RFphase_angle = 0; RFphase_add_thisTR = RFphase_increment;
            for bloch_iter = 1:Nex_LUT
                
                B1_mag = [0, 0, flips_radian(1)/RFduration_sec/gamma/2/pi, 0]'; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
                B1 = B1_mag * exp(1i*RFphase_angle); % complex number to reflect nonzero RF phase; single coil
                % NOTE: as described on page 49 in the Handbook of Pulse Seq (Figure 2.6), with b1 denoting an RF pulse with flip=theta and phase=phi: B1_vector(b1) = b1*cos(phi) + j*b1*sin(phi) = b1*exp(j*phi) using Euler's relation
%                 [mxsol,mysol,mzsol] = blochCim(B1,[gx,gy,gz],dt,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol,mysol,mzsol); % calls the irt version of Hargreaves' simulator
                                
                    %% Hazar edits
%                     B1_mag_rect_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod), zeros(1,dt(4)*1e6/rf_samplingPeriod)];
%                     B1_allSamples = B1_mag_rect_allSamples * exp(1i*RFphase_angle);
% 
                    B1_mag_sinc_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), msinc(dt(3)*1e6/rf_samplingPeriod,TBW_sinc/4), zeros(1,dt(4)*1e6/rf_samplingPeriod)];        
                    B1_mag_sinc_allSamples = B1_mag_sinc_allSamples * (flips_radian(1)*180/pi) / (gamma*2*pi*sum(B1_mag_sinc_allSamples*rf_samplingPeriod/1e6) * 180/pi);                    
                    B1_allSamples = B1_mag_sinc_allSamples * exp(1i*RFphase_angle);

                    [mxsol,mysol,mzsol] = blochCim(B1_allSamples,[gx_allSamples,gy_allSamples,gz_allSamples],t_allSamples,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol,mysol,mzsol); % calls the irt version of Hargreaves' simulator                
                    
                    %Assume perfect spoiling (Hazar)
%                     mxsol = 0*mxsol; mysol = 0*mysol;
                    %% Hazar ends               
                
                RFphase_angle = RFphase_angle + RFphase_add_thisTR;
                RFphase_add_thisTR = RFphase_add_thisTR + RFphase_increment;
            end
            Msig_bloch_fa1 = abs( mean( (mxsol + 1i*mysol) ) ); % this ignores the phase term, that is, exp(-i*RFphase_angle) (see spgrsignal.m by B. Hargreaves)

            % -- bloch equation (Nex_LUT repetitions) for flip angle #2
            mxsol = 0; mysol = 0; mzsol = 1; % initial magnetization is normalized thermal equilibrium
            RFphase_angle = 0; RFphase_add_thisTR = RFphase_increment;
            for bloch_iter = 1:Nex_LUT
                B1_mag = [0, 0, flips_radian(2)/RFduration_sec/gamma/2/pi, 0]'; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
                B1 = B1_mag * exp(1i*RFphase_angle); % complex number to reflect nonzero RF phase; single coil
                % NOTE: as described on page 49 in the Handbook of Pulse Seq (Figure 2.6), with b1 denoting an RF pulse with flip=theta and phase=phi: B1_vector(b1) = b1*cos(phi) + j*b1*sin(phi) = b1*exp(j*phi) using Euler's relation
%                 [mxsol,mysol,mzsol] = blochCim(B1,[gx,gy,gz],dt,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol,mysol,mzsol); % calls the irt version of Hargreaves' simulator
                
                    %% Hazar edits
%                     B1_mag_rect_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod), zeros(1,dt(4)*1e6/rf_samplingPeriod)];
%                     B1_allSamples = B1_mag_rect_allSamples * exp(1i*RFphase_angle);
% % 
                    B1_mag_sinc_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), msinc(dt(3)*1e6/rf_samplingPeriod,TBW_sinc/4), zeros(1,dt(4)*1e6/rf_samplingPeriod)];        
                    B1_mag_sinc_allSamples = B1_mag_sinc_allSamples * (flips_radian(2)*180/pi) / (gamma*2*pi*sum(B1_mag_sinc_allSamples*rf_samplingPeriod/1e6) * 180/pi);
                    B1_allSamples = B1_mag_sinc_allSamples * exp(1i*RFphase_angle);

                    [mxsol,mysol,mzsol] = blochCim(B1_allSamples,[gx_allSamples,gy_allSamples,gz_allSamples],t_allSamples,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol,mysol,mzsol); % calls the irt version of Hargreaves' simulator                
                    
                    %Assume perfect spoiling (Hazar)
%                     mxsol = 0*mxsol; mysol = 0*mysol;
                    %% Hazar ends
                
                RFphase_angle = RFphase_angle + RFphase_add_thisTR;
                RFphase_add_thisTR = RFphase_add_thisTR + RFphase_increment;
            end
            Msig_bloch_fa2 = abs( mean( (mxsol + 1i*mysol) ) ); % this ignores the phase term, that is, exp(-i*RFphase_angle) (see spgrsignal.m by B. Hargreaves)

            % -- DESPOT1 fitting for computing T1; this linear matrix eq. setup is based on Deoni's DESPOT and matches the formulation in Liberman et al. JMRI 2014; DOI: 10.1002/jmri.24373
            meas_spgr = [abs(Msig_bloch_fa1), abs(Msig_bloch_fa2)]; 
            y_axis_values = meas_spgr ./ sin(flips_radian);
            x_axis_values = meas_spgr ./ tan(flips_radian); 
            C_matrix = [x_axis_values.' ones(size(x_axis_values.'))];
            % the followng solves a linear matrix equation subject to a nonnegattivity constraint for the solution vector:
            solution_vec = lsqnonneg(C_matrix,y_axis_values.');

            LUT_T1apparent_bloch(indx_T1) = -TR/log(solution_vec(1));   % in milliseconds 
            M0_estimated = solution_vec(2)/(1 - solution_vec(1));    


        %if mod(indx_T1,10)==0, disp([' T1 index = ' num2str(indx_T1) ' done (T1 = ' num2str(round(T1_sec*1e3)) ' ms); time elapsed since start: ' num2str(toc) ' [sec]']); end

        end % indx_T1
        indx_T1 = length(LUT_T1range_true);
        T1_sec = 1e-3 * LUT_T1range_true(indx_T1); 
        disp([' T1 index = ' num2str(indx_T1) ' done (T1 = ' num2str(round(T1_sec*1e3)) ' ms); time elapsed since start: ' num2str(toc) ' [sec]']);
        disp(' '); disp(' '); disp(' '); 



        %% Dynamic T1 tracking loop (apprent T1 --> LUT-corrected T1) 

        RFphase_angle = 0; RFphase_add_thisTR = RFphase_increment;
        mxsol_fa1 = 0; mysol_fa1 = 0; mzsol_fa1 = 1; % initial magnetization is normalized thermal equilibrium
        mxsol_fa2 = 0; mysol_fa2 = 0; mzsol_fa2 = 1; % initial magnetization is normalized thermal equilibrium
        Msig_bloch_fa1 = zeros(num_TRs_to_simulate, 1); 
        Msig_bloch_fa2 = zeros(num_TRs_to_simulate, 1); 
        Mz_bloch_fa1 = zeros(num_TRs_to_simulate, 1); 
        Mz_bloch_fa2 = zeros(num_TRs_to_simulate, 1); 
        T1_true_bloch = zeros(num_TRs_to_simulate, 1); % ms 
        T1_apparent_BLOCH = zeros(num_TRs_to_simulate, 1); % ms 
        indexTR = 1;


        % -- LOOP with dynamic T1 changes every (Nex_block * TR) milliseconds 
        while ( indexTR <= num_TRs_to_simulate )

            T1_sec = T1_vec(indexTR)*1e-3; % sec
            T2_sec = T2_vec(indexTR)*1e-3; % sec

            % -- bloch equation (Nex_block repetitions) with fixed T1
            for bloch_iter = 1:Nex_block

                B1_mag = [0, 0, flips_radian(1)/RFduration_sec/gamma/2/pi, 0]'; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
                B1 = B1_mag * exp(1i*RFphase_angle); % complex number to reflect nonzero RF phase; single coil
                % NOTE: as described on page 49 in the Handbook of Pulse Seq (Figure 2.6), with b1 denoting an RF pulse with flip=theta and phase=phi: B1_vector(b1) = b1*cos(phi) + j*b1*sin(phi) = b1*exp(j*phi) using Euler's relation
%                 [mxsol_fa1,mysol_fa1,mzsol_fa1] = blochCim(B1,[gx,gy,gz],dt,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol_fa1,mysol_fa1,mzsol_fa1); % calls the irt version of Hargreaves' simulator
                
                    %% Hazar edits
%                     B1_mag_rect_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod), zeros(1,dt(4)*1e6/rf_samplingPeriod)];
%                     B1_allSamples = B1_mag_rect_allSamples * exp(1i*RFphase_angle);
% 
                    B1_mag_sinc_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), msinc(dt(3)*1e6/rf_samplingPeriod,TBW_sinc/4), zeros(1,dt(4)*1e6/rf_samplingPeriod)];        
                    B1_mag_sinc_allSamples = B1_mag_sinc_allSamples * (flips_radian(1)*180/pi) / (gamma*2*pi*sum(B1_mag_sinc_allSamples*rf_samplingPeriod/1e6) * 180/pi);
                    B1_allSamples = B1_mag_sinc_allSamples * exp(1i*RFphase_angle);

                    [mxsol_fa1,mysol_fa1,mzsol_fa1] = blochCim(B1_allSamples,[gx_allSamples,gy_allSamples,gz_allSamples],t_allSamples,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol_fa1,mysol_fa1,mzsol_fa1); % calls the irt version of Hargreaves' simulator                
                    
                    %Assume perfect spoiling (Hazar)
%                     mxsol = 0*mxsol; mysol = 0*mysol;
                    %% Hazar ends               
                
                Msig_bloch_fa1(indexTR) = mean( (mxsol_fa1 + 1i*mysol_fa1) ); % this ignores the phase term, that is, exp(-i*RFphase_angle) (see spgrsignal.m by B. Hargreaves)
                Mz_bloch_fa1(indexTR) = mean(mzsol_fa1);

                B1_mag = [0, 0, flips_radian(2)/RFduration_sec/gamma/2/pi, 0]'; % B1 field magnitude in gauss units (vector of size dt); B1 is only nonzero during RFduration_sec
                B1 = B1_mag * exp(1i*RFphase_angle); % complex number to reflect nonzero RF phase; single coil
                % NOTE: as described on page 49 in the Handbook of Pulse Seq (Figure 2.6), with b1 denoting an RF pulse with flip=theta and phase=phi: B1_vector(b1) = b1*cos(phi) + j*b1*sin(phi) = b1*exp(j*phi) using Euler's relation
%                 [mxsol_fa2,mysol_fa2,mzsol_fa2] = blochCim(B1,[gx,gy,gz],dt,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol_fa2,mysol_fa2,mzsol_fa2); % calls the irt version of Hargreaves' simulator
                
                    % %% Hazar edits
%                     B1_mag_rect_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), max(B1_mag)*ones(1,dt(3)*1e6/rf_samplingPeriod), zeros(1,dt(4)*1e6/rf_samplingPeriod)];
%                     B1_allSamples = B1_mag_rect_allSamples * exp(1i*RFphase_angle);
% 
                    B1_mag_sinc_allSamples = [zeros(1,int32(dt(1)*1e6/rf_samplingPeriod)), zeros(1,dt(2)*1e6/rf_samplingPeriod), msinc(dt(3)*1e6/rf_samplingPeriod,TBW_sinc/4), zeros(1,dt(4)*1e6/rf_samplingPeriod)];        
                    B1_mag_sinc_allSamples = B1_mag_sinc_allSamples * (flips_radian(2)*180/pi) / (gamma*2*pi*sum(B1_mag_sinc_allSamples*rf_samplingPeriod/1e6) * 180/pi);
                    B1_allSamples = B1_mag_sinc_allSamples * exp(1i*RFphase_angle);

                    [mxsol_fa2,mysol_fa2,mzsol_fa2] = blochCim(B1_allSamples,[gx_allSamples,gy_allSamples,gz_allSamples],t_allSamples,T1_sec,T2_sec,off_resonance_freq,[xx(:),yy(:),zz(:)],mode,sens,mxsol_fa2,mysol_fa2,mzsol_fa2); % calls the irt version of Hargreaves' simulator                
                    
                    %Assume perfect spoiling (Hazar)
%                     mxsol = 0*mxsol; mysol = 0*mysol;                    
                    % %% Hazar ends                
                
                Msig_bloch_fa2(indexTR) = mean( (mxsol_fa2 + 1i*mysol_fa2) ); % this ignores the phase term, that is, exp(-i*RFphase_angle) (see spgrsignal.m by B. Hargreaves)
                Mz_bloch_fa2(indexTR) = mean(mzsol_fa2);

                % -- DESPOT1 fitting for computing T1; this linear matrix eq. setup is based on Deoni's DESPOT and matches the formulation in Liberman et al. JMRI 2014; DOI: 10.1002/jmri.24373
                meas_spgr = [abs(Msig_bloch_fa1(indexTR)), abs(Msig_bloch_fa2(indexTR))];
                y_axis_values = meas_spgr ./ sin(flips_radian);
                x_axis_values = meas_spgr ./ tan(flips_radian);
                C_matrix = [x_axis_values.' ones(size(x_axis_values.'))];
                % the followng solves a linear matrix equation subject to a nonnegattivity constraint for the solution vector:
                solution_vec = lsqnonneg(C_matrix,y_axis_values.');

                T1_apparent_BLOCH(indexTR) = -TR/log(solution_vec(1));   % in milliseconds
                M0_estimated = solution_vec(2)/(1 - solution_vec(1));

                T1_true_bloch(indexTR) = T1_sec*1e3;
                RFphase_angle = RFphase_angle + RFphase_add_thisTR;
                RFphase_add_thisTR = RFphase_add_thisTR + RFphase_increment;        
                indexTR = indexTR + 1;

            end    

            if mod(indexTR,100)==0, disp([' .. bloch eq simulated TR # ' num2str(indexTR) ' .. ']); end
        end
        disp(' '); disp(' '); disp(' '); 


        %% mapping "T1_apparent_BLOCH" --> "T1_corrected_BLOCH" using the look-up table (LUT)

        % LUT_T1range_true: unique values of T1 in the variable-T1 simulation setting [ms]
        % LUT_T1apparent_bloch: solution of the LUT 

        T1_corrected_BLOCH = zeros(num_TRs_to_simulate, 1);

        parfor indexTR = 1:num_TRs_to_simulate

            T1_app_thisTR = T1_apparent_BLOCH(indexTR); 

            if (T1_app_thisTR >= min(LUT_T1apparent_bloch)) && (T1_app_thisTR <= max(LUT_T1apparent_bloch))

                dummy_vec = abs( T1_app_thisTR - LUT_T1apparent_bloch ); 
                T1app_match_index = find( dummy_vec == min(dummy_vec) ); T1app_match_index = T1app_match_index(1); 
                T1_corrected_BLOCH(indexTR) = LUT_T1range_true( T1app_match_index ); 

            else
                T1_corrected_BLOCH(indexTR) = T1_app_thisTR; 
            end

        end

        %% displaying results

        time_vector = TE + ( [0:num_TRs_to_simulate-1] * TR );  % milliseconds 
        TR_start_zoom = max(num_TRs_init_fixed, 600); 

        estimated_T1s(T1_idx) = mean(T1_corrected_BLOCH(end-50:end));

        
    end %end of different T1 base values
end

%% Plot the results

MBV_groundTruth = 100*(1/T1_bases(4) - 1/T1_bases(5)) / (1/T1_bases(1) - 1/T1_bases(end));
MBV_imperfectProfile = 100*(1/estimated_T1s(4) - 1/estimated_T1s(5)) / (1/estimated_T1s(1) - 1/estimated_T1s(end));

figure, 
subplot(1,2,1)
plot(T1_bases,T1_bases,'--'); hold on; plot(T1_bases,estimated_T1s,'-o')
legend('ground truth','imperfect slice profile')
xlabel('Ground truth T1 [ms]'); ylabel('Estimated T1 [ms]')
subplot(1,2,2); 
plot([1:20],MBV_groundTruth*ones(1,20),'--'); hold on; plot(TBW_sinc,MBV_imperfectProfile,'o')
xlabel('Time bandwidth product'); ylabel('Estimated MBV [%]')
legend('ground truth','imperfect slice profile')
ylim([5 7])
