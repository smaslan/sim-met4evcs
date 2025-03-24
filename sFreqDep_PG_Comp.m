function [Uc1, first, last] = sFreqDep_PG_Comp(U1, fft_size, CmpVector1)
% Frequency Dependant Phase and Gain Compensation
%
% Input arguments:
% U1         : Sampling buffer data array (uncompensated) 
% fft_size   : FFT size
% CmpVector1 : Complex vector holding Gain and Phase compensation

% Output: 
% Uc1         : Compensated output data array
% first, last : Index of which part of the input buffer is used for output
%
% This is part of the PWRTDI - TimeDomainIntegration power alg.
% Part of TWM tool: https://github.com/smaslan/TWM
% (c) 2018, Kristian Ellingsberg
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%
% Modified: Stanislav Maslan (smaslan@cmi.cz)                
%

    if ~isequal(fft_size, 2^nextpow2(fft_size))  
        error('Compensation vector must be of size 2^N!');
    end
    N = 2^nextpow2(fft_size); % Size of FFT
    Hw = HannFcompMask(N)/2;  % Masking for overlapping processing windows
    
    % set startconditions
    C1 = zeros(1,N/4); 
    Uc1 = []; % Start with empty output array
    [Bs Be Tp]=PackMan(size(U1,2),N,N/4,1);           
    %[Bs Be Tp]
    for frame = 1:Tp
        % find frame:
        [Bs Be Tp]=PackMan(size(U1,2),N,N/4,frame);  %Find Frame Position         
        %[frame Bs Be]
        
        Uo1 = FDcomp(U1(Bs:Be),CmpVector1); % Run compensation over frame data:
        % Windowing the result 
        Uo1 = Uo1.*Hw;  
                
        if frame ~=1,
            % pick data for sum and concatenation
            Uc1 = [Uc1,C1+Uo1(N*1/4+1:N*2/4)]; %data for morhping with prew. frame
        end
        C1 = Uo1(N*2/4+1:N*3/4); %data for morhping with next frame
    end
    
    first =N/2+1;
    last = Be-N/2;

end


function  maskWindow = HannFcompMask(N)
% This is part of Frequency Dependant Phase and Gain Compensation alg. Generates padded Hann window.
% Note the window itself is normalized.
%
% N = fft data buffer size
%
% This is part of the PWRTDI - TimeDomainIntegration power alg.
% (c) 2018, Kristian Ellingsberg
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%
% Modified: Stanislav Maslan (smaslan@cmi.cz)                
%    
    
    NFFTh = 2^nextpow2(N)/4; % Next power of 2 from length of y
    
    % normalized hanning window (periodic):
    w = hanning(NFFTh*2+1);
    w = w(1:end-1)/mean(w(1:end-1));

    % generate padded window:
    maskWindow = [zeros(1,NFFTh) w.' zeros(1,NFFTh)];

end

function [FrameStart,FrameEnd,Steps]=PackMan(Data_size,Buf_size,Buf_overlap,index)
% This is part of Frequency Dependant Phase and Gain Compensation alg.
% Data_size   : The Array-length of the input data buffer
% Buf_size    : sub-array-size, who can be the FFT-size
% Buf_overlap : Num of samples shift pr. frame. Buf_size/2(pow-calc) or Buf_size/4(FDcomp)
%
% This is part of the PWRTDI - TimeDomainIntegration power alg.
% (c) 2018, Kristian Ellingsberg
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%

    FrameStart = (index - 1)*Buf_overlap + 1;
    FrameEnd = FrameStart + Buf_size - 1;
    Steps = ceil((Data_size - Buf_size)/Buf_overlap);

end

function [Uo] = FDcomp(Ui, PhicompVector)
% This is part of Frequency Dependant Phase and Gain Compensation alg.
% Compensation of Frequency-dependent gain and phase errors.
% Raw time-frequency-time domain. Windowing is done Pre. and Post. this function
% ------------------------------------------
% Ui1 Sampled data
% PhicompVecto: Complex array, containing phase and gain correction over the spectrum
%
% This is part of the PWRTDI - TimeDomainIntegration power alg.
% (c) 2018, Kristian Ellingsberg
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%
% Modified: Stanislav Maslan (smaslan@cmi.cz)                
%

    % Check Input arguments
    if ~isequal(size(Ui,1), size(PhicompVector,1),1)  
        error('In function ''FDcomp(Ui,PhicompVector)'', Vector(s) not 1-domentional!');
    end
    if ~isequal(size(Ui,2), size(PhicompVector,2))  
        error('In function ''FDcomp(Ui1,PhicompVector)'', function was called with vector length mismatch!');
    end
    
    % samples count
    fft_size=size(Ui,2);
    
    % Check Input array-length
    if ~isequal(fft_size, 2^nextpow2(fft_size))  
        error('In function ''FDcomp(Ui1,Ui2,PhicompVector)'', Array lengt not exactly 2^N. Length mismatch!');
    end
    
    % Time-to-Frequenzy of sampled chunk Current-channel
    F_domain = fft(Ui,fft_size);  
    
    % Phase and Gain is corrected according to PhiCompVector
    Fcmpt = F_domain.* PhicompVector;  
    
    % ###todo: this should be somehow fixed, because in FFT filtering it is not possible to make phase
    %          correction of the nyquist component! So far just removing imaginary part...
    Fcmpt(fft_size/2+1) = real(Fcmpt(fft_size/2+1));
    
    % reconstruct time-domain:
    %###note: changed for Octave compatibility
    %Uo = ifft(Fcmpt,'symmetric'); % Uo2=ifft(YF,fft_size);
    Uo = real(ifft(Fcmpt)); % Uo2=ifft(YF,fft_size);

end

