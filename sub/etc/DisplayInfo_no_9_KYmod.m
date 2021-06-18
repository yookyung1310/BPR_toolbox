function DeviceSelected = loadDisplays(varargin)
%
% RYU 2011/08/09
%


DispDev(1).name          = 'Test window (800x600)';
DispDev(1).widthPx       = 800;
DispDev(1).heightPx      = 600;
DispDev(1).widthCm       = 34.5;
DispDev(1).heightCm      = 26;
DispDev(1).distanceCm    = 87;
DispDev(1).refRateHz     = 60;
DispDev(1).degPerPx      = (2*atan(DispDev(1).heightCm/2/DispDev(1).distanceCm) * 180/pi) / DispDev(1).heightPx ;
DispDev(1).pxPerDeg      = 1/DispDev(1).degPerPx;
DispDev(1).calibFileName = [];
DispDev(1).calib         = [];

DispDev(2).name          = 'SNUBIC projector with ND filter (0022_SNUBIC_cannonPJ_NDFilter_110726.mat)';
DispDev(2).widthPx       = 1400;
DispDev(2).heightPx      = 1050;
DispDev(2).widthCm       = 34.5;
DispDev(2).heightCm      = 26;
DispDev(2).distanceCm    = 87;
DispDev(2).refRateHz     = 60;
% DispDev(2).degPerPx      = (2*(atan(DispDev(2).heightCm)/2/DispDev(2).distanceCm) * 180/pi) / DispDev(2).heightPx ;   % RYU 2011/12/08
% DispDev(2).pxPerDeg      = DispDev(2).heightPx / ((atan(DispDev(2).heightCm/DispDev(2).distanceCm) * 180/pi )) ;      % RYU 2011/12/08
DispDev(2).degPerPx      = (2*atan(DispDev(2).heightCm/2/DispDev(2).distanceCm) * 180/pi) / DispDev(2).heightPx ;
DispDev(2).pxPerDeg      = 1/DispDev(2).degPerPx;
DispDev(2).calibFileName ='0022_SNUBIC_cannonPJ_NDFilter_110726.mat';
DispDev(2).calib         = [];

DispDev(3).name          = 'SNUBIC projector (0001_SNUBIC_canonPJ_new_100809.mat)'; % VA: 16' x 12'
DispDev(3).widthPx       = 1400;
DispDev(3).heightPx      = 1050;
DispDev(3).widthCm       = 34.5;
DispDev(3).heightCm      = 26;
DispDev(3).distanceCm    = 87;
DispDev(3).refRateHz     = 60;
DispDev(3).degPerPx      = (2*atan(DispDev(3).heightCm/2/DispDev(3).distanceCm) * 180/pi) / DispDev(3).heightPx ;
DispDev(3).pxPerDeg      = 1/DispDev(3).degPerPx;
DispDev(3).calibFileName ='0001_SNUBIC_canonPJ_new_100809.mat';
DispDev(3).calib         = [];

DispDev(4).name          = 'SNUBIC tmp goggle'; % VA: 30' x 22.5'
DispDev(4).widthPx       = 800;
DispDev(4).heightPx      = 600;
DispDev(4).widthCm       = 53.2;
DispDev(4).heightCm      = 39.9;
DispDev(4).distanceCm    = 100;
DispDev(4).refRateHz     = 60;
DispDev(4).degPerPx      = (2*atan(DispDev(4).heightCm/2/DispDev(4).distanceCm) * 180/pi) / DispDev(4).heightPx ;
DispDev(4).pxPerDeg      = 1/DispDev(4).degPerPx;
DispDev(4).calibFileName ='LC-XG110-calib_012308.mat';
DispDev(4).calib         = [];

DispDev(5).name          = 'CBI@NYU Projector-Screen';
DispDev(5).widthPx       = 1024;
DispDev(5).heightPx      = 768;
DispDev(5).widthCm       = 31;
DispDev(5).heightCm      = 23;
DispDev(5).distanceCm    = 57;
DispDev(5).refRateHz     = 60;
DispDev(5).degPerPx      = ((atan(DispDev(5).widthCm)/DispDev(5).distanceCm) * 180/pi) / DispDev(5).widthPx ;
DispDev(5).pxPerDeg      = DispDev(5).widthPx / ((atan(DispDev(5).widthCm/DispDev(5).distanceCm) * 180/pi )) ;
DispDev(5).calibFileName = 'LC-XG110-calib_012308.mat';
DispDev(5).calib         = [];

DispDev(6).name          = 'MRI suite @Ajou Univ';
DispDev(6).widthPx       = 1280;
DispDev(6).heightPx      = 1024;
DispDev(6).widthCm       = 7/1024*1280;%14;
DispDev(6).heightCm      = 7;
DispDev(6).distanceCm    = 15;
DispDev(6).refRateHz     = 60;
DispDev(6).degPerPx      = (2*atan(DispDev(6).heightCm/2/DispDev(6).distanceCm) * 180/pi) / DispDev(6).heightPx ;
DispDev(6).pxPerDeg      = 1/DispDev(6).degPerPx;
DispDev(6).calibFileName ='LC-XG110-calib_012308.mat';
DispDev(6).calib         = [];


DispDev(7).name          = 'SNUBIC projector with ND filter (0001_Cannon_XEED_SX60_NewLamp_NDfilter_plato_170710.mat)';
DispDev(7).widthPx       = 1400;
DispDev(7).heightPx      = 1050;
DispDev(7).widthCm       = 39.5;
DispDev(7).heightCm      = 29.5;
DispDev(7).distanceCm    = 98;
DispDev(7).refRateHz     = 60;
% DispDev(7).degPerPx      = (2*(atan(DispDev(7).heightCm)/2/DispDev(7).distanceCm) * 180/pi) / DispDev(7).heightPx ;   % RYU 2011/12/08
% DispDev(7).pxPerDeg      = DispDev(7).heightPx / ((atan(DispDev(7).heightCm/DispDev(7).distanceCm) * 180/pi )) ;      % RYU 2011/12/08
DispDev(7).degPerPx      = (2*atan(DispDev(7).heightCm/2/DispDev(7).distanceCm) * 180/pi) / DispDev(7).heightPx ;
DispDev(7).pxPerDeg      = 1/DispDev(7).degPerPx;
DispDev(7).calibFileName ='0001_Cannon_XEED_SX60_NewLamp_NDfilter_plato_170710.mat';
DispDev(7).calib         = [];

DispDev(8).name          = 'Dell UP2414Q 4K for Dali (0001_UP2414Q_Bright00_Contrast40_Dali_170913)';
DispDev(8).widthPx       = 1920;
DispDev(8).heightPx      = 1080;
DispDev(8).widthCm       = 52.704;
DispDev(8).heightCm      = 29.646;
DispDev(8).distanceCm    = 90;
DispDev(8).refRateHz     = 60;
% DispDev(8).degPerPx      = (2*(atan(DispDev(8).heightCm)/2/DispDev(8).distanceCm) * 180/pi) / DispDev(8).heightPx ;   % RYU 2011/12/08
% DispDev(8).pxPerDeg      = DispDev(8).heightPx / ((atan(DispDev(8).heightCm/DispDev(8).distanceCm) * 180/pi )) ;      % RYU 2011/12/08
DispDev(8).degPerPx      = (2*atan(DispDev(8).heightCm/2/DispDev(8).distanceCm) * 180/pi) / DispDev(8).heightPx ;
DispDev(8).pxPerDeg      = 1/DispDev(8).degPerPx;
DispDev(8).calibFileName ='0001_UP2414Q_Bright00_Contrast40_Dali_170913.mat';
DispDev(8).calib         = [];

DispDev(9).name          = 'LG Flatron L1954TP for Rm649 (0001_FLATRON-L1954TP_Rm649_Bright40_Contrast50_Gamma0_180110)';
DispDev(9).widthPx       = 1280;
DispDev(9).heightPx      = 1024;
DispDev(9).widthCm       = 37.8;
DispDev(9).heightCm      = 30;
DispDev(9).distanceCm    = 68;
DispDev(9).refRateHz     = 60;
% DispDev(9).degPerPx      = (2*(atan(DispDev(9).heightCm)/2/DispDev(9).distanceCm) * 180/pi) / DispDev(9).heightPx ;   % RYU 2011/12/08
% DispDev(9).pxPerDeg      = DispDev(9).heightPx / ((atan(DispDev(9).heightCm/DispDev(9).distanceCm) * 180/pi )) ;      % RYU 2011/12/08
DispDev(9).degPerPx      = (2*atan(DispDev(9).heightCm/2/DispDev(9).distanceCm) * 180/pi) / DispDev(9).heightPx ;
DispDev(9).pxPerDeg      = 1/DispDev(9).degPerPx;
DispDev(9).calibFileName ='0001_FLATRON-L1954TP_Rm649_Bright40_Contrast50_Gamma0_180110.mat';
DispDev(9).calib         = [];
% 

%% Choose device & load gamma table

if nargin == 1
    dispDevNo = 9;
    
elseif nargin == 0
    nDispDev = length(DispDev);
    
    fprintf(1,'-----------------------------------------------\n');
    fprintf(1,'	Choose the Display Device from the List \n');
    fprintf(1,'-----------------------------------------------\n');
    
    for iD=1:nDispDev
        fprintf(1,' %d. %s\n',iD,DispDev(iD).name);
    end
    fprintf(1,'-----------------------------------------------\n');
    
    while 1
        dispDevNo = input('Please enter the DispDev number from the list: ');
        if isnumeric(dispDevNo) && (1<= dispDevNo && dispDevNo <= nDispDev)
            break;
        end
    end
else
    
end

% load(DispDev(dispDevNo).calibFileName)
% DispDev(dispDevNo).calib = calib;

DeviceSelected = DispDev(dispDevNo);

return;