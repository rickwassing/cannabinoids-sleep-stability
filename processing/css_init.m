function css_init(varargin)

if nargin > 0
    domain = varargin{1};
else
    domain = 'local';
end

switch domain
    case 'sleepdrive'
        if exist('eeglab') == 0 %#ok<EXIST>
            addpath('/Volumes/sleep/Sleep/SleepSoftware/eeglab/latest'); eeglab; close all
        end
        if exist('ft_defaults') == 0 %#ok<EXIST>
            addpath('/Volumes/sleep/Sleep/SleepSoftware/fieldtrip/latest'); ft_defaults;
        end
        if exist('EEG_Processor') == 0 %#ok<EXIST>
            addpath(genpath('/Volumes/sleep/Sleep/SleepSoftware/EEG_Processor/develop'));
        end
    case 'local'
        if exist('eeglab') == 0 %#ok<EXIST>
            addpath('~/Local/eeglab/latest'); eeglab; close all
        end
        if exist('ft_defaults') == 0 %#ok<EXIST>
            addpath('~/Local/fieldtrip/latest'); ft_defaults;
        end
        if exist('EEG_Processor') == 0 %#ok<EXIST>
            addpath(genpath('~/Local/EEG_Processor/develop'));
        end
    case 'superpc'
        if exist('eeglab') == 0 %#ok<EXIST>
            addpath('S:\Sleep\SleepSoftware\eeglab\latest'); eeglab; close all
        end
        if exist('ft_defaults') == 0 %#ok<EXIST>
            addpath('S:\Sleep\SleepSoftware\fieldtrip\latest'); ft_defaults;
        end
        if exist('EEG_Processor') == 0 %#ok<EXIST>
            addpath(genpath('S:\Sleep\SleepSoftware\EEG_Processor\develop'));
        end
end

end
