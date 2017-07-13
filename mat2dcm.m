function varargout = mat2dcm(QSM, folder_output, varargin)
% Write 3D matrix as DICOM files.
%   Examples
%   --------
%
%
%


% Parse input arguments. 
%
n_ex = 0;
n_se = 0;
siu = dicomuid;
str_desc = 'QSM';
str_sample = '';
if size(varargin,2)>0
    for k=1:size(varargin,2)
        if strcmpi(varargin{k},'exam')
            n_ex=varargin{k+1};
        end
        if strcmpi(varargin{k},'series')
            n_se=varargin{k+1};
        end
        if strcmpi(varargin{k},'siu')
            siu=varargin{k+1};
        end  
        if strcmpi(varargin{k},'description')
            str_desc=varargin{k+1};
        end     
        if strcmpi(varargin{k},'sample')
            str_sample=varargin{k+1};
        end
    end
end

% Default info
if isempty(str_sample)
    info = dicominfo_default();
else
    info = dicominfo(str_sample);
    siu = info.StudyInstanceUID;
end

% UID = dicomuid;
info.StudyInstanceUID = siu;
info.SeriesInstanceUID = dicomuid;
info.StudyID = sprintf('%05d', n_ex);
info.SeriesNumber = n_se;
info.Width = size(QSM,1);
info.Height = size(QSM,2);
info.Private_0021_104f = size(QSM,3);
info.Rows = size(QSM,2);
info.Columns = size(QSM,1);
info.WindowCenter = 0;
info.WindowWidth = 500;
info.SeriesDescription = str_desc;

QSM = int16(QSM*1000);
n_sl = size(QSM,3);

% folder_output = [folder_output '/' 'DICOM'];
mkdir(folder_output);
for idx_sl = 1:n_sl
        info.MediaStorageSOPInstanceUID = dicomuid;
        info.SOPInstanceUID = info.MediaStorageSOPInstanceUID;
        info.InstanceNumber = idx_sl;
        info.SliceLocation = idx_sl*info.SliceThickness;
        info.Filename = sprintf('ex%05dse%05dsl%05d.dcm', n_ex, n_se, idx_sl);
        dicomwrite(QSM(:,:,idx_sl),[folder_output, '/', info.Filename], info);
end