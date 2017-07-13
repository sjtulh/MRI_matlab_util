%%

cd(folder_load)


%% Load parameter from medi.log

if exist('medi.log', 'file')
    
    fid = fopen('medi.log', 'r');
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
    %     disp(tline)

        m = regexp(tline, 'Matrix Sizes: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['matrix_size = [', m{1}{1}, ']'])
            n_echo = matrix_size(4)
            matrix_size = matrix_size(1:3)
        end

        m = regexp(tline, 'Voxel Sizes: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['voxel_size = [', m{1}{1}, ']'])
        end

        m = regexp(tline, 'TEs( \(sec\))?: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['TE = [', m{1}{2}, ']'])
        end   

        m = regexp(tline, 'TR( \(sec\))?: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['TR = [', m{1}{2}, ']'])
        end    

        m = regexp(tline, 'Delta TE( \(sec\))?: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['delta_TE = [', m{1}{2}, ']'])
        end  

        m = regexp(tline, 'Field Direction: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['B0_dir = [', m{1}{1}, ']'])
        end   

        m = regexp(tline, 'Center Frequency( \(Hz\))?: (.*)$', 'tokens');
        if ~isempty(m)
            eval(['CF = [', m{1}{2}, ']'])
        end       

    end
    fclose(fid);

end

n_echo_load = length(TE);

%% Load data from C MEDI
cd temporary_dir_for_intermediate_files

if exist('iField.bin', 'file')
    fid = fopen('iField.bin','rb');
    iField_C = fread(fid, inf, '*float');
    fclose(fid);
    iField_C = reshape(iField_C, [2,n_echo_load,matrix_size]);
    iField_C = squeeze(iField_C(1,:,:,:,:) + sqrt(-1)*iField_C(2,:,:,:,:));
    iField_C = permute(iField_C, [2,3,4,1]);
end

if exist('iFreq_raw.bin', 'file')
    fid = fopen('iFreq_raw.bin','rb');
    iFreq_raw_C = fread(fid, inf, '*float');
    fclose(fid);
    iFreq_raw_C = reshape(iFreq_raw_C, matrix_size);
end

if exist('iFreq.bin', 'file')
    fid = fopen('iFreq.bin','rb');
    iFreq_C = fread(fid, inf, '*float');
    fclose(fid);
    iFreq_C = reshape(iFreq_C, matrix_size);
end

if exist('RDF.bin', 'file')
    fid = fopen('RDF.bin','rb');
    RDF_C = fread(fid, inf, '*float');
    fclose(fid);
    RDF_C = reshape(RDF_C, matrix_size);
end

if exist('iMag.bin', 'file')
    fid = fopen('iMag.bin','rb');
    iMag_C = fread(fid, inf, '*float');
    fclose(fid);
    iMag_C = reshape(iMag_C, matrix_size);
end

if exist('R2star.bin', 'file')
    fid = fopen('R2star.bin','rb');
    R2s_C = fread(fid, inf, '*float');
    fclose(fid);
    R2s_C = reshape(R2s_C, matrix_size);
end

if exist('P.bin', 'file')
    fid = fopen('P.bin','rb');
    P_C = fread(fid, inf, '*float');
    fclose(fid);
    P_C = reshape(P_C, matrix_size);
end

if exist('Mask.bin', 'file')
    fid = fopen('Mask.bin','rb');
    Mask_C = fread(fid, inf, '*int');
    fclose(fid);
    Mask_C = single(reshape(Mask_C, matrix_size));
end

if exist('N_std.bin', 'file')
    fid = fopen('N_std.bin','rb');
    N_std_C = fread(fid, inf, '*float');
    fclose(fid);
    N_std_C = reshape(N_std_C, matrix_size);
end

if exist(sprintf('recon_QSM_%02d.bin', n_echo_load), 'file')
    fid = fopen(sprintf('recon_QSM_%02d.bin', n_echo_load),'rb');
    QSM_C = fread(fid, inf, '*float');
    fclose(fid);
    QSM_C = reshape(QSM_C, matrix_size);
    QSM_C(QSM_C <= -32.75) = 0;
end


cd ..
