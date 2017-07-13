function [data,kx,ky,nr,d,hdr] = read_pfile(filename)
%function [data,kx,ky,nr,d,hdr] = read_pfile(filename)  data raw kspace= zeros(xres,yres,ncoi,nech,nsli)
%
% Program to read in a p-file from the GE scanner
% If it's UTE then the kspace file is read to get
% trajectory, density etc.
%
% data is the raw kspace data
% kx is the x-trajectory
% ky is the y-trajectory
% nr is the matrix size
% d is the density weighting
% hdr is the entire GE header

% filename
if nargin==0
	[filename curr_path] = uigetfile('P*.7');
	if filename==0; return; end
	filename = [curr_path filename];
end

% read header (function below)
fid = fopen(filename,'r');
if fid==-1
    error(['File not found (' cd '\' filename ')'])
end
hdr = read_gehdr(fid);
hdr_rev = hdr.rdbm_rev;

% size of header (depends on revision)
% ==== Try different offset, by Zhe Liu ==== %
if hdr_rev >= 24
    hdr_size = 157276;      % Get from Pascal's code
elseif hdr_rev >= 20
    hdr_size = 145908+3880;
elseif hdr_rev >= 14
    hdr_size = 145908;
else
    disp('read_pfile: warning old pfile version, may not work');
    if hdr_rev < 7
        hdr_size = 39940;
    elseif hdr_rev < 8
        hdr_size = 39984;
    elseif hdr_rev < 9
        hdr_size = 60464;
    elseif hdr_rev < 11
        hdr_size = 61464;
    elseif hdr_rev < 14
        hdr_size = 66072;
    end
end

% some size parameters
nex  = hdr.navs
nex = 1
nsli = hdr.nslices
nech = hdr.nechoes
yres = hdr.nframes
yacq = hdr.da_yres; % includes initial baseline
xres = hdr.frame_size
edr  = hdr.point_size*2; % bytes per data point
nr   = double(hdr.rc_xres); % to avoid class mismatches
ncoi = hdr.dab_stop_rcv(1)+1

% storage
if edr==4
    precision = 'int16=>float';
    data = zeros(xres,yres,ncoi,nech,nsli,'single');
else
    precision = 'int32=>double';
    data = zeros(xres,yres,ncoi,nech,nsli,'double');
end

% read data
for s = 1:nsli
	for e = 1:nech

		% start at beginning of file, skip header
		fseek(fid,hdr_size,'bof');
        
		% jump to slice
		fseek(fid,(s-1)*xres*yacq*nech*edr,'cof');
        
		% jump to echo
		fseek(fid,(e-1)*xres*yacq*edr,'cof');
        
        % images from all coils
		for c = 1:ncoi

            % skip initial baseline
            fseek(fid,xres*edr,'cof');
            
			% read in data
            tmp = fread(fid,[2 double(xres*yres)],precision);
            
            % convert to complex
            tmp = complex(tmp(1,:),tmp(2,:));

            % save data
            data(:,:,c,e,s) = reshape(tmp,xres,yres);
            
			% next coil (-1 since fread moved us +1)
			fseek(fid,(nsli*nex*nech-1)*xres*yacq*edr,'cof');

        end
	end
end
fclose(fid);

% read trajectory and density (UTE files only)
fid = fopen([filename '.kspace'],'r');
if fid==-1
    kx = [];
    ky = [];
    d = [];
else
    traj = fread(fid,[3 double(xres*yres)],'float32','ieee-le');
    fclose(fid);
    kx = reshape(traj(1,:),xres,yres) * double(nr);
    ky = reshape(traj(2,:),xres,yres) * double(nr);
    d  = reshape(traj(3,:),xres,yres);
end

% ----------------------------------------------------------------- % 

function hdr = read_gehdr(id)

% Reads in a GE pfile header
%  id:    file identifier obtained from fopen()
%  hdr:   structure RDB_HEADER_REC from GE file rdbm.h
%
%  Script generated using 'gehdr2matlab' written by DB Clayton
%  Hacked by O Unal for 14.3F
%  Hacked by M Bydder for read_pfile()

fseek(id, 0, 'bof');  % go to start of raw data file
hdr.rdbm_rev = fread(id, 1, '*float');  
hdr.run_int = fread(id, 1, '*long');  % Rdy pkt Run Number 
hdr.scan_seq = fread(id, 1, 'short=>int');  % Rdy pkt Sequence Number 
hdr.run_char = fread(id, 6, '*char')';  % Rdy pkt Run no in char 
hdr.scan_date = fread(id, 10, '*char')';  % 
hdr.scan_time = fread(id, 8, '*char')';  % 
hdr.logo = fread(id, 10, '*char')';  % rdbm used to verify file 

hdr.file_contents = fread(id, 1, 'short=>int');  % Data type 0=emp 1=nrec 2=rw 0, 1, 2 
hdr.lock_mode = fread(id, 1, 'short=>int');  % unused 
hdr.dacq_ctrl = fread(id, 1, 'short=>int');  % rhdacqctrl bit mask 15 bits 
hdr.recon_ctrl = fread(id, 1, 'short=>int');  % rhrcctrl bit mask 15 bits 
hdr.exec_ctrl = fread(id, 1, 'short=>int');  % rhexecctrl bit mask 15 bits 
hdr.scan_type = fread(id, 1, 'short=>int');  % bit mask 15 bits 
hdr.data_collect_type = fread(id, 1, 'short=>int');  % rhtype bit mask 15 bits 
hdr.data_format = fread(id, 1, 'short=>int');  % rhformat bit mask 15 bits 
hdr.recon = fread(id, 1, 'short=>int');  % rhrecon proc-a-son recon 0 - 100 
hdr.datacq = fread(id, 1, 'short=>int');  % rhdatacq proc-a-son dacq 

hdr.npasses = fread(id, 1, 'short=>int');  % rhnpasses passes for a scan 0 - 256 
hdr.npomp = fread(id, 1, 'short=>int');  % rhnpomp pomp group slices 1,2 
hdr.nslices = fread(id, 1, 'short=>int');  % rhnslices slices in a pass 0 - 256 
hdr.nechoes = fread(id, 1, 'short=>int');  % rhnecho echoes of a slice 1 - 32 
hdr.navs = fread(id, 1, 'short=>int');  % rhnavs num of excitiations 1 - 32727 
hdr.nframes = fread(id, 1, 'short=>int');  % rhnframes yres 0 - 1024 
hdr.baseline_views = fread(id, 1, 'short=>int');  % rhbline baselines 0 - 1028 
hdr.hnover = fread(id, 1, 'short=>int');  % rhhnover overscans 0 - 1024 
hdr.frame_size = fread(id, 1, 'ushort=>int');  % rhfrsize xres 0 - 32768 
hdr.point_size = fread(id, 1, 'short=>int');  % rhptsize 2 - 4 

hdr.vquant = fread(id, 1, 'short=>int');  % rhvquant 3d volumes 1 

hdr.cheart = fread(id, 1, 'short=>int');  % RX Cine heart phases 1 - 32 
hdr.ctr = fread(id, 1, '*float');  % RX Cine TR in sec 0 - 3.40282e38
hdr.ctrr = fread(id, 1, '*float');  % RX Cine RR in sec 0 - 30.0 

hdr.initpass = fread(id, 1, 'short=>int');  % rhinitpass allocate passes 0 - 32767 
hdr.incrpass = fread(id, 1, 'short=>int');  % rhincrpass tps autopauses 0 - 32767 

hdr.method_ctrl = fread(id, 1, 'short=>int');  % rhmethod 0=recon, 1=psd 0, 1 
hdr.da_xres = fread(id, 1, 'ushort=>int');  % rhdaxres 0 - 32768 
hdr.da_yres = fread(id, 1, 'short=>int');  % rhdayres 0 - 2049 
hdr.rc_xres = fread(id, 1, 'short=>int');  % rhrcxres 0 - 1024 
hdr.rc_yres = fread(id, 1, 'short=>int');  % rhrcyres 0 - 1024 
hdr.im_size = fread(id, 1, 'short=>int');  % rhimsize 0 - 512 
hdr.rc_zres = fread(id, 1, '*long');  % power of 2 > rhnslices 0 - 128 

hdr.raw_pass_size = fread(id, 1, '*ulong');  % rhrawsize 0 - 2147483647
hdr.sspsave = fread(id, 1, '*ulong');  % rhsspsave 0 - 2147483647
hdr.udasave = fread(id, 1, '*ulong');  % rhudasave 0 - 2147483647

hdr.fermi_radius = fread(id, 1, '*float');  % rhfermr fermi radius 0 - 3.40282e38
hdr.fermi_width = fread(id, 1, '*float');  % rhfermw fermi width 0 - 3.40282e38
hdr.fermi_ecc = fread(id, 1, '*float');  % rhferme fermi excentiricty 0 - 3.40282e38
hdr.clip_min = fread(id, 1, '*float');  % rhclipmin 4x IP limit +-16383 
hdr.clip_max = fread(id, 1, '*float');  % rhclipmax 4x IP limit +-16383 
hdr.default_offset = fread(id, 1, '*float');  % rhdoffset default offset = 0 +-3.40282e38 
hdr.xoff = fread(id, 1, '*float');  % rhxoff scroll img in x +-256 
hdr.yoff = fread(id, 1, '*float');  % rhyoff scroll img in y +-256 
hdr.nwin = fread(id, 1, '*float');  % rhnwin hecho window width 0 - 256 
hdr.ntran = fread(id, 1, '*float');  % rhntran hecho trans width 0 - 256 
hdr.scalei = fread(id, 1, '*float');  % PS rhscalei +-3.40282e38 
hdr.scaleq = fread(id, 1, '*float');  % PS rhscaleq def = 0 +-3.40282e38 
hdr.rotation = fread(id, 1, 'short=>int');  % RX 0 90 180 270 deg 0 - 3 
hdr.transpose = fread(id, 1, 'short=>int');  % RX 0, 1 n / y transpose 0 - 1
hdr.kissoff_views = fread(id, 1, 'short=>int');  % rhblank zero image views 0 - 512 
hdr.slblank = fread(id, 1, 'short=>int');  % rhslblank slice blank 3d 0 - 128 
hdr.gradcoil = fread(id, 1, 'short=>int');  % RX 0=off 1=Schnk 2=Rmr 0 - 2 
hdr.ddaover = fread(id, 1, 'short=>int');  % rhddaover unused 

hdr.sarr = fread(id, 1, 'short=>int');  % SARR bit mask 15 bits 
hdr.fd_tr = fread(id, 1, 'short=>int');  % SARR feeder timing info 
hdr.fd_te = fread(id, 1, 'short=>int');  % SARR feeder timing info 
hdr.fd_ctrl = fread(id, 1, 'short=>int');  % SARR control of feeder 
hdr.algor_num = fread(id, 1, 'short=>int');  % SARR df decimation ratio 
hdr.fd_df_dec = fread(id, 1, 'short=>int');  % SARR which feeder algor 

buff = fread(id, 8, 'short=>int');  % kluge for type RDB_MULTI_RCV_TYPE
hdr.dab_start_rcv = buff(1:2:end);  % kluge for type RDB_MULTI_RCV_TYPE
hdr.dab_stop_rcv = buff(2:2:end);  % kluge for type RDB_MULTI_RCV_TYPE

hdr.user0 = fread(id, 1, '*float');  % rhuser0 +-3.40282e38 
hdr.user1 = fread(id, 1, '*float');  % rhuser1 +-3.40282e38 
hdr.user2 = fread(id, 1, '*float');  % rhuser2 +-3.40282e38 
hdr.user3 = fread(id, 1, '*float');  % rhuser3 +-3.40282e38 
hdr.user4 = fread(id, 1, '*float');  % rhuser4 +-3.40282e38 
hdr.user5 = fread(id, 1, '*float');  % rhuser5 +-3.40282e38 
hdr.user6 = fread(id, 1, '*float');  % rhuser6 +-3.40282e38 
hdr.user7 = fread(id, 1, '*float');  % rhuser7 +-3.40282e38 
hdr.user8 = fread(id, 1, '*float');  % rhuser8 +-3.40282e38 
hdr.user9 = fread(id, 1, '*float');  % rhuser9 +-3.40282e38 
hdr.user10 = fread(id, 1, '*float');  % rhuser10 +-3.40282e38 
hdr.user11 = fread(id, 1, '*float');  % rhuser11 +-3.40282e38 
hdr.user12 = fread(id, 1, '*float');  % rhuser12 +-3.40282e38 
hdr.user13 = fread(id, 1, '*float');  % rhuser13 +-3.40282e38 
hdr.user14 = fread(id, 1, '*float');  % rhuser14 +-3.40282e38 
hdr.user15 = fread(id, 1, '*float');  % rhuser15 +-3.40282e38 
hdr.user16 = fread(id, 1, '*float');  % rhuser16 +-3.40282e38 
hdr.user17 = fread(id, 1, '*float');  % rhuser17 +-3.40282e38 
hdr.user18 = fread(id, 1, '*float');  % rhuser18 +-3.40282e38 
hdr.user19 = fread(id, 1, '*float');  % rhuser19 +-3.40282e38 

hdr.v_type = fread(id, 1, '*long');  % rhvtype bit mask 31 bits 
hdr.v_coefxa = fread(id, 1, '*float');  % RX x flow direction control 0 - 4 
hdr.v_coefxb = fread(id, 1, '*float');  % RX x flow direction control 0 - 4 
hdr.v_coefxc = fread(id, 1, '*float');  % RX x flow direction control 0 - 4 
hdr.v_coefxd = fread(id, 1, '*float');  % RX x flow direction control 0 - 4 
hdr.v_coefya = fread(id, 1, '*float');  % RX y flow direction control 0 - 4 
hdr.v_coefyb = fread(id, 1, '*float');  % RX y flow direction control 0 - 4 
hdr.v_coefyc = fread(id, 1, '*float');  % RX y flow direction control 0 - 4 
hdr.v_coefyd = fread(id, 1, '*float');  % RX y flow direction control 0 - 4 
hdr.v_coefza = fread(id, 1, '*float');  % RX z flow direction control 0 - 4 
hdr.v_coefzb = fread(id, 1, '*float');  % RX z flow direction control 0 - 4 
hdr.v_coefzc = fread(id, 1, '*float');  % RX z flow direction control 0 - 4 
hdr.v_coefzd = fread(id, 1, '*float');  % RX z flow direction control 0 - 4 
hdr.vm_coef1 = fread(id, 1, '*float');  % RX weight for mag image 1 0 - 1 
hdr.vm_coef2 = fread(id, 1, '*float');  % RX weight for mag image 2 0 - 1 
hdr.vm_coef3 = fread(id, 1, '*float');  % RX weight for mag image 3 0 - 1 
hdr.vm_coef4 = fread(id, 1, '*float');  % RX weight for mag image 4 0 - 1 
hdr.v_venc = fread(id, 1, '*float');  % RX vel encodeing cm / sec 0.001 - 5000 

hdr.spectral_width = fread(id, 1, '*float');  % specwidth filter width kHz 500 - 3355432 
hdr.csi_dims = fread(id, 1, 'short=>int');  % spectro 
hdr.xcsi = fread(id, 1, 'short=>int');  % rhspecrescsix 2 - 64 
hdr.ycsi = fread(id, 1, 'short=>int');  % rhspecrescsiy 2 - 64 
hdr.zcsi = fread(id, 1, 'short=>int');  % spectro 
hdr.roilenx = fread(id, 1, '*float');  % RX x csi volume dimension 
hdr.roileny = fread(id, 1, '*float');  % RX y csi volume dimension 
hdr.roilenz = fread(id, 1, '*float');  % RX z csi volume dimension 
hdr.roilocx = fread(id, 1, '*float');  % RX x csi volume center 
hdr.roilocy = fread(id, 1, '*float');  % RX y csi volume center 
hdr.roilocz = fread(id, 1, '*float');  % RX z csi volume center 
hdr.numdwell = fread(id, 1, '*float');  % specdwells 0 - 3.40282e38

hdr.ps_command = fread(id, 1, '*long');  % PS internal use only 
hdr.ps_mps_r1 = fread(id, 1, '*long');  % PS MPS R1 setting 1 - 7 
hdr.ps_mps_r2 = fread(id, 1, '*long');  % PS MPS R2 setting 1 - 30 
hdr.ps_mps_tg = fread(id, 1, '*long');  % PS MPS Transmit gain setting 0 - 200
hdr.ps_mps_freq = fread(id, 1, '*long');  % PS MPS Center frequency hz +-3.40282e38 
hdr.ps_aps_r1 = fread(id, 1, '*long');  % PS APS R1 setting 1 - 7 
hdr.ps_aps_r2 = fread(id, 1, '*long');  % PS APS R2 setting 1 - 30 
hdr.ps_aps_tg = fread(id, 1, '*long');  % PS APS Transmit gain setting 0 - 200
hdr.ps_aps_freq = fread(id, 1, '*long');  % PS APS Center frequency hz +-3.40282e38 
hdr.ps_scalei = fread(id, 1, '*float');  % PS rational scaling +-3.40282e38 
hdr.ps_scaleq = fread(id, 1, '*float');  % PS unused 
hdr.ps_snr_warning = fread(id, 1, '*long');  % PS noise test 0=16 1=32 bits 0, 1 
hdr.ps_aps_or_mps = fread(id, 1, '*long');  % PS prescan order logic 0 - 5 
hdr.ps_mps_bitmap = fread(id, 1, '*long');  % PS bit mask 4 bits
hdr.ps_powerspec = fread(id, 256, '*char')';  % PS 
hdr.ps_filler1 = fread(id, 1, '*long');  % PS filler 
hdr.ps_filler2 = fread(id, 1, '*long');  % PS filler 
hdr.obsolete1 = fread(id, 16, '*float');  % PS mean noise each receiver +-3.40282e38 
hdr.obsolete2 = fread(id, 16, '*float');  % PS noise calc for muti rec +-3.40282e38 

hdr.halfecho = fread(id, 1, 'short=>int');  % spectro full, half echo 0, 1 

hdr.im_size_y = fread(id, 1, 'short=>int');  % rh???? 0 - 512 
hdr.data_collect_type1 = fread(id, 1, '*long');  % rh???? bit mask 31 bits 
hdr.freq_scale = fread(id, 1, '*float');  % rh???? freq k-space step +-3.40282e38 
hdr.phase_scale = fread(id, 1, '*float');  % rh???? freq k-space step +-3.40282e38 

hdr.ovl = fread(id, 1, 'short=>int');  % rhovl - overlaps for MOTSA 
hdr.pclin = fread(id, 1, 'short=>int');  % Linear Corr. 0:off, 1:linear, 2:polynomial 
hdr.pclinnpts = fread(id, 1, 'short=>int');  % fit number of points 
hdr.pclinorder = fread(id, 1, 'short=>int');  % fit order 
hdr.pclinavg = fread(id, 1, 'short=>int');  % linear phase corr avg 0:off, 1:on 
hdr.pccon = fread(id, 1, 'short=>int');  % Const Corr. 0:off, 1:Ky spec., 2:polyfit(2/ilv), 3:polyfit(1/ilv) 
hdr.pcconnpts = fread(id, 1, 'short=>int');  % fit number of points 
hdr.pcconorder = fread(id, 1, 'short=>int');  % fit order 
hdr.pcextcorr = fread(id, 1, 'short=>int');  % external correction file 0:don't use, 1: use 
hdr.pcgraph = fread(id, 1, 'short=>int');  % Phase Correction coef. image 0:off, 1:linear & constant 
hdr.pcileave = fread(id, 1, 'short=>int');  % Interleaves to use for correction: 0=all, 1=only first 
hdr.hdbestky = fread(id, 1, 'short=>int');  % bestky view for fractional Ky scan 
hdr.pcctrl = fread(id, 1, 'short=>int');  % phase correction research control 
hdr.pcthrespts = fread(id, 1, 'short=>int');  % 2..512 adjacent points 
hdr.pcdiscbeg = fread(id, 1, 'short=>int');  % 0..512 beginning point to discard 
hdr.pcdiscmid = fread(id, 1, 'short=>int');  % 0..512 middle point to discard 
hdr.pcdiscend = fread(id, 1, 'short=>int');  % 0..512 ending point to discard 
hdr.pcthrespct = fread(id, 1, 'short=>int');  % Threshold percentage 
hdr.pcspacial = fread(id, 1, 'short=>int');  % Spacial best ref scan index 0..512 
hdr.pctemporal = fread(id, 1, 'short=>int');  % Temporal best ref scan index 0..512 
hdr.pcspare = fread(id, 1, 'short=>int');  % spare for phase correction 
hdr.ileaves = fread(id, 1, 'short=>int');  % Number of interleaves 
hdr.kydir = fread(id, 1, 'short=>int');  % Ky traversal dircetion 0: top-down, 1:center out 
hdr.alt = fread(id, 1, 'short=>int');  % Alt read sign 0=no, 1=odd/even, 2=pairs 
hdr.reps = fread(id, 1, 'short=>int');  % Number of scan repetitions 
hdr.ref = fread(id, 1, 'short=>int');  % Ref Scan 0: off 1: on 

hdr.pcconnorm = fread(id, 1, '*float');  % Constant S term normalization factor 
hdr.pcconfitwt = fread(id, 1, '*float');  % Constant polyfit weighting factor 
hdr.pclinnorm = fread(id, 1, '*float');  % Linear S term normalization factor 
hdr.pclinfitwt = fread(id, 1, '*float');  % Linear polyfit weighting factor 

hdr.pcbestky = fread(id, 1, '*float');  % Best Ky location 

hdr.vrgf = fread(id, 1, '*long');  % control word for VRG filter 
hdr.vrgfxres = fread(id, 1, '*long');  % control word for VRGF final x resolution 

hdr.bp_corr = fread(id, 1, '*long');  % control word for bandpass asymmetry 
hdr.recv_freq_s = fread(id, 1, '*float');  % starting frequency (+62.5) 
hdr.recv_freq_e = fread(id, 1, '*float');  % ending frequency (-62.5) 

hdr.hniter = fread(id, 1, '*long');  % Selects the number of (continued...)

hdr.fast_rec = fread(id, 1, '*long');  % Added for homodyne II, tells if (continued...)

hdr.refframes = fread(id, 1, '*long');  % total # of frames for ref scan 
hdr.refframep = fread(id, 1, '*long');  % # of frames per pass for a ref scan 
hdr.scnframe = fread(id, 1, '*long');  % total # of frames for a entire scan 
hdr.pasframe = fread(id, 1, '*long');  % # of frames per pass 

hdr.user_usage_tag = fread(id, 1, '*ulong');  % for spectro 
hdr.user_fill_mapMSW = fread(id, 1, '*ulong');  % for spectro 
hdr.user_fill_mapLSW = fread(id, 1, '*ulong');  % for Spectro 

hdr.user20 = fread(id, 1, '*float');  % all following usercv are for spectro 
hdr.user21 = fread(id, 1, '*float');  
hdr.user22 = fread(id, 1, '*float');  
hdr.user23 = fread(id, 1, '*float');  
hdr.user24 = fread(id, 1, '*float');  
hdr.user25 = fread(id, 1, '*float');  
hdr.user26 = fread(id, 1, '*float');  
hdr.user27 = fread(id, 1, '*float');  
hdr.user28 = fread(id, 1, '*float');  
hdr.user29 = fread(id, 1, '*float');  
hdr.user30 = fread(id, 1, '*float');  
hdr.user31 = fread(id, 1, '*float');  
hdr.user32 = fread(id, 1, '*float');  
hdr.user33 = fread(id, 1, '*float');  
hdr.user34 = fread(id, 1, '*float');  
hdr.user35 = fread(id, 1, '*float');  
hdr.user36 = fread(id, 1, '*float');  
hdr.user37 = fread(id, 1, '*float');  
hdr.user38 = fread(id, 1, '*float');  
hdr.user39 = fread(id, 1, '*float');  
hdr.user40 = fread(id, 1, '*float');  
hdr.user41 = fread(id, 1, '*float');  
hdr.user42 = fread(id, 1, '*float');  
hdr.user43 = fread(id, 1, '*float');  
hdr.user44 = fread(id, 1, '*float');  
hdr.user45 = fread(id, 1, '*float');  
hdr.user46 = fread(id, 1, '*float');  
hdr.user47 = fread(id, 1, '*float');  
hdr.user48 = fread(id, 1, '*float');  

hdr.pcfitorig = fread(id, 1, 'short=>int');  % Adjust view indexes if set so bestky view = 0 
hdr.pcshotfirst = fread(id, 1, 'short=>int');  % First view within an echo group used for fit 
hdr.pcshotlast = fread(id, 1, 'short=>int');  % Last view within an echo group used for fit 
hdr.pcmultegrp = fread(id, 1, 'short=>int');  % If = 1, force pts from other egrps to be used 
hdr.pclinfix = fread(id, 1, 'short=>int');  % If = 2, force slope to be set to pclinslope 

hdr.pcconfix = fread(id, 1, 'short=>int');  % If = 2, force slope to be set to pcconslope 

hdr.pclinslope = fread(id, 1, '*float');  % Value to set lin slope to if forced 
hdr.pcconslope = fread(id, 1, '*float');  % Value to set con slope to if forced 
hdr.pccoil = fread(id, 1, 'short=>int');  % If 1,2,3,4, use that coil's results for all 

hdr.vvsmode = fread(id, 1, 'short=>int');  % Variable view sharing mode 
hdr.vvsaimgs = fread(id, 1, 'short=>int');  % number of original images 
hdr.vvstr = fread(id, 1, 'short=>int');  % TR in microseconds 
hdr.vvsgender = fread(id, 1, 'short=>int');  % gender: male or female 

hdr.zip_factor = fread(id, 1, 'short=>int');  % Slice ZIP factor: 0=OFF, 2, or 4 

hdr.maxcoef1a = fread(id, 1, '*float');  % Coefficient A for flow image 1 
hdr.maxcoef1b = fread(id, 1, '*float');  % Coefficient B for flow image 1 
hdr.maxcoef1c = fread(id, 1, '*float');  % Coefficient C for flow image 1 
hdr.maxcoef1d = fread(id, 1, '*float');  % Coefficient D for flow image 1 
hdr.maxcoef2a = fread(id, 1, '*float');  % Coefficient A for flow image 2 
hdr.maxcoef2b = fread(id, 1, '*float');  % Coefficient B for flow image 2 
hdr.maxcoef2c = fread(id, 1, '*float');  % Coefficient C for flow image 2 
hdr.maxcoef2d = fread(id, 1, '*float');  % Coefficient D for flow image 2 
hdr.maxcoef3a = fread(id, 1, '*float');  % Coefficient A for flow image 3 
hdr.maxcoef3b = fread(id, 1, '*float');  % Coefficient B for flow image 3 
hdr.maxcoef3c = fread(id, 1, '*float');  % Coefficient C for flow image 3 
hdr.maxcoef3d = fread(id, 1, '*float');  % Coefficient D for flow image 3 

hdr.ut_ctrl = fread(id, 1, '*long');  % System utility control variable 
hdr.dp_type = fread(id, 1, 'short=>int');  % EPI II diffusion control cv 

hdr.arw = fread(id, 1, 'short=>int');  % Arrhythmia rejection window(percentage:1-100)

hdr.vps = fread(id, 1, 'short=>int');  % View Per Segment for FastCine 

hdr.mcReconEnable = fread(id, 1, 'short=>int');  % N-Coil recon map 
hdr.fov = fread(id, 1, '*float');  % Auto-NCoil 

hdr.te = fread(id, 1, '*long');  % TE for first echo 
hdr.te2 = fread(id, 1, '*long');  % TE for second and later echoes 
hdr.dfmrbw = fread(id, 1, '*float');  % BW for navigator frames 
hdr.dfmctrl = fread(id, 1, '*long');  % Control flag for dfm (0=off, other=on)
hdr.raw_nex = fread(id, 1, '*long');  % Uncombined NEX at start of recon 
hdr.navs_per_pass = fread(id, 1, '*long');  % Max. navigator frames in a pass 
hdr.dfmxres = fread(id, 1, '*long');  % xres of navigator frames 
hdr.dfmptsize = fread(id, 1, '*long');  % point size of navigator frames 
hdr.navs_per_view = fread(id, 1, '*long');  % Num. navigators per frame (tag table) 
hdr.dfmdebug = fread(id, 1, '*long');  % control flag for dfm debug 
hdr.dfmthreshold = fread(id, 1, '*float');  % threshold for navigator correction 

hdr.grid_control = fread(id, 1, 'short=>int');  % bit settings controlling gridding 
hdr.b0map = fread(id, 1, 'short=>int');  % B0 map enable and map size 
hdr.grid_tediff = fread(id, 1, 'short=>int');  % TE difference between b0 map arms 
hdr.grid_motion_comp = fread(id, 1, 'short=>int');  % flag to apply motion compensation 
hdr.grid_radius_a = fread(id, 1, '*float');  % variable density transition 
hdr.grid_radius_b = fread(id, 1, '*float');  % variable density transition 
hdr.grid_max_gradient = fread(id, 1, '*float');  % Max gradient amplitude 
hdr.grid_max_slew = fread(id, 1, '*float');  % Max slew rate 
hdr.grid_scan_fov = fread(id, 1, '*float');  % Rx scan field of view 
hdr.grid_a2d_time = fread(id, 1, '*float');  % A to D sample time microsecs 
hdr.grid_density_factor = fread(id, 1, '*float');  % change factor for variable density 
hdr.grid_display_fov = fread(id, 1, '*float');  % Rx display field of view 

hdr.fatwater = fread(id, 1, 'short=>int');  % for Fat and Water Dual Recon 
hdr.fiestamlf = fread(id, 1, 'short=>int');  % MFO FIESTA recon control bit 16bits 

hdr.app = fread(id, 1, 'short=>int');  % Auto Post-Processing opcode 
hdr.rhncoilsel = fread(id, 1, 'short=>int');  % Auto-Ncoil 
hdr.rhncoillimit = fread(id, 1, 'short=>int');  % Auto-Ncoil 
hdr.app_option = fread(id, 1, 'short=>int');  % Auto Post_processing options 
hdr.grad_mode = fread(id, 1, 'short=>int');  % Gradient mode in Gemini project 
hdr.pfile_passes = fread(id, 1, 'short=>int');  % Num passes stored in a multi-pass Pfile (0 means 1 pass) 

hdr.asset = fread(id, 1, '*int');  
hdr.asset_calthresh = fread(id, 1, '*int');  
hdr.asset_R = fread(id, 1, '*float');  
hdr.coilno = fread(id, 1, '*int');  
hdr.asset_phases = fread(id, 1, '*int');  
hdr.scancent = fread(id, 1, '*float');  % Table position 
hdr.position = fread(id, 1, '*int');  % Patient position 
hdr.entry = fread(id, 1, '*int');  % Patient entry 
hdr.lmhor = fread(id, 1, '*float');  % Landmark 
hdr.last_slice_num = fread(id, 1, '*int');   
hdr.asset_slice_R = fread(id, 1, '*float');  % Slice reduction factor 
hdr.asset_slabwrap = fread(id, 1, '*float');  

hdr.dwnav_coeff = fread(id, 1, '*float');  % Coeff for amount of phase correction 
hdr.dwnav_cor = fread(id, 1, 'short=>int');  % Navigator echo correction 
hdr.dwnav_view = fread(id, 1, 'short=>int');  % Num of views of nav echoes 
hdr.dwnav_corecho = fread(id, 1, 'short=>int');  % Num of nav echoes for actual correction 
hdr.dwnav_sview = fread(id, 1, 'short=>int');  % Start view for phase correction process 
hdr.dwnav_eview = fread(id, 1, 'short=>int');  % End view for phase correction process 
hdr.dwnav_sshot = fread(id, 1, 'short=>int');  % Start shot for delta phase estimation in nav echoes 
hdr.dwnav_eshot = fread(id, 1, 'short=>int');  % End shot for delta phase estimation in nav echoes 

hdr.win3d_type = fread(id, 1, 'short=>int');  % 0 = Modified Hanning, 1 = modified Tukey 
hdr.win3d_apod = fread(id, 1, '*float');  % degree of apodization; 0.0 = boxcar, 1.0=hanning 
hdr.win3d_q = fread(id, 1, '*float');  % apodization at ends, 0.0 = max, 1.0 = boxcar 

hdr.ime_scic_enable = fread(id, 1, 'short=>int');  % Surface Coil Intensity Correction: 1 if enabled 
hdr.clariview_type = fread(id, 1, 'short=>int');  % Type of Clariview/Name of Filter 
hdr.ime_scic_edge = fread(id, 1, '*float');  % Edge paramaters for Enhanced Recon 
hdr.ime_scic_smooth = fread(id, 1, '*float');  % Smooth paramaters for Enhanced Recon 
hdr.ime_scic_focus = fread(id, 1, '*float');  % Focus paramaters for Enhanced Recon 
hdr.clariview_edge = fread(id, 1, '*float');  % Edge paramaters for clariview 
hdr.clariview_smooth = fread(id, 1, '*float');  % Smooth paramaters for clariview 
hdr.clariview_focus = fread(id, 1, '*float');  % Focus paramaters for clariview 
hdr.scic_reduction = fread(id, 1, '*float');  % Reduction paramater for SCIC 
hdr.scic_gauss = fread(id, 1, '*float');  % Gauss paramater for SCIC 
hdr.scic_threshold = fread(id, 1, '*float');  % Threshold paramater for SCIC 

hdr.ectricks_no_regions = fread(id, 1, '*long');  % Total no of regions acquired by PSD 
hdr.ectricks_input_regions = fread(id, 1, '*long');  % Total no of input regions for reordering 

hdr.psc_reuse = fread(id, 1, 'short=>int');  % Header field for smart prescan 

hdr.left_blank = fread(id, 1, 'short=>int');  
hdr.right_blank = fread(id, 1, 'short=>int');  

hdr.acquire_type = fread(id, 1, 'short=>int');  % Acquire type information from CV 

hdr.retro_control = fread(id, 1, 'short=>int');  % Retrosective FSE phase correction control flag. (continued...)

hdr.etl = fread(id, 1, 'short=>int');  % Added for Retrospective FSE phase correction. This (continued...)

hdr.pcref_start = fread(id, 1, 'short=>int');  % 1st view to use for dynamic EPI phase correction. 
hdr.pcref_stop = fread(id, 1, 'short=>int');  % Last view to use for dynamic EPI phase correction. 
hdr.ref_skip = fread(id, 1, 'short=>int');  % Number of passes to skip for dynamic EPI phase correction. 
hdr.extra_frames_top = fread(id, 1, 'short=>int');  % Number of extra frames at top of K-space 
hdr.extra_frames_bot = fread(id, 1, 'short=>int');  % Number of extra frames at bottom of K-space 
hdr.multiphase_type = fread(id, 1, 'short=>int');  % 0 = INTERLEAVED , 1 = SEQUENTIAL 
hdr.nphases = fread(id, 1, 'short=>int');  % Number of phases in a multiphase scan 
hdr.pure = fread(id, 1, 'short=>int');  % PURE flag from psd 
hdr.pure_scale = fread(id, 1, '*float');  % Recon scale factor ratio for cal scan 
hdr.off_data = fread(id, 1, '*int');  % Byte offset to start of raw data (i.e size of POOL_HEADER) 
hdr.off_per_pass = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_per_pass of POOL_HEADER 
hdr.off_unlock_raw = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_unlock_raw of POOL_HEADER 
hdr.off_data_acq_tab = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_data_acq_tab of POOL_HEADER 
hdr.off_nex_tab = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_nex_tab of POOL_HEADER 
hdr.off_nex_abort_tab = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_nex_abort_tab of POOL_HEADER 
hdr.off_tool = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_tool of POOL_HEADER 
hdr.off_exam = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_exam of POOL_HEADER 
hdr.off_series = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_series of POOL_HEADER 
hdr.off_image = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_image of POOL_HEADER 
hdr.off_ps = fread(id, 1, '*int');  % Byte offset to start of rdb_hdr_ps of POOL_HEADER 
hdr.off_spare_b = fread(id, 1, '*int');  % spare 
hdr.new_wnd_level_flag = fread(id, 1, '*int');  % New WW/WL algo enable/disable flag 
hdr.wnd_image_hist_area = fread(id, 1, '*int');  % Image Area % 
hdr.wnd_high_hist = fread(id, 1, '*float');  % Histogram Area Top 
hdr.wnd_lower_hist = fread(id, 1, '*float');  % Histogram Area Bottom 
hdr.pure_filter = fread(id, 1, 'short=>int');  % PURE noise reduction on=1/off=0 
hdr.cfg_pure_filter = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_fit_order = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_kernelsize_z = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_kernelsize_xy = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_weight_radius = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_intensity_scale = fread(id, 1, 'short=>int');  % PURE cfg file value 
hdr.cfg_pure_noise_threshold = fread(id, 1, 'short=>int');  % PURE cfg file value 

hdr.wienera = fread(id, 1, '*float');  % NB maintain alignment of floats 
hdr.wienerb = fread(id, 1, '*float');  
hdr.wienert2 = fread(id, 1, '*float');  
hdr.wieneresp = fread(id, 1, '*float');  
hdr.wiener = fread(id, 1, 'short=>int');  
hdr.flipfilter = fread(id, 1, 'short=>int');  
hdr.dbgrecon = fread(id, 1, 'short=>int');  
hdr.ech2skip = fread(id, 1, 'short=>int');  

hdr.tricks_type = fread(id, 1, '*int');  % 0 = Subtracted, 1 = Unsubtracted 

hdr.lcfiesta_phase = fread(id, 1, '*float');  % LC Fiesta 
hdr.lcfiesta = fread(id, 1, 'short=>int');  % LC Fiesta 
hdr.herawflt = fread(id, 1, 'short=>int');  % Half echo raw data filter 
hdr.herawflt_befnwin = fread(id, 1, 'short=>int');  % Half echo raw data filter 
hdr.herawflt_befntran = fread(id, 1, 'short=>int');  % Half echo raw data filter 
hdr.herawflt_befamp = fread(id, 1, '*float');  % Half echo raw data filter 
hdr.herawflt_hpfamp = fread(id, 1, '*float');  % Half echo raw data filter 
hdr.heover = fread(id, 1, 'short=>int');  % Half echo over sampling 

hdr.pure_correction_threshold = fread(id, 1, 'short=>int');  % PURE Correction threshold 

hdr.swiftenable = fread(id, 1, '*int');  % SWIFT enable/disable flag 
hdr.numslabs = fread(id, 1, 'short=>int');  % Number of slabs to be used by TRICKS 
hdr.swiftcoilnos = fread(id, 1, 'ushort=>int');  % Number of coils to SWIFT between 

hdr.ps_autoshim_status = fread(id, 1, '*int');   

hdr.dynaplan_numphases = fread(id, 1, '*int');  % Number of phases for Dynamic Plan 

hdr.excess = fread(id, 216, 'short=>int');  % free space for later expansion 
