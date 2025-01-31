function [] = KDecon_v4(varargin)

    try
        disp(' ');
        disp('KDecon v4 - KINTSUGI'); 
        disp(' ');
        % Check the number of input arguments
        if nargin ~= 16
            showinfo();
            if isdeployed
                exit(1);
            end
            return
        end

        % Read command line parameters
        inpath = varargin{1};
        dxy = varargin{2};
        dz = varargin{3};
        numit = varargin{4};
        NA = varargin{5};
        rf = varargin{6};
        lambda_ex = varargin{7};
        lambda_em = varargin{8};
        fcyl = varargin{9};
        slitwidth = varargin{10};
        damping = varargin{11};
        clipval = varargin{12};
        stop_criterion = varargin{13};
        mem_percent = varargin{14};
        device = varargin{15};
        blocksize = varargin{16};

        %convert command line parameters from string to double
        if isdeployed
            dxy = str2double(dxy);      
            dz = str2double(dz);
            numit = str2double(numit);
            NA = str2double(NA);
            rf = str2double(rf);
            lambda_ex = str2double(lambda_ex);
            lambda_em = str2double(lambda_em);
            fcyl = str2double(fcyl);
            slitwidth = str2double(slitwidth);
            damping = str2double(damping);
            clipval = str2double(clipval);
            stop_criterion = str2double(stop_criterion);
            mem_percent = str2double(mem_percent);
            blocksize = str2double(blocksize);
        end

        %start data processing
        [CPUmem, GPUmem, deviceNames] = getmemory;

        if strcmpi(device, 'CPU')
            mem = CPUmem * mem_percent;
            deviceUsed = 'CPU';
        else
            mem = GPUmem * mem_percent;
            deviceUsed = 'GPU';
        end
        if isnan(mem) || isinf(mem) || mem <= 0
            error('Invalid memory calculation. Check memory percent value.');
        end
        fprintf('Available devices and their memory:\n');
        for i = 1:length(deviceNames)
            sprintf('%s\n', deviceNames{i});
        end
        fprintf('Using %s with %.2f GB of memory\n', deviceUsed, mem / 1e9);

        [nx, ny, nz] = getstackinfo(inpath);
        if nx * ny * nz == 0
            error('No valid TIFF files could be found!');
        end

        [tx, ty, tz] = autosplit(nx, ny, nz, mem, deviceUsed, blocksize);
        if tx * ty * tz > 1
            disp(['deconvolution split into ' num2str(tx) ' x ' num2str(ty) ' x ' num2str(tz) ' blocks.']);
        end    

        process(inpath, tx, ty, tz, dxy, dz, numit, NA, rf, lambda_ex, lambda_em, fcyl, slitwidth, damping, clipval, stop_criterion, deviceUsed);
        if isdeployed
            exit(0);
        end
    catch ME 
        % Error handling
        text = getReport(ME, 'extended', 'hyperlinks', 'off');
        disp(text);
        if isdeployed
            exit(1);
        end
    end 
end

function [CPU_mem, GPU_mem, deviceNames] = getmemory
    deviceNames = {};
    if ispc
        [~, m] = memory;
        CPU_mem = m.PhysicalMemory.Available;
        fprintf('CPU supported with %.2f GB memory\n', CPU_mem / 1e9);
        deviceNames{end+1} = fprintf('CPU: %.2f GB available\n', CPU_mem / 1e9);
    else
        [~, result] = system('vmstat -s -S M | grep "free memory"');
        CPU_mem = str2double(regexp(result, '\d+', 'match')) * 1e6;
        deviceNames{end+1} = fprintf('CPU: %.2f GB available', CPU_mem / 1e9);
    end

    try
        check = gpuDevice; %check if CUDA is available
        if (check.DeviceSupported == 1)
            g = gpuDevice;
            GPU_mem = g.AvailableMemory;
            fprintf('GPU supported with %.2f GB memory\n', GPU_mem / 1e9);
            deviceNames{end+1} = fprintf('GPU: %s\n', g.Name);
        else
            GPU_mem = 0.0;
        end
    catch
        GPU_mem = 0.0;
        deviceNames{end+1} = 'No supported GPU found';
    end
end

function showinfo()    
    disp('Usage: LsDeconv TIFDIR DELTAXY DELTAZ nITER NA RI LAMBDA_EX LAMBDA_EM FCYL SLITWIDTH DAMPING HISOCLIP STOP_CRIT MEM_PERCENT');
    disp(' ');
    disp('TIFFDIR: Directory containing the 2D-tiff files to be deconvolved(16-bit 0r 32bit float grayscale images supported).');
    disp('The images are expected to be formated with numerical endings, as e.g. xx00001.tif, xx00002.tif, ... xx00010.tif....');    
    disp(' ');
    disp('DELTAXY DELTAZ: xy- and z-Size of a voxel in nanometer. Choosing e.g. 250 500 means that a voxel is 250 nm x 250 nm wide in x- and y-direction');
    disp('and 500 nm wide in z-direction (vertical direction). Values depend on the total magnification of the microscope and the camera chip.');
    disp(' ');
    disp('nITER: max. number of iterations for Lucy-Richardson algorithm. Deconvolution stops before if stop_crit is reached.');
    disp(' ');
    disp('NA: numerical aperture of the objective.');
    disp(' ');
    disp('RI: refractive index of imaging medium and sample');
    disp(' ');
    disp('LAMBDA_EX: fluorescence excitation wavelength in nanometer.');
    disp(' ');
    disp('LAMBDA_EM: fluorescence emmision wavelength in nanometer.');
    disp(' ');
    disp('FCYL: focal length f in millimeter of the cylinder lens used for light sheet generation.');
    disp(' ');
    disp('SLITWIDTH: full width w in millimeter of the slit aperture placed in front of the cylinder lens. The NA of the light sheet');
    disp('generator system is calculated as NaLs = sin(arctan(w/(2*f))');
    disp(' ');
    disp('DAMPING: parameter between 0% and 10%. Increase value for images that are noisy. For images with');
    disp('good signal to noise ratio the damping parameter should be set to zero (no damping at all)'); 
    disp(' ');
    disp('HISTOCLIP: percent value between 0 and 5 percent. If HISTOCLIP is set e.g. to 0.01% then the histogram of the deconvolved');
    disp('stack is clipped at the 0.01 and the 99.99 percentile and the remaininginte intensity values are scaled to the full range (0...65535');
    disp('in case of of 16 bit images, and 0...Imax in case of 32 bit float images, where Imax is the highest intensity value');
    disp('occuring in the source stack');
    disp(' ');
    disp('STOP_CRIT: I the pecentual change to the last iteration step becomes lower than STOP_CRIT the deconvolution of the current');
    disp('block is finished. If STOP_CRIT is e.g. set to 2% then the iteration stops, if there is less than 2% percent');
    disp('change compared to the last iteration step.');
    disp(' ');
    disp('MEM_PERCENT: percent of RAM (or GPU memory, respectively, that can maximally by occopied by a data block. If the size of the image stack');
    disp('is larger than RAMSIZE  *  MEM_PERCENT / 100, the data set is split into blocks that are deconvolved sequentially and then stitched.');
    disp('A value of 4% usually is a good choice when working on CPU, a value of 50% when using the GPU. Decrease this value if other memory consuming');
	disp('programs are active.');    
    disp(' ');
    disp("device: 'CPU' = peform convolutions on CPU, 'GPU' = perform convolutions on GPU");
end

%determine the required number of blocks that are deconvolved sequentially
function [tx, ty, tz] = autosplit(npix_x, npix_y, npix_z, maxblocksize, deviceUsed, blocksize)  
    
    % Input validation
    if npix_x <= 0 || npix_y <= 0 || npix_z <= 0
        error('All image dimensions must be positive');
    end
    if maxblocksize <= 0
        error('Memory allocation must be positive');
    end
    disp(['Image dimensions: x=' num2str(npix_x) ' y=' num2str(npix_y) ' z=' num2str(npix_z) ' blocksize' num2str(blocksize)]);
    
    tx = 1; ty = 1; tz = 1;

    bytes = npix_x * npix_y * npix_z * 4;

    fprintf('Total pixels: %.2f GB Max block size: %.2f GB\n', bytes / 1e9, (maxblocksize / 1e9)*0.8);
    
    useGPU = strcmpi(deviceUsed, 'GPU');

    while true
        % Calculate current block dimensions
        block_x = ceil(npix_x / tx);
        block_y = ceil(npix_y / ty);
        
        if useGPU
            block_x = ceil(block_x/blocksize);
            block_y = ceil(block_y/blocksize);
       
            block_bytes = (block_x * blocksize) * (block_y* blocksize) * npix_z * 4;
            if block_bytes <= (maxblocksize*0.8) && blocksize > 1
                disp(['block_x ' num2str(block_x) ' block_y ' num2str(block_y) ' block_z ' num2str(tz)]);
                tx=block_x;
                ty=block_y;
                return
            end
        end

        % Increment splits based on largest dimension
        if npix_x / tx >= npix_y / ty
            tx = tx + 1;
        else
            ty = ty + 1;
        end
        
        % Check maximum splits
        if tx >= 5000 || ty >= 5000           
            error('Calculated at least 5000 blocks necessary. Aborting.');
        end
       
        bytes = ceil(npix_x / tx) * ceil(npix_y / ty) * npix_z * 4;
        if bytes <= maxblocksize
            
            return
        end
    end
end

function  process(inpath, tx, ty, tz, dxy, dz, numit, NA, rf, lambda_ex, lambda_em, fcyl, slitwidth, damping, clipval, stop_criterion, deviceUsed)
    startTime = datetime('now');
    tic;      
    %generate PSF
    Rxy = 0.61 * lambda_em / NA;
    dxy_corr = min(dxy, Rxy / 3);
    [psf, nxy, nz, FWHMxy, FWHMz] = LsMakePSF(dxy_corr, dz, NA, rf, lambda_ex, lambda_em);

    %start deconvolution
    scal  = deconvolve(inpath, psf, numit, damping, tx, ty, tz, clipval, stop_criterion, deviceUsed);
    outpath = fullfile(inpath, 'deconvolved');

    %write parameter info file
    disp('generating info file...');
    fid = fopen(fullfile(outpath, 'DECONV_parameters.txt'), 'w');
    fprintf(fid, '%s\r\n',['deconvolution finished at ' char(datetime('now'))]);
    fprintf(fid, '%s\r\n',['data processed on: ' deviceUsed]);
    elapsedTime = duration(datetime('now') - startTime);
    elapsedTime.Format = 'hh:mm';
    fprintf(fid, '%s\r\n',['elapsed time: ' char(elapsedTime)]);
    fprintf(fid, '%s\r\n', '');
    fprintf(fid, '%s\r\n',['highest intensity value in deconvolved data: ' num2str(scal)]);
    fprintf(fid, '%s\r\n',['focal length of cylinder lens (mm): ' num2str(fcyl)]);
    fprintf(fid, '%s\r\n',['width of slit aperture (mm): ' num2str(slitwidth)]);
    fprintf(fid,'%s\r\n', ['histogram clipping value (%): ' num2str(clipval)]);
    fprintf(fid, '%s\r\n',['numerical aperture: ' num2str(NA)]);
    fprintf(fid, '%s\r\n',['excitation wavelength (nm): ' num2str(lambda_ex)]);
    fprintf(fid, '%s\r\n',['emission wavelength (nm): ' num2str(lambda_em)]);
    fprintf(fid, '%s\r\n',['refractive index: ' num2str(rf)]);
    fprintf(fid, '%s\r\n',['max. iterations: ' num2str(numit)]);
    fprintf(fid, '%s\r\n',['damping factor (%): ' num2str(damping)]);
    fprintf(fid, '%s\r\n',['stop criterion (%): ' num2str(stop_criterion)]);
    fprintf(fid, '%s\r\n', '');
    fprintf(fid, '%s\r\n',['source data folder: ' inpath]);
    fprintf(fid, '%s\r\n',['number of blocks (x y z): ' num2str(tx)  ' x ' num2str(ty) ' x ' num2str(tz)]);
    fprintf(fid, '%s\r\n', '');
    fprintf(fid, '%s\r\n',['voxel size : ' num2str(dxy)  ' nm x ' num2str(dxy) ' nm x ' num2str(dz) ' nm']);
    fprintf(fid, '%s\r\n',['size of PSF (pixel): ' num2str(nxy)  ' x ' num2str(nxy) ' x ' num2str(nz)]);
    fprintf(fid, '%s\r\n',['FWHHM of PSF lateral (nm): ' num2str(FWHMxy)]);
    fprintf(fid, '%s\r\n',['FWHHM of PSF axial (nm): ' num2str(FWHMz)]);
    Rxy = 0.61 * lambda_em / NA;
    Rz = (2 * lambda_em * rf) / NA^2;
    fprintf(fid, '%s\r\n',['Rayleigh range of objective lateral (nm): ' num2str(Rxy)]);
    fprintf(fid, '%s\r\n',['Rayleigh range of objective axial (nm): ' num2str(Rz)]);
    fclose(fid);
    
    disp(['deconvolution of: ' inpath ' finished successfully']);
    disp(['elapsed time: ' char(elapsedTime)]);
    disp('----------------------------------------------------------------------------------------');
end

function scal = deconvolve(inpath, psf, numit, damping, numxblocks, numyblocks, numzblocks, clipval, stop_criterion, deviceUsed)

    try
        [info.x, info.y, info.z] = getstackinfo(inpath);
        fprintf('x %d y %d z %d\n', info.x, info.y, info.z);
        [p1, p2] = split(info, numxblocks, numyblocks, numzblocks);
        
        blocklist = strings(size(p1, 1), 1);
        rawmax = 0; deconvmax = 0; deconvmin = Inf;
        x = 1; y = 2; z = 3;

        totalBlocks = numxblocks * numyblocks * numzblocks;
        fprintf('Processing %d blocks...\n', totalBlocks);
    
        for i = 1 : size(blocklist, 1)

            %begin processing next block     
            disp( ['loading block ' num2str(i) ' from ' num2str(size(blocklist, 1))]);
            
            %load next block of data
            startp = p1(i, :); endp = p2(i, :);

            %load next block into memory     
            x1 = startp(x); x2 = endp(x);
            y1 = startp(y); y2 = endp(y);
            z1 = startp(z); z2 = endp(z);
            bl = load_block(inpath, x1, x2, y1, y2, z1, z2, deviceUsed);
                    
            %get min-max of raw data stack
            temp = max(bl(:));
            if temp > rawmax
                rawmax = temp;
            end
            
            %deconvolve current block of data
            disp(['Processing block ' num2str(i) ' from ' num2str(size(blocklist, 1))]);        
            bl = process_block(bl, psf, numit, damping, stop_criterion, deviceUsed);        
                    
            %find maximum value in processed data
            temp = max(bl(:));
            if temp > deconvmax
                deconvmax = temp;
            end
            temp = min(bl(:));
            if temp < deconvmin
                deconvmin = temp;
            end
            
            if totalBlocks > 1 % there is more than one block
                disp(['Saving block ' num2str(i) ' from ' num2str(size(blocklist, 1))]);
                blocklist(i) = [tempdir 'bl' num2str(i) '_temp.mat'];
                save(blocklist(i), 'bl', '-v7.3', '-nocompression');                

                if strcmpi(deviceUsed, 'GPU')
                    gpuDevice(1);
                end
            end
             % Progress update
             fprintf('Block %d/%d complete (%.1f%%)\n', i, totalBlocks, 100*i/totalBlocks);
        end
        
        if clipval > 0
            %estimate the global histogram and upper and lower clipping values
            nbins = 1e6;
            binwidth = deconvmax / nbins;
            bins = 0 : binwidth : deconvmax;
            
            %calculate cumulative histogram by scanning all blocks
            disp('calculating histogram...');
            
            chist = cumsum(histcounts(bl, bins)); %histogram from last block
            if totalBlocks > 1
                for i = size(blocklist, 1) - 1 : -1 : 1
                    S = load(blocklist(i), 'bl');
                    chist = chist + cumsum(histcounts(S.bl, bins));
                end
                clear S;
            end
            
            %determine upper and lower histogram clipping values
            chist = chist / max(chist) * 100; %normalize cumulative histogram to 0..100%
            
            low_clip = findClosest(chist, clipval) * binwidth;
            high_clip = findClosest(chist, 100-clipval) * binwidth;
        end
        
        %make folder for results and write tiff-images
        outpath = fullfile(inpath, 'deconvolved');
        
        if isfolder(outpath) == 0
            disp(['Creating results folder ' outpath]);
            mkdir(outpath);
        else
            delete([outpath '*.tif']);
        end
        
        %mount data and save data layeer by layer
        blnr = 1; imagenr = 1;
        for i = 1 : numzblocks
            disp(['mounting layer ' num2str(i) ' from ' num2str(numzblocks)]);
            
            %load and mount next layer of images
            if numxblocks * numyblocks * numzblocks > 1
                R = zeros(info.x, info.y, p2(blnr, z)- p1(blnr, z)+1, 'single');
                for j = 1 : numxblocks * numyblocks
                    S = load(blocklist(blnr), 'bl');
                   
                    delete(convertStringsToChars(blocklist(blnr)));
                    
                    R( p1(blnr, x) : p2(blnr, x), p1(blnr, y) : p2(blnr, y), :) = S.bl;
                    blnr = blnr + 1;
                end
                clear S;
            else %if there is only one block of data
                R = bl;
                clear bl;
            end
            
            %rescale deconvolved data
            if rawmax <= 65535
                scal = 65535;
            else
                scal = rawmax; %scale to maximum of input data
            end
            
            if clipval > 0
                %perform histogram clipping
                R(R < low_clip) = low_clip;
                R(R > high_clip) = high_clip;
                R = (R - low_clip) ./ (high_clip - low_clip) .* scal;
            else
                %otherwise scale to min..max
                R = (R - deconvmin) ./ (deconvmax - deconvmin) * scal;
            end
            
            %write images to folder
            disp('saving images...');
            
            for k = 1 : size(R, 3)
                s = num2str(imagenr);
                while length(s) < 3
                    s = strcat('0', s);
                end
                fullname = fullfile(outpath, ['deconv_' s '.tif']);
                
                if rawmax <= 65535 %16bit data
                    im = uint16(squeeze(R(:, :, k)));
                    imwrite(im', fullname);
                else %32 bit data
                    im = squeeze(R(:, :, k));
                    writeTiff32(im', fullname) %R must be single;
                end
                imagenr = imagenr + 1;
            end
        end

    catch ME
        cleanup_blocks(blocklist);
        if strcmpi(deviceUsed, 'GPU')
            gpuDevice(1);
        end
        rethrow(ME);
    end
end

function cleanup_blocks(blocklist)
    for i = 1:numel(blocklist)
        if ~isempty(blocklist(i)) && isfile(blocklist(i))
            delete(blocklist(i));
        end
    end
end

function bl = load_block(inpath, start_x, end_x, start_y, end_y, start_z, end_z, deviceUsed)	
    try

        filelist = dir(fullfile(inpath, '*.tif'));
        
        nx = end_x - start_x;
        ny = end_y - start_y;
        nz = end_z - start_z;
        
        bl = zeros(nx+1, ny+1, nz+1, 'single');
        for k = 1 : nz
            imgPath = fullfile(inpath, filelist((k-1)+start_z).name); 
            im = im2single(imread(imgPath, 'PixelRegion', {[start_y,end_y], [start_x,end_x]}));
            bl(:,:,k) = im'; % No transpose
        end	
    catch ME
        if exist('bl', 'var')
            clear bl;
        end
        if strcmpi(deviceUsed, 'GPU')
            gpuDevice(1); % Reset GPU on error
        end
        rethrow(ME);
    end
end


function block = process_block(block, psf, numit, damping, stopcrit, deviceUsed)
    try
        %for efficiency of FFT pad data in a way that the largest prime factor becomes <= 5
        blx = size(block, 1); 
        bly = size(block, 2); 
        blz = size(block, 3);

        pad_x = 0.5 * (findGoodFFTLength(blx + 4 * size(psf, 1)) - blx);
        pad_y = 0.5 * (findGoodFFTLength(bly + 4 * size(psf, 2)) - bly);
        pad_z = 0.5 * (findGoodFFTLength(blz + 4 * size(psf, 3)) - blz);
        
        block = padarray(block, [floor(pad_x) floor(pad_y) floor(pad_z)], 'pre', 'symmetric');
        block = padarray(block, [ceil(pad_x) ceil(pad_y) ceil(pad_z)], 'post', 'symmetric');
        
        %deconvolve block using Lucy-Richardson algorithm
        if strcmpi(deviceUsed, 'GPU')
            gpuDevice(1);
            block = deconGPU(block, psf, numit, damping, stopcrit);
            gpuDevice(1);
        else
            block = deconCPU(block, psf, numit, damping, stopcrit);
        end
    
        %remove padding
        block = block(floor(pad_x) : end-ceil(pad_x)-1, floor(pad_y) : end-ceil(pad_y)-1, floor(pad_z) : end-ceil(pad_z)-1);
    catch ME
        if strcmpi(deviceUsed, 'GPU')
            gpuDevice(1); % Reset GPU on error
        end
        rethrow(ME);
    end
end

%provides coordinates of subblocks after splitting
function [p1, p2] = split(info, nx, ny, nz)
    xw = ceil(info.x / nx);
    yw = ceil(info.y / ny);
    zw = ceil(info.z / nz);
    
    p1 = zeros(nx*ny*nz, 3);
    p2 = zeros(nx*ny*nz, 3);
    
    n = 0;
    for i = 0 : nz-1
        zs = i * zw + 1;
        for j = 0 : ny-1
            ys = j * yw + 1;
            for k = 0 : nx-1
                xs = k * xw + 1;
                n = n + 1;
                p1(n, 1) = xs;
                p2(n, 1) = min([xs + xw - 1, info.x]);
                
                p1(n, 2) = ys;
                p2(n, 2) = min([ys + yw - 1, info.y]);
                
                p1(n, 3) = zs;
                p2(n, 3) = min([zs + zw - 1, info.z]);
            end
        end
    end
end

%writes 32bit float tiff-imges
function writeTiff32(img, fname)
    t = Tiff(fname, 'w');
    tag.ImageHeight = size(img, 2);
    tag.ImageWidth = size(img, 1);
    tag.Compression = Tiff.Compression.LZW;
    tag.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tag.Photometric = Tiff.Photometric.MinIsBlack;
    tag.BitsPerSample = 32;
    tag.SamplesPerPixel = 1;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tag);
    t.write(img);
    t.close();
end

function index = findClosest(data, x)
    [~,index] = min(abs(data-x));
end

function x = findGoodFFTLength(x)
    while max(factor(x)) > 5
        x = x + 1;
    end
end

%Lucy-Richardson deconvolution
function deconvolved = deconCPU(stack, psf, niter, lambda, stop_criterion) 
    try

        deconvolved = stack;        
        OTF = single(psf2otf(psf, size(stack)));      
        R = 1/26 * ones(3, 3, 3, 'single'); 
        R(2,2,2) = single(0);

        deconvolved_new = zeros(size(stack), 'single');
        deltaL = Inf;
    
        for i = 1 : niter

            denom = convFFT(deconvolved, OTF);
            denom(denom < eps('single')) = eps('single'); % protect against division by zero
            
            if lambda == 0                       
                deconvolved_new = convFFT(stack ./ denom, conj(OTF)) .* deconvolved;
            else           
                deconvolved_new = (1 - lambda) .* convFFT(stack ./ denom, conj(OTF)) .* deconvolved ...
                    + lambda .* convn(deconvolved, R, 'same');
            end
            
            %estimate quality criterion
            delta = sqrt(sum((deconvolved(:) - deconvolved_new(:)).^2));
            if i == 1
                delta_rel = 0;
            else
                delta_rel = (deltaL - delta) / deltaL * 100;
            end
            
            deconvolved = deconvolved_new;
            deltaL = delta;
            
            disp(['iteration: ' num2str(i), ' delta: ' num2str(delta_rel, 3)]);   
            
            if i > 1 && delta_rel <= stop_criterion
                disp('stop criterion reached. Finishing iterations.');           
                break
            end              
        end
        
        %get rid of imaginary artifacts
        deconvolved = abs(deconvolved);

    catch ME
        % Memory cleanup
        clear deconvolved_new
        rethrow(ME)
    end
end


function deconvolved = deconGPU(stack, psf, niter, lambda, stop_criterion) 
    try

        deconvolved = gpuArray(stack);        
        psf_inv = gpuArray(psf(end:-1:1,end:-1:1, end:-1:1)); % spatially reversed psf 
        R = gpuArray(1/26 * ones(3, 3, 3, 'single')); 
        R(2,2,2) = 0;

        % Pre-allocate
        deconvolved_new = gpuArray.zeros(size(stack), 'single');
        deltaL = Inf;
    
        for i = 1 : niter 

            denom = convGPU(deconvolved, psf);
            denom(denom < eps('single')) = eps('single'); % protect against division by zero
            
            if lambda == 0                       
                deconvolved_new = convGPU(stack ./ denom, psf_inv) .* deconvolved;
            else           
                deconvolved_new = (1 - lambda) .* convGPU(stack ./ denom, psf_inv) .* deconvolved ...
                    + lambda .* convn(deconvolved, R, 'same');
            end
            
            %estimate quality criterion
            delta = sqrt(sum((deconvolved(:) - deconvolved_new(:)).^2));
            if i == 1
                delta_rel = 0;
            else
                delta_rel = (deltaL - delta) / deltaL * 100;
            end
            
            deconvolved = deconvolved_new;
            deltaL = delta;
            
            disp(['iteration: ' num2str(i), ' delta: ' num2str(delta_rel, 3)]);   
            
            if i > 1 && delta_rel <= stop_criterion
                disp('stop criterion reached. Finishing iterations.');           
                break
            end              
        end
        
        %get rid of imaginary artifacts
        deconvolved = gather(abs(deconvolved));
    catch ME
        gpuDevice(1);
        rethrow(ME);
    end
end

%deconvolve with OTF
function R = convFFT(data , otf)      
    R = ifftn(otf .* fftn(data));   
end

function R = convGPU(data, psf)
    R = gather(convn(gpuArray(data), psf, 'same'));
end

%calculates a theoretical point spread function
function [psf, nxy, nz, FWHMxy, FWHMz] = LsMakePSF(dxy, dz, NA, nf, lambda_ex, lambda_em)        
    disp('calculating PSF...');
    [nxy, nz, FWHMxy, FWHMz] = DeterminePSFsize(dxy, dz, NA, nf, lambda_ex, lambda_em);
      
    psf = samplePSF(dxy, dz, nxy, nz, NA, nf, lambda_ex, lambda_em);   
    disp('ok');
end


%determine the required grid size (xyz) for psf sampling
function [nxy, nz, FWHMxy, FWHMz] = DeterminePSFsize(dxy, dz, NA, nf, lambda_ex, lambda_em)      
    %Size of PSF grid is gridsize (xy z) times FWHM     
    gridsizeXY = 2;    
    gridsizeZ = 2;    
      
    halfmax = 0.5 .* LsPSFeq(0, 0, 0, NA, nf, lambda_ex, lambda_em);
    
    %find zero crossings
    fxy = @(x)LsPSFeq(x, 0, 0, NA, nf, lambda_ex, lambda_em) - halfmax;
    fz = @(x)LsPSFeq(0, 0, x, NA, nf, lambda_ex, lambda_em) - halfmax;    
    FWHMxy = 2 * abs(fzero(fxy, 100));
    FWHMz = 2 * abs(fzero(fz, 100));
      
    Rxy = 0.61 * lambda_em / NA;
    dxy_corr = min(dxy, Rxy / 3);
    
    nxy = ceil(gridsizeXY * FWHMxy / dxy_corr);
    nz = ceil(gridsizeZ * FWHMz / dz);
       
    %ensure that the grid dimensions are odd
    if mod(nxy, 2) == 0
        nxy = nxy + 1;
    end
    if mod(nz, 2) == 0
        nz = nz + 1;
    end    
end

function psf = samplePSF(dxy, dz, nxy, nz, NA, rf, lambda_ex, lambda_em)

    if mod(nxy, 2) == 0 || mod(nz, 2) == 0
        error('function samplePSF: nxy and nz must be odd!');
    end      
    
    psf = zeros((nxy - 1) / 2 + 1, (nxy - 1) / 2 + 1, (nz - 1) / 2 + 1, 'single');                        
    for z = 0 : (nz - 1) / 2
        for y = 0 : (nxy - 1) / 2
            for x = 0 : (nxy - 1) / 2
               psf(x+1, y+1, z+1) = LsPSFeq(x*dxy, y*dxy, z*dz, NA, rf, lambda_ex, lambda_em);   
            end        
        end
    end 
        
    %Since the PSF is symmetrical around all axes only the first Octand
    %is calculated for computation efficiency. The other 7 Octands are obtained by mirroring around
    %the respecitve axes
    psf = mirror8(psf);   
    
    %normalize psf to integral one
    psf = psf ./ sum(psf(:));
end

function R = mirror8(p1)    
    %mirrors the content of the first quadrant to all other quadrants to
    %obtain the complete PSF.
    
    sx = 2 * size(p1, 1) - 1; sy = 2 * size(p1, 2) - 1; sz = 2 * size(p1, 3) - 1;        
    cx = ceil(sx / 2); cy = ceil(sy / 2); cz = ceil(sz / 2);
    
    R = zeros(sx, sy, sz, 'single');
    R(cx:sx, cy:sy, cz:sz) = p1;
    R(cx:sx, 1:cy, cz:sz) = flip3D(p1, 0, 1 ,0);
    R(1:cx, 1:cy, cz:sz) = flip3D(p1, 1, 1, 0);
    R(1:cx, cy:sy, cz:sz) = flip3D(p1, 1, 0, 0);     
    R(cx:sx, cy:sy, 1:cz) = flip3D(p1, 0, 0, 1);
    R(cx:sx, 1:cy, 1:cz) =  flip3D(p1, 0, 1 ,1);
    R(1:cx, 1:cy, 1:cz) =  flip3D(p1, 1, 1, 1);
    R(1:cx, cx:sy, 1:cz) =  flip3D(p1, 1, 0, 1);
end

%utility function for mirror8
function R = flip3D(data, x, y, z)        
    R = data;
    if x
        R = flip(R, 1);
    end
    if y
        R = flip(R, 2);
    end
    if z
        R = flip(R, 3);
    end
end

%calculates PSF at point (x,y,z)
function R = LsPSFeq(x, y, z, NA, n, lambda_ex, lambda_em)                        
    R = PSF(x, y, z, NA, n, lambda_ex) .* PSF(x, y, z, NA, n, lambda_em);    
end

%utility function for LsPSFeq
function R = PSF(x, y, z, NA, n, lambda)                                
    f2 = @(p)f1(p, x, y, z, lambda, NA, n);    
    R = 4 .* abs(integral(f2, 0, 1, 'AbsTol',1e-3)).^2;
end

%utility function for LsPSFeq
function R = f1(p, x, y, z, lambda, NA, n)            
    R = besselj(0, 2 .* pi .* NA .* sqrt(x.^2 + y.^2) .* p ./ (lambda .* n))...
        .* exp(1i .* (-pi .* p.^2 .* z .* NA.^2) ./ (lambda .* n.^2)) .* p;    
end

function [x, y, z] = getstackinfo(datadir)
    filelist = dir(fullfile(datadir, '*.tif'));
    if numel(filelist) == 0
        x = 0; y = 0; z = 0;
    else
        test = imread(fullfile(datadir, filelist(1).name));
        x = size(test, 2);
        y = size(test, 1);
        z = numel(filelist);
    end
end