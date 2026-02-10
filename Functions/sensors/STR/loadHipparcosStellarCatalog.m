function catalog = loadHipparcosStellarCatalog(magnitude_limit, force_reload)
%==========================================================================
% loadHipparcosStellarCatalog: Load and filter Hipparcos star catalog with
%                              automatic caching in Data/ folder.
%
% Inputs:
%   magnitude_limit - Absolute magnitude upper limit for filtering.
%   force_reload    - Boolean flag to bypass cached .mat file and force
%                     re-parsing of hip_main.dat (default: false).
%
% Outputs:
%   catalog         - Structure containing filtered star data:
%                     .HIP_ID      - Nx1 Hipparcos identifiers (integer).
%                     .RA_deg      - Nx1 Right Ascension, J2000 [deg].
%                     .DEC_deg     - Nx1 Declination, J2000 [deg].
%                     .magnitude   - Nx1 Visual magnitude [mag].
%                     .pmRA        - Nx1 Proper motion in RA*cos(δ) [mas/yr].
%                     .pmDEC       - Nx1 Proper motion in DEC [mas/yr].
%                     .r_ECI       - 3xN unit vectors in ECI frame (J2000).
%
% Method:
%   First call: parses hip_main.dat (VizieR I/239 fixed-width format),
%   filters by magnitude, converts to ECI unit vectors, and saves to
%   Data/hipparcos_mag_X.X.mat. Subsequent calls load the .mat file.
%
% Notes:
%   - Requires hip_main.dat (50.8 MB) in Data/ or root dir.
%   - Download: http://cdsarc.u-strasbg.fr/ftp/cats/I/239/hip_main.dat
%   - Cached files: Data/hipparcos_mag_X.X.mat (magnitude-specific).
%   - Creates Data/ folder automatically if missing.
%==========================================================================

    % Validate and set default inputs
    if nargin < 2
        force_reload = false;
    end

    % Ensure Data/ directory exists
    data_dir = [pwd, '\Data'];
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
        fprintf('Created directory: %s/\n', data_dir);
    end

    % Define cached file path in Data/ folder
    mat_file = fullfile(data_dir, sprintf('hipparcos_mag_%.1f.mat', magnitude_limit));

    % Attempt to load pre-filtered catalog from cache
    if exist(mat_file, 'file') && ~force_reload
        fprintf('=== Loading cached Hipparcos catalog ===\n');
        fprintf('File: %s\n', mat_file);
        
        load(mat_file, 'catalog', 'N_stars', 'magnitude_limit_saved');
        
        % Validate magnitude consistency
        if magnitude_limit_saved ~= magnitude_limit
            warning('Cached magnitude limit (%.1f) mismatch. Re-parsing...', magnitude_limit_saved);
        else
            fprintf('Loaded: %d stars with M <= %.1f\n', N_stars, magnitude_limit);
            return;
        end
    end

    % Cache miss: locate and parse raw hip_main.dat
    fprintf('=== Parsing hip_main.dat (VizieR I/239) ===\n');
    
    % Search for hip_main.dat in Data/ folder first, then root
    hip_file_root = 'hip_main.dat';
    hip_file_data = fullfile(data_dir, hip_file_root);
    
    if exist(hip_file_data, 'file')
        hip_file = hip_file_data;
    elseif exist(hip_file_root, 'file')
        hip_file = hip_file_root;
    else
        error(['File not found: hip_main.dat\n' ...
               'Expected locations:\n' ...
               '  - %s\n' ...
               '  - %s\n' ...
               'Download from: http://cdsarc.u-strasbg.fr/ftp/cats/I/239/hip_main.dat'], ...
               hip_file_data, hip_file_root);
    end
    
    fprintf('Reading from: %s\n', hip_file);
    
    fid = fopen(hip_file, 'r');
    if fid == -1
        error('Failed to open file: %s', hip_file);
    end
    
    fprintf('Parsing fixed-width ASCII...\n');
    t_start = tic;

    % Pre-allocate buffers for full catalog (118,218 entries)
    buffer_size = 120000;
    HIP_ID_buf  = zeros(buffer_size, 1);
    RA_deg_buf  = zeros(buffer_size, 1);
    DEC_deg_buf = zeros(buffer_size, 1);
    Amag_buf    =   nan(buffer_size, 1);
    pmRA_buf    = zeros(buffer_size, 1);
    pmDEC_buf   = zeros(buffer_size, 1);

    n_read = 0;
    while ~feof(fid)
        line = fgetl(fid);
        
        if ischar(line) && length(line) >= 150
            n_read = n_read + 1;
            
            % Parse fixed-width columns (VizieR I/239 byte alignment)
            HIP_str  = strtrim(line(9:14));   % HIP identifier
            Amag_str = strtrim(line(42:46));  % Absolute magnitude
            RA_str   = strtrim(line(52:63));  % RA (ICRS, J2000)
            DEC_str  = strtrim(line(65:76));  % DEC (ICRS, J2000)
            pmRA_str = strtrim(line(88:95));  % PM in RA*cos(δ)
            pmDE_str = strtrim(line(97:104)); % PM in DEC
            
            % Convert to numeric (skip if empty/invalid)
            if ~isempty(HIP_str),  HIP_ID_buf(n_read)  = str2double(HIP_str);  end
            if ~isempty(Amag_str), Amag_buf(n_read)    = str2double(Amag_str); end
            if ~isempty(RA_str),   RA_deg_buf(n_read)  = str2double(RA_str);   end
            if ~isempty(DEC_str),  DEC_deg_buf(n_read) = str2double(DEC_str);  end
            if ~isempty(pmRA_str), pmRA_buf(n_read)    = str2double(pmRA_str); end
            if ~isempty(pmDE_str), pmDEC_buf(n_read)   = str2double(pmDE_str); end
        end
        
        % Progress indicator every 20k entries
        if mod(n_read, 20000) == 0
            fprintf('  %d entries parsed...\n', n_read);
        end
    end
    
    fclose(fid);
    t_parse = toc(t_start);
    
    fprintf('Parsing complete: %d entries in %.2f sec\n', n_read, t_parse);

    % Trim buffers to actual read size
    HIP_ID_buf  =  HIP_ID_buf(1:n_read);
    RA_deg_buf  =  RA_deg_buf(1:n_read);
    DEC_deg_buf = DEC_deg_buf(1:n_read);
    Amag_buf    =    Amag_buf(1:n_read);
    pmRA_buf    =    pmRA_buf(1:n_read);
    pmDEC_buf   =   pmDEC_buf(1:n_read);

    % Filter by magnitude limit (exclude NaN and out-of-range)
    fprintf('Filtering Vmag <= %.1f...\n', magnitude_limit);
    
    valid_mask = ~isnan(Amag_buf) & (Amag_buf <= magnitude_limit);
    
    catalog.HIP_ID    = HIP_ID_buf(valid_mask);
    catalog.RA_deg    = RA_deg_buf(valid_mask);
    catalog.DEC_deg   = DEC_deg_buf(valid_mask);
    catalog.magnitude = Amag_buf(valid_mask);
    catalog.pmRA      = pmRA_buf(valid_mask);
    catalog.pmDEC     = pmDEC_buf(valid_mask);
    
    N_stars = length(catalog.HIP_ID);
    fprintf('Filtered catalog: %d stars\n', N_stars);

    % Convert spherical coordinates to ECI unit vectors (J2000)
    % r_ECI = [cos(δ)*cos(α), cos(δ)*sin(α), sin(δ)]^T
    catalog.r_ECI = [cosd(catalog.DEC_deg) .* cosd(catalog.RA_deg), ...
                     cosd(catalog.DEC_deg) .* sind(catalog.RA_deg), ...
                     sind(catalog.DEC_deg)]';

    % Save filtered catalog to .mat in Data/ folder
    fprintf('Saving to cache: %s\n', mat_file);
    magnitude_limit_saved = magnitude_limit;
    save(mat_file, 'catalog', 'N_stars', 'magnitude_limit_saved', '-v7.3');
    
    file_info = dir(mat_file);
    fprintf('Cache size: %.2f MB\n', file_info.bytes / 1e6);

    % Print catalog statistics
    fprintf('\n--- Catalog Statistics ---\n');
    fprintf('  Magnitude range:  [%6.2f, %6.2f]\n',     min(catalog.magnitude),  max(catalog.magnitude));
    fprintf('  Magnitude mean:   %6.2f ± %.2f\n',       mean(catalog.magnitude), std(catalog.magnitude));
    fprintf('  RA coverage:      [%7.2f, %7.2f] deg\n', min(catalog.RA_deg),     max(catalog.RA_deg));
    fprintf('  DEC coverage:     [%7.2f, %7.2f] deg\n', min(catalog.DEC_deg),    max(catalog.DEC_deg));
    
    % Proper motion statistics (exclude outliers |PM| > 1000 mas/yr)
    pmRA_valid  = catalog.pmRA( abs( catalog.pmRA) < 1000);
    pmDEC_valid = catalog.pmDEC(abs(catalog.pmDEC) < 1000);
    fprintf('  PM RA (typical):  %6.2f ± %.2f mas/yr\n',  mean(pmRA_valid), std(pmRA_valid));
    fprintf('  PM DEC (typical): %6.2f ± %.2f mas/yr\n', mean(pmDEC_valid), std(pmDEC_valid));

    % Magnitude distribution histogram
    mag_bins   = 0:1:7;
    mag_counts = histcounts(catalog.magnitude, mag_bins);
    fprintf('\n  Magnitude distribution:\n');
    for i = 1:length(mag_counts)
        fprintf('    [%d, %d) mag: %5d stars\n', mag_bins(i), mag_bins(i+1), mag_counts(i));
    end

    % Verify known bright stars (sanity check)
    fprintf('\n--- Verification: Bright Stars ---\n');
    verifyBrightStar(catalog, 60718, 'Acrux');        % Alpha Crucis
    verifyBrightStar(catalog, 21421, 'Aldebaran');    % Alpha Tauri
    verifyBrightStar(catalog, 677,   'Alpheratz');    % Alpha Andromedae
    verifyBrightStar(catalog, 97649, 'Altair');       % Alpha Aquilae
    verifyBrightStar(catalog, 62434, 'Becrux');       % Beta Crucis (Mimosa)
    verifyBrightStar(catalog, 25336, 'Bellatrix');    % Gamma Orionis
    verifyBrightStar(catalog, 27989, 'Betelgeuse');   % Alpha Orionis
    verifyBrightStar(catalog, 61084, 'Gacrux');       % Gamma Crucis
    verifyBrightStar(catalog, 59747, 'Imai');         % Delta Crucis
    verifyBrightStar(catalog, 11767, 'Polaris');      % Alpha Ursae Minoris
    verifyBrightStar(catalog, 70890, 'Proxima Cen');  % Proxima Centauri
    verifyBrightStar(catalog, 32349, 'Sirius');       % Alpha Canis Majoris
    verifyBrightStar(catalog, 91262, 'Vega');         % Alpha Lyrae
    verifyBrightStar(catalog, 63608, 'Vindemiatrix'); % Epsilon Virginis

    fprintf('\n=== Catalog ready ===\n');
end

%==========================================================================
% verify_bright_star: Check presence and properties of a known star.
%==========================================================================
function verifyBrightStar(catalog, hip_id, star_name)
    idx = find(catalog.HIP_ID == hip_id, 1);
    if ~isempty(idx)
        fprintf('  %-12s (HIP %6d): RA = %7.2f°, DEC = %7.2f°, M = %5.2f\n', ...
                star_name, hip_id, ...
                catalog.RA_deg(idx), catalog.DEC_deg(idx), catalog.magnitude(idx));
    else
        fprintf('  %-12s (HIP %6d): Not in catalog (M > limit)\n', star_name, hip_id);
    end
end