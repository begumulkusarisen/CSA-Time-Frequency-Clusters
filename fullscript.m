%1. SETUP (Revised)
clear; clc;

% --- CONFIGURATION ---
% Use an absolute base path to avoid "Folder Not Found" errors
base_dir = '/Users/begumulku/Desktop/coh_files'; 
output_folder = fullfile(base_dir, 'Processed_Network_Matrices'); 

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('Created output folder at: %s\n', output_folder);
end

target_networks = {'Default', 'Cont', 'DorsAttn'}; 


% --- FILE SEARCH ---
target_directory = '/Users/begumulku/Desktop/coh_files'; 
file_pattern = fullfile(target_directory, '*connectn_cohere_*'); 
files = dir(file_pattern); 

if isempty(files)
    error('No connectivity files found in the specified directory: %s', target_directory); 
end
fprintf('Found %d connectivity files in %s.\n', length(files), target_directory);

% Initialize Group Data containers
Group_Data = struct(); 
for n = 1:length(target_networks)
    Group_Data.(target_networks{n}) = [];
end
Subject_Map = struct();

%% 2. MAIN LOOP (Load -> Extract -> Save -> Stack)
for s = 1:length(files)
    
    full_path = fullfile(files(s).folder, files(s).name);
    
    % --- A. ID Extraction ---
    tokens = regexp(full_path, 'sub-(\d+)', 'tokens');
    if ~isempty(tokens)
        subj_id_str = ['sub' tokens{1}{1}];
    else
        [~, name_part] = fileparts(files(s).name);
        subj_id_str = name_part; 
    end
    
    fprintf('Processing %d/%d: %s ... \n', s, length(files), subj_id_str);
    
    % --- B. LOAD DATA ---
    data = load(full_path);
    if ~isfield(data, 'TF')
        warning('   Variable "TF" missing in this file. Skipping.');
        continue;
    end
    TF_Matrix = data.TF; 
    
    % 2. Check for Atlas/Labels
    if isfield(data, 'Atlas') && isfield(data.Atlas, 'Scouts')
        roi_labels = {data.Atlas.Scouts.Label};
    elseif isfield(data, 'RowNames')
        roi_labels = data.RowNames'; 
        if size(roi_labels, 1) > 1, roi_labels = roi_labels'; end
    else
        warning('   No ROI labels found (Atlas or RowNames). Skipping.');
        continue;
    end
    
    % --- C. NETWORK MAPPING ---
    network_names = cell(length(roi_labels), 1);
    for i = 1:length(roi_labels)
        lbl = roi_labels{i};
        idx = find(lbl == '_', 1, 'first');
        if isempty(idx)
            network_names{i} = lbl; 
        else
            network_names{i} = lbl(1:idx-1); 
        end
    end
    [unique_networks, ~, network_id_vector] = unique(network_names);
    
    % --- D. RECONSTRUCT & EXTRACT 3D MATRICES ---
    Ntime = size(TF_Matrix, 2);
    num_rois = length(roi_labels);
    mask_full_triu = triu(true(num_rois)); 
    
    Network_Within_Dyna_Connectivity = struct();
    
    for t = 1:Ntime
        R_roi = zeros(num_rois);
        vec = TF_Matrix(:, t);
        
        if length(vec) ~= sum(mask_full_triu(:))
            warning('   Vector size mismatch at time %d. Skipping.', t);
            continue;
        end
        
        R_roi(mask_full_triu) = vec;
        R_roi = R_roi + triu(R_roi, 1)'; 
        
        for net_idx = 1:length(unique_networks)
            net_name = unique_networks{net_idx};
            
            if ~ismember(net_name, target_networks), continue; end
            
            idx_rois = (network_id_vector == net_idx);
            sub_matrix = R_roi(idx_rois, idx_rois);
            
            if t == 1
                N_net = size(sub_matrix, 1);
                Network_Within_Dyna_Connectivity.(net_name) = NaN(N_net, N_net, Ntime);
            end
            
            Network_Within_Dyna_Connectivity.(net_name)(:, :, t) = sub_matrix;
        end
    end
    
    % --- E. SAVE INDIVIDUAL SUBJECT FILE ---
    save_filename = fullfile(output_folder, ['Network_Within_Dyna_Connectivity_' subj_id_str '.mat']);
    save(save_filename, 'Network_Within_Dyna_Connectivity');
    
    % --- F. FLATTEN FOR K-MEANS ---
    Subject_Map(s).Name = subj_id_str;
    Subject_Map(s).Points = Ntime;
    
    for n = 1:length(target_networks)
        net_name = target_networks{n};
        if ~isfield(Network_Within_Dyna_Connectivity, net_name), continue; end
        
        Matrix_3D = Network_Within_Dyna_Connectivity.(net_name);
        [n_sub_rois, ~, n_time] = size(Matrix_3D);
        n_edges = n_sub_rois * (n_sub_rois - 1) / 2;
        
        Flat_Features = zeros(n_time, n_edges);
        for t = 1:n_time
            current_slice = Matrix_3D(:, :, t);
            mask_tri = triu(true(n_sub_rois), 1); 
            Flat_Features(t, :) = abs(current_slice(mask_tri)); 
        end
        Group_Data.(net_name) = [Group_Data.(net_name); Flat_Features];
    end
    
    clear data TF_Matrix Network_Within_Dyna_Connectivity; 
end
fprintf('\nProcessing Complete. Starting Clustering...\n');


%% 3. GROUP CLUSTERING (K-MEANS)
Results = struct();
k = 5; 
for n = 1:length(target_networks)
    net_name = target_networks{n};
    X = Group_Data.(net_name);
    
    if isempty(X)
        fprintf('No data collected for network: %s\n', net_name);
        continue;
    end
    
    fprintf('Running K-Means (k=%d) for %s (Input: %d samples)...\n', k, net_name, size(X,1));
    % OPTIMIZATION: Reduced Replicates to 2 for faster clustering
    [idx_group, Centroids] = kmeans(X, k, 'Distance', 'sqeuclidean', 'Replicates', 2);
    
    Results.(net_name).All_Labels = idx_group;
    Results.(net_name).Centroids = Centroids;
    
    % Split back to Subjects
    curr_idx = 1;
    for s = 1:length(Subject_Map)
        if isempty(Subject_Map(s).Name), continue; end
        
        n_points = Subject_Map(s).Points;
        subj_labels = idx_group(curr_idx : curr_idx + n_points - 1);
        
        safe_name = regexprep(Subject_Map(s).Name, '[^a-zA-Z0-9]', '');
        Results.(net_name).Subjects.(['s_' safe_name]) = subj_labels;
        
        curr_idx = curr_idx + n_points;
    end
end
fprintf('\nDONE! Check "Results" variable and "Processed_Network_Matrices" folder.\n');



%% 4 & 5. TRANSITION PROBABILITY (TPM) & DWELL TIME
transition_rates = struct();
fprintf('\nStarting TPM and Dwell Time calculation for K=%d states...\n', k);

for n = 1:length(target_networks)
    net_name = target_networks{n};
    if ~isfield(Results, net_name), continue; end
    
    NetResults = Results.(net_name);
    subject_fields = fieldnames(NetResults.Subjects);
    num_subjects = length(subject_fields);
    
    % Initialize storage
    transition_rates.(net_name).TransitionMatrix = zeros(k, k, num_subjects); 
    transition_rates.(net_name).DwellTime = zeros(num_subjects, k); % Subj x States
    transition_rates.(net_name).SubjectIDs = cell(num_subjects, 1);
    
    for s = 1:num_subjects
        subj_field = subject_fields{s};
        subj_id = strrep(subj_field, 's_', ''); 
        transition_rates.(net_name).SubjectIDs{s} = subj_id;
        
        state_labels = NetResults.Subjects.(subj_field);
        
        % --- 1. Transition Probability Matrix (TPM) ---
        from_states = state_labels(1:end-1);
        to_states = state_labels(2:end);
        linear_indices = sub2ind([k, k], from_states, to_states);
        TPM_counts = accumarray(linear_indices, 1, [k*k, 1], @sum, 0);
        TPM = reshape(TPM_counts, k, k); 
        row_sums = sum(TPM, 2); 
        TPM = TPM ./ row_sums;
        TPM(isnan(TPM)) = 0;
        transition_rates.(net_name).TransitionMatrix(:, :, s) = TPM;
        
        % --- 2. Dwell Time (Empirical Calculation) ---
        % We find the length of consecutive identical state occurrences
        for state_i = 1:k
            % Create a binary vector for current state
            is_state = (state_labels == state_i);
            if ~any(is_state)
                transition_rates.(net_name).DwellTime(s, state_i) = 0;
                continue;
            end
            
            % Identify switches (where state starts/ends)
            % diff([0; is_state; 0]) finds 1 for start and -1 for end
            changes = diff([0; is_state; 0]);
            run_starts = find(changes == 1);
            run_ends = find(changes == -1);
            run_lengths = run_ends - run_starts;
            
            % Average duration in this state
            transition_rates.(net_name).DwellTime(s, state_i) = mean(run_lengths);
        end
    end
end
%% 6. COLORFUL TPM HEATMAPS
fprintf('\n--- Transition Rate Mapping (TPMs) ---\n');

network_names = fieldnames(transition_rates);

for i = 1:length(network_names)
    net_name = network_names{i};
    NetData = transition_rates.(net_name);

    if ~isfield(NetData, 'TransitionMatrix')
        continue;
    end

    TPM_3D = NetData.TransitionMatrix;
    SubjectIDs = NetData.SubjectIDs;

    fprintf('\nNetwork: %s\n', net_name);

end


%% 7. GROUP-AVERAGED TPM MAPS (across all subjects)
fprintf('\n--- Group-Averaged Transition Probability Maps ---\n');

network_names = fieldnames(transition_rates);

for i = 1:length(network_names)
    net_name = network_names{i};
    NetData = transition_rates.(net_name);

    if ~isfield(NetData, 'TransitionMatrix')
        warning('No TPM data found for network %s. Skipping.', net_name);
        continue;
    end

    TPM_3D = NetData.TransitionMatrix;   % k x k x subjects
    num_subjects = size(TPM_3D, 3);

    if num_subjects == 0
        warning('No subjects found for network %s.', net_name);
        continue;
    end

    % --- Compute group mean TPM ---
    GroupAvgTPM = mean(TPM_3D, 3);

    fprintf('Averaged TPM computed for network: %s (%d subjects)\n', ...
            net_name, num_subjects);

    % --- Plot heatmap ---
    figure;
    h = heatmap(GroupAvgTPM);
    h.Colormap = parula;
    h.ColorLimits = [0 1];
    h.Title = sprintf('Group-Averaged TPM - %s Network', net_name);
    h.XLabel = 'TO State';
    h.YLabel = 'FROM State';

    % State labels
    h.XDisplayLabels = string(1:k);
    h.YDisplayLabels = string(1:k);
end

%%
fprintf('\n--- Individual Subject Dwell Times (by Network) ---\n');

for n = 1:length(target_networks)
    net_name = target_networks{n};
    if ~isfield(transition_rates, net_name), continue; end
    
    % Retrieve Data
    dwell_data = transition_rates.(net_name).DwellTime;
    sub_ids = transition_rates.(net_name).SubjectIDs;
    
    fprintf('\nNETWORK: %s\n', upper(net_name));
    fprintf('%-25s | %-10s | %-10s | %-10s\n', 'Subject ID', 'State 1', 'State 2', 'State 3');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for s = 1:length(sub_ids)
        fprintf('%-25s | %-10.2f | %-10.2f | %-10.2f\n', ...
            sub_ids{s}, dwell_data(s, 1), dwell_data(s, 2), dwell_data(s, 3));
    end
    
    % Optional: Create a MATLAB Table for easy export
    T = table(sub_ids, dwell_data(:,1), dwell_data(:,2), dwell_data(:,3), ...
        'VariableNames', {'SubjectID', 'State1_Dwell', 'State2_Dwell', 'State3_Dwell'});
    
    % Save table to the structure for later use
    transition_rates.(net_name).DwellTable = T;
    
    % Uncomment the line below if you want to save to CSV automatically:
    % writetable(T, fullfile(output_folder, [net_name, '_DwellTimes.csv']));
end

%%
excel_filename = fullfile(output_folder, 'Subject_Dwell_Times_Results.xlsx');

% Check if file exists and delete it to ensure a fresh export
if exist(excel_filename, 'file')
    delete(excel_filename);
end

for n = 1:length(target_networks)
    net_name = target_networks{n};
    
    if ~isfield(transition_rates, net_name)
        continue;
    end
    
    % Retrieve Data
    dwell_data = transition_rates.(net_name).DwellTime;
    sub_ids = transition_rates.(net_name).SubjectIDs;
    
    % Create Column Names for the Table
    colNames = {'SubjectID'};
    for state_idx = 1:k
        colNames{end+1} = sprintf('State_%d_DwellTime', state_idx);
    end
    
    % Create Table
    T_export = table(sub_ids, 'VariableNames', {'SubjectID'});
    
    % Concatenate the dwell data columns
    for state_idx = 1:k
        T_export.(sprintf('State_%d_DwellTime', state_idx)) = dwell_data(:, state_idx);
    end
    
    % Write to Excel: Each network gets its own sheet
    writetable(T_export, excel_filename, 'Sheet', net_name);
    
    fprintf('   - Successfully exported %s network to sheet: %s\n', net_name, net_name);
end