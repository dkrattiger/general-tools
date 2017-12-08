function [e_PHI,e_f,PHI,f] = eig_error(PHI_true,f_true,PHI,f)

[n_eig] = length(f_true);

% %% sort modes by frequency
% % ======================================================================= %
% [f_true,i_sort_true] = sort(f_true);
% PHI_true = PHI_true(:,i_sort);
% 
% [f,i_sort_] = sort(f);
% PHI = PHI(:,i_sort);

%% normalize modes to unit length
% ======================================================================= %
mag_true = sqrt(diag(PHI_true'*PHI_true));
PHI_true = PHI_true*diag(1./mag_true);

mag = sqrt(diag(PHI'*PHI));
PHI = PHI*diag(1./mag);

%% Sort test eigen-pairs to optimally match true eigen-pairs
% ======================================================================= %
% use linear programming algorithm to find optimal matching
Cmat = -abs(PHI'*PHI_true); % cost matrix
[i_sort] = match_points(Cmat); % dynamic programming matching
% i_sort = 1:length(f);
PHI = PHI(:,i_sort);
f = f(i_sort);

% abs(diag(PHI_true'*PHI))
% figure(121212);clf;
testvec1 = abs(diag(PHI_true'*PHI));

%% resolve degenerate modes to best matching modes
% ======================================================================= %

% frequency matching tolerance
tol_f = max(f_true)*5e-2;

% initialize counter
i_remain = 1:n_eig;

deg_label = zeros(n_eig,1);
count = 1;

while ~isempty(i_remain)
    
    i_match1 = find(abs(f_true(i_remain(1))-f_true(i_remain))<tol_f);
    i_match2 = find(abs(f_true(i_remain(i_match1(end)))-f_true(i_remain))<tol_f);
    i_match = union(i_match1,i_match2);
    
    
    i_deg = i_remain(i_match);i_remain(i_match) = [];
    n_deg = length(i_deg);
    
    % use Moore-Penrose inverse to combine degenerate modes for
    % least squares matching with true modes    
    C = ((PHI_true(:,i_deg)'*PHI_true(:,i_deg))\PHI_true(:,i_deg)')*PHI(:,i_deg);
    PHI_true(:,i_deg) = PHI_true(:,i_deg)*C;

    % normalize modes to unit length
    mag = sqrt(diag(PHI_true(:,i_deg)'*PHI_true(:,i_deg)));
    PHI_true(:,i_deg) = PHI_true(:,i_deg)*diag(1./mag);

    deg_label(i_deg) = count;
    count = count + 1;
end




%% Sort test eigen-pairs to optimally match true eigen-pairs
% ======================================================================= %

% preallocate arrays
e_PHI = zeros(size(f));

for i = 1:n_eig

    % full mode shape
    PHI_1 = PHI_true(:,i);
    PHI_1 = PHI_1/sqrt(PHI_1'*PHI_1);

    % BMS mode shape
    PHI_2 = PHI(:,i);
    PHI_2 = PHI_2/sqrt(PHI_2'*PHI_2);

    % error in test modes
    e_PHI(i) = 1-abs(PHI_1'*PHI_2);
    
end

% error in test frequencies
e_f = abs(f-f_true)./f_true;

% % return error vectors in original sort
% e_f(i_sort) = e_f; 
% e_PHI(i_sort) = e_PHI; 
