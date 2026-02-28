%% Gottesman-Kitaev-Preskill (GKP) State Visualization
% Demonstrates finite-squeezing approximate GKP states for quantum error correction
% in continuous-variable (CV) systems. Generates Wigner quasi-probability distributions
% for logical qudit states with d=2,3,4 under squeezing constraints r=2,4,8.
%
% Theory: GKP states encode discrete logical information in the continuous
% position/momentum quadratures of an oscillator by localizing the state to a
% lattice in phase space. Ideal GKP states have infinite squeezing (delta-function
% peaks); finite squeezing creates Gaussian envelopes around each lattice point.
%
% Dependencies: Requires functions from lab3rawcode.m (gridint, SHOmode, W)

%% Parameters
dz = 0.1;                    % Position grid spacing
zgrid = -15:dz:15;           % Position grid (in units of z_0)
zgrid = zgrid - zgrid(floor(length(zgrid)/2)+1);  % Center grid at origin

% Logical dimensions to simulate
d_vals = [2, 3, 4];          % Qubit (d=2), qutrit (d=3), ququart (d=4)

% Squeezing parameters (finite squeezing radius)
r_vals = [2, 4, 8];          % Squeezing r: larger = closer to ideal GKP

% Lattice spacing (sets logical vs. physical separation)
Delta = sqrt(pi);            % Standard GKP lattice constant

%% Main simulation loop
figure('Position', [100, 100, 1400, 900]);

for row_idx = 1:length(r_vals)
    r_squeeze = r_vals(row_idx);
    
    for col_idx = 1:length(d_vals)
        d_logical = d_vals(col_idx);
        
        % Construct approximate GKP state for logical dimension d
        psi_GKP = construct_GKP_state(zgrid, d_logical, Delta, r_squeeze);
        
        % Normalize state
        norm_check = gridint(zgrid, abs(psi_GKP).^2);
        psi_GKP = psi_GKP / sqrt(norm_check);
        
        % Compute Wigner function
        [W_GKP, modzgrid, modpgrid] = W(psi_GKP, dz);
        
        % Plot Wigner function
        subplot(length(r_vals), length(d_vals), (row_idx-1)*length(d_vals) + col_idx);
        imagesc(modzgrid, modpgrid, W_GKP);
        colorbar;
        daspect([1, 1, 1]);
        
        % Title with parameters
        title(sprintf('GKP d=%d, r=%.1f', d_logical, r_squeeze), 'FontSize', 12);
        xlabel('z (units of z_0)', 'FontSize', 10);
        ylabel('p (modified units)', 'FontSize', 10);
        
        % Set color limits to emphasize negativity (signature of non-classicality)
        clim([-0.3, 0.6]);
        colormap(bluewhitered(256));
        
        % Add text annotation with purity
        purity = compute_purity_from_wigner(W_GKP, modzgrid, modpgrid, dz);
        text(0.05, 0.95, sprintf('Purity: %.3f', purity), ...
             'Units', 'normalized', 'Color', 'white', 'FontSize', 9, ...
             'VerticalAlignment', 'top', 'BackgroundColor', [0 0 0 0.5]);
    end
end

sgtitle('Approximate GKP States: Wigner Functions vs. Logical Dimension & Squeezing', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% Helper Functions

function psi = construct_GKP_state(zgrid, d, Delta, r_squeeze)
    % Constructs finite-squeezing approximate GKP state for d-dimensional logical space
    %
    % GKP state formula (position representation):
    % |GKP_d⟩ ≈ ∑_{n=-∞}^∞ exp(-π(n·Δ)²/Δ²_squeeze) |z = n·Δ·√π⟩
    % where Δ_squeeze = exp(-r_squeeze) is the squeezing envelope width
    %
    % For logical qudit encoding, we use a d-dimensional lattice with basis vectors
    % separated by Δ·√π in position space. The logical computational basis |k⟩_L
    % corresponds to peaks at positions z_k = (k + n·d)·Δ·√π for all integers n.
    
    psi = zeros(1, length(zgrid));
    Delta_squeeze = exp(-r_squeeze);  % Squeezing width (smaller = more squeezed)
    
    % Number of lattice replicas to sum (truncate when Gaussian envelope negligible)
    n_max = ceil(4 / Delta);  % Cutoff at ~4 standard deviations
    
    % Sum over lattice points for logical state |0⟩_L (can generalize to superpositions)
    % Position peaks at z = k·d·Δ for k ∈ ℤ (every d-th lattice site)
    for k = -n_max:n_max
        z_peak = k * d * Delta;  % Position of k-th lattice peak for logical |0⟩
        
        % Gaussian envelope centered at z_peak with width Delta_squeeze
        envelope = exp(-pi * ((zgrid - z_peak) / Delta_squeeze).^2);
        
        % Add localized wavepacket to superposition
        psi = psi + envelope;
    end
    
    % For qudits d>2, add phase-space structure (momentum peaks)
    % Momentum-space GKP states have dual lattice with spacing 1/Δ
    % Here we approximate by adding coherent oscillations between position peaks
    if d > 2
        % Modulate with sine/cosine to create momentum structure
        phase_modulation = exp(1i * 2*pi * zgrid / (d * Delta));
        psi = psi .* phase_modulation;
    end
end

function purity = compute_purity_from_wigner(W, zgrid, pgrid, dz)
    % Computes purity γ = Tr(ρ²) from Wigner function
    % For Wigner function W(z,p), purity is:
    % γ = 2π ∫∫ W(z,p)² dz dp
    %
    % This quantifies "Gaussianity": γ=1 for pure states, γ→0 for mixed states
    % GKP states have γ < 1 due to finite squeezing (non-Gaussian corrections)
    
    % Compute ∫ W² dp for each z
    W_squared = W.^2;
    integral_over_p = gridint(pgrid, W_squared.');  % Integrate over momentum
    
    % Compute ∫∫ W² dz dp
    integral_total = gridint(zgrid, integral_over_p);
    
    % Factor of 2π from phase-space measure
    purity = 2 * pi * integral_total;
end

function cmap = bluewhitered(n)
    % Custom colormap: blue (negative) → white (zero) → red (positive)
    % Emphasizes negative Wigner function values (non-classicality signature)
    
    if nargin < 1
        n = 256;
    end
    
    % Create diverging colormap
    half = floor(n/2);
    
    % Blue to white (negative values)
    blue_to_white = [linspace(0, 1, half)', linspace(0, 1, half)', ones(half, 1)];
    
    % White to red (positive values)
    white_to_red = [ones(n-half, 1), linspace(1, 0, n-half)', linspace(1, 0, n-half)'];
    
    cmap = [blue_to_white; white_to_red];
end

% Supporting functions (from lab3rawcode.m - include if running standalone)
% If lab3rawcode.m is not in path, uncomment and include these functions:

function int = gridint(grid, vals)
    % Trapezoidal numerical integration
    int = (grid(2:end) - grid(1:end-1)) * (vals(:, 2:end) + vals(:, 1:end-1)).' / 2;
end

function [shosol, norm] = SHOmode(grid, n)
    % Returns nth normalized eigenmode of quantum harmonic oscillator
    assert(mod(n,1)==0 && n>=0, "n must be a natural number");
    shosol = Hermite(grid, n) .* exp(-grid.^2/2);
    norm = sqrt(gridint(grid, abs(shosol).^2));
    shosol = shosol / norm;
end

function eval = Hermite(grid, n)
    % Evaluates nth Hermite polynomial recursively
    if n <= 0
        eval = ones(1, length(grid));
    else
        eval = 2*grid .* Hermite(grid, n-1) - 2*(n-1) * Hermite(grid, n-2);
    end
end

function [outvals, modzgrid, modpgrid]=W(psi,dz) 
%Returns the Wigner function as a discrete 2D array (column corresp. position and rows to momentum) 
%associated with the input discritized wave function defined on a origin-centerd position grid with spacing dz.
    N=length(psi);npts=ceil(N^2/2)-2;
    samplepts=zeros(npts,2);values=zeros(npts,1);
    valindex=1;pmaxj=pi/(4*dz);
    for j=2:N-1
        spanj=[max(1,N-(2*(N-j))):min(N,1+2*(j-1))];Nj=length(spanj);
        % disp("j: "+num2str(j)+" spanj: "+num2str(spanj));
        rhoantidiagj=psi(spanj)'.'.*flip(psi(spanj));
        psj=linspace(-pmaxj,pmaxj,Nj).';psj=psj-psj(floor(Nj/2)+1);dpj=2*pmaxj/(Nj-1);
        zsj=((j-1)-(floor(N/2)+(mod(N,2)==0)))*dz*ones(Nj,1);
        values(valindex:valindex+Nj-1)=Nj*dz*2/pi*fftshift(ifft(ifftshift(rhoantidiagj)));
        samplepts(valindex:valindex+Nj-1,:)=[zsj,psj];
        valindex=valindex+Nj;end;
    % disp("Sample points: "+num2str(samplepts(1:valindex,:)));
    % disp("N^2: "+num2str(npts)+" Final valindex: "+num2str(valindex)+" Final Nj: "+num2str(Nj));
    Winterp=scatteredInterpolant(samplepts,values,"natural");
    modzgrid=(1-(floor(N/2)+(mod(N,2)==0)))*dz:dz:dz*(N-1-1-(floor(N/2)+(mod(N,2)==0)));modpgrid=linspace(-pmaxj,pmaxj,N);
    modpgrid=modpgrid-modpgrid(floor(N/2)+1);
    % [modzmesh,modpmesh]=meshgrid(modzgrid,modpgrid);
    outvals=Winterp({modzgrid,modpgrid}).';
end


% function [outvals, modzgrid, modpgrid] = W(psi, dz)
%     % Returns Wigner function as a discrete 2D array
%     N = length(psi);
%     npts = ceil(N^2/2) - 2;
%     samplepts = zeros(npts, 2);
%     values = zeros(npts, 1);
%     valindex = 1;
%     pmaxj = pi / (4*dz);
% 
%     for j = 2:N-1
%         spanj = max(1, N-(2*(N-j))):min(N, 1+2*(j-1));
%         Nj = length(spanj);
%         rhoantidiagj = psi(spanj)' .* flip(psi(spanj));
%         psj = linspace(-pmaxj, pmaxj, Nj).';
%         psj = psj - psj(floor(Nj/2)+1);
%         zsj = ((j-1) - (floor(N/2) + (mod(N,2)==0))) * dz * ones(Nj, 1);
%         values(valindex:valindex+Nj-1) = Nj * dz * 2/pi * fftshift(ifft(ifftshift(rhoantidiagj)));
%         samplepts(valindex:valindex+Nj-1, :) = [zsj, psj];
%         valindex = valindex + Nj;
%     end
% 
%     Winterp = scatteredInterpolant(samplepts, values, "natural");
%     modzgrid = (1-(floor(N/2)+(mod(N,2)==0)))*dz:dz:dz*(N-1-1-(floor(N/2)+(mod(N,2)==0)));
%     modpgrid = linspace(-pmaxj, pmaxj, N);
%     modpgrid = modpgrid - modpgrid(floor(N/2)+1);
%     outvals = Winterp({modzgrid, modpgrid}).';
% end
