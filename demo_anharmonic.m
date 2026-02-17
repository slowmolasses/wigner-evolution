% Requires: wignerHO.m (for function definitions) to be on the MATLAB path,
% or copy the four functions (gridint, Hermite, SHOmode, W) into the same file.

dz = 0.1;
zgrid = -15:dz:15;
firstncoh = 8;

% Build Hamiltonian eigenstates
[shomodes, ~] = deal(zeros(firstncoh, length(zgrid)));
for i = 1:firstncoh
    shomodes(i,:) = SHOmode(zgrid, i-1);
end

% Anharmonic Hamiltonian (finite difference)
symize = @(M)(M + M');
hmat = diag((1/dz)^2 + 0.5*zgrid.^2) ...
     - (1/dz)^2/2 * symize([zeros(length(zgrid)-1,1) eye(length(zgrid)-1); zeros(1,length(zgrid))]);
hmat_anho = hmat + 0.02 * diag(zgrid.^4);

[anhovecs, anhovalmat] = eigs(hmat_anho, firstncoh, 'smallestabs');
firstnanhovals = diag(anhovalmat);
anhovecnorms = arrayfun(@(v) sqrt(gridint(zgrid, abs(v{:}.').^2)), num2cell(anhovecs, 1));
normalized_anhovecs = anhovecs ./ anhovecnorms;

% Initial displaced Gaussian
cohstate = SHOmode(zgrid - 2, 0);
coeffsanhocoh = gridint(zgrid, cohstate .* normalized_anhovecs');

% Time evolution
wtanhogrid = 0:0.5:250;
anhocohmovie = coeffsanhocoh .* exp(-1i * wtanhogrid.' * firstnanhovals.') * normalized_anhovecs.';

% Plot Wigner function at revival time ωt = 154  (frame index 309)
[Wrevival, modzgrid, modpgrid] = W(anhocohmovie(310,:), dz);
imagesc(modzgrid, modpgrid, Wrevival);
colorbar;  daspect([1 1 1]);
title({'W(z,p) for a displaced Gaussian state', ...
       'evolved by a quantum anharmonic oscillator', ...
       'after dissipating and reviving (\omega{}t=154)'});
xlabel('z (in units of z_0 m^{-1})');
ylabel('p (in modified units)');