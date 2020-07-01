close all
figure()
for i=1:8
    [final_mesh, Rg] = adaptive_mesh_refinement1D(10, 0, 1, 10^(-i), 'fd1d_bvp_test03');
    Rg = [0; Rg; 0];
    semilogy(final_mesh, Rg, '-*', 'linewidth', 4)
    hold on
    set(gca, 'fontsize', 30)
    xlabel('x coordinate')
    ylabel('local norm of the gradient of the error')
    title('Advection-diffusion-reaction problem with space-dependent coefficients')
end
legend('tolerance = 1e-1', 'tolerance = 1e-2', 'tolerance = 1e-3', 'tolerance = 1e-4', 'tolerance = 1e-5' ,...
    'tolerance = 1e-6', 'tolerance = 1e-7', 'tolerance = 1e-8')