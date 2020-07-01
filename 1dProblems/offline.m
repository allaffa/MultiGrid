%% Offline creation of the dataset

num_experiments = 50000;

n0 = 201;
x0 = 0; 
x1 = 1;
epsilon = 1e-8;

grid = linspace(x0, x1, n0)'; 

diffusion_coeff_range = [0 10];
advection_coeff_range = [-10 10];
reaction_coeff_range = [-10 10];
rhs_range = [-10 10];

% gradients_cell = cell(num_experiments);
% meshes_cell = cell(num_experiments);
gradients = [];
meshes = [];
solutions = [];

for iter = 1:num_experiments
    
    diff = rand * (diffusion_coeff_range(2)-diffusion_coeff_range(1)); 
    diff_x = rand * (diffusion_coeff_range(2)-diffusion_coeff_range(1)); 
    advection = rand * (advection_coeff_range(2)-advection_coeff_range(1)); 
    advection_x = rand * (advection_coeff_range(2)-advection_coeff_range(1)); 
    reaction = rand * (reaction_coeff_range(2)-reaction_coeff_range(1)); 
    reaction_x = rand * (reaction_coeff_range(2)-reaction_coeff_range(1)); 
    rhs = rand * (rhs_range(2)-rhs_range(1));
    rhs_x = rand * (rhs_range(2)-rhs_range(1));
    
    [sol, original_gradient, final_mesh] = dynamic_mesh_refinement1D_2(n0, x0, x1, epsilon, diff, diff_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x);
%     gradients_cell{iter} = original_gradient; 
%     meshes_cell{iter} = final_mesh; 
    gradients = [gradients; original_gradient'];
    meshes = [meshes; final_mesh'];
    solutions = [solutions; sol']; 
    
end

dlmwrite('meshes.txt',meshes,'delimiter','\t')
dlmwrite('gradients.txt',gradients,'delimiter','\t')
dlmwrite('solutions.txt',solutions,'delimiter','\t')






