%%
clear;
clc;

%% Adjustable parameters
% General
Simulation_time = 35000; % In seconds - overall simulation time
Time_step = 0.25; % In seconds - used time step for updating concentrations (lower value gives more accurate results)
Update_delay = 100; % In seconds - frequency of figures update
Num_pts = 300; % Grid size (Num_pts x Num_pts)
Side_length = 15;% In mm - dimensions of reaction media (SIde_length x Side_length)

%Concentrations of the components in the CSTR (mM)
C0_S = 30; % Substrate
C0_A = 0; % Activator
C0_I = 5; % Fast inhibitor

% Switch of the fast inhibitor concentration
Switch_time = 9500; % In seconds - time when the concentration of Fast Inhibitor changes
C0_I_end = 7.8; % New concentrition for Fast Inhibitor

% Related to gel matrix
Ct_G = 300; % Concentration of activator-binding groups in hydrogel matrix - constant total concentration (G + GA) in mM
Kd = 50; % Dissosiation constant ([G]*[A]/[GA]) in mM

%Diffusion coefficients (mm^2/s)
DC_S = 4*10^(-4); % Sunstrate
DC_A = 4*10^(-4); % Activator
DC_I = 2*10^(-4); % Inhibitor

%Reaction kinetics
k0 = 0.0006; % 1/s - Supply rate of the componenets from the CSTR into the reaction media
k1 = 1*10^(-3); % 1/mM*s - Autocatalysis rate constant (S + A/GA -> A/GA + A)
k2 = 0.3; % 1/mM*s - Fast inhibition rate constant(A/GA + I -> x)
k3 = 5*10^(-3); % 1/s - Slow inhibition rate constant (A/GA -> x)
k4 = 3*10^(-5); % 1/s - Substrate hydrolysis rate constant  S -> A/GA

%% Non-adjustable parameters
% General
Timer = 0; % In seconds - elapsed time
Point_dist = Side_length/Num_pts; % In mm/pixel - 1 pixel represents x mm of reaction media
Time_steps = [0.001, 0.01, Time_step]; % s - first, second, and third time step (Initially, small time steps are required due to large concentration gradients)
End_times = [1, 10, Simulation_time]; % s - end time of first, second, and third period

% Initial concentrations
C_S = zeros(Num_pts);
C_A = zeros(Num_pts);
C_I = zeros(Num_pts);
C_S(:,:) = C0_S;
C_A(:,:) = C0_A;
C_I(:,:) = C0_I;

% Derived variables for concentration calculations
Sigma_numerator_matrix = zeros(Num_pts) + Ct_G*Kd;
Sigma_uppercase = (1 + Ct_G./(C_A + Kd));
Ct_A = C_A.*Sigma_uppercase; % Total concentration of activator - free and bound (in mM)

% Applying initial multi-spot activation through template
Template_pic = imread("Template.jpg");
Template_matrix = imresize(Template_pic(:,:,1), [Num_pts, Num_pts]);
Template_matrix(Template_matrix<125) = 0;
Template_matrix = ~logical (Template_matrix);
C_A(Template_matrix) = 1;
C_I(Template_matrix) = 0;

% Calculting diffusion increment
IC_S = DC_S/((Point_dist)^2); % For substrate
IC_A = DC_A/((Point_dist)^2); % For activator
IC_I = DC_I/((Point_dist)^2); % For fast inhibitor

%% Plotting
First_figure = figure('Name','Concentrations', 'Position', [0 0 1080 1080]);

subplot(2,2,1);
S_graph = image(C_S,'CDataMapping','scaled');
set(gca, 'YTickLabel', {}, 'Position', [0.02 0.58 0.47 0.4]);
colorbar
% clim([0 C0_S]);
xlabel('Substrate, mM');

subplot(2,2,2);
A_graph = image(C_A,'CDataMapping','scaled');
set(gca, 'YTickLabel', {}, 'Position', [0.52 0.58 0.47 0.4]);
colorbar
% clim([0 C0_S]);
xlabel('Free activator, mM');

subplot(2,2,3);
I_graph = image(C_I, 'CDataMapping','scaled');
set(gca, 'YTickLabel', {},'Position', [0.02 0.10 0.47 0.4]);
colorbar
xlabel('Fast Inhibitor, mM');

subplot(2,2,4);
FullA_graph = image(Ct_A, 'CDataMapping','scaled');
set(gca, 'YTickLabel', {},'Position', [0.52 0.10 0.47 0.4]);
colorbar
xlabel('Activator - Free + Bound, mM');

time_str = strcat('Time passed:', 32, num2str(Timer), 32, 's');
First_annot = annotation('textbox', [0.85, 0.45, 0.1, 0.1], 'string', time_str);

Second_figure = figure('Name','Concentrations', 'Position', [0 0 1080 1080], 'color', 'white');

FullA_graph_second = image(Ct_A, 'CDataMapping','scaled');
set(gca,'YTickLabel', {}, 'XTick', [0:(Num_pts/5):Num_pts], 'XTickLabel', {'0', '3', '6', '9', '12', '15'}, 'Position', [0.15 0.15 0.8 0.8]);
colorbar
% clim([0 25]);
xlabel('Total autocatalyst concentration, mM');
g_annot = annotation('textbox', [0.8, 0.02, 0.1, 0.1], 'string', time_str);

%% Recording
Movie_rec = VideoWriter("Movies", 'MPEG-4');
Movie_rec.FrameRate = 24;
open(Movie_rec);
writeVideo(Movie_rec,getframe(First_figure));

%% Calculation + Animation
Past_overall_time = 0; % s

for i=1:length(Time_steps)
      Current_time_step = Time_steps(i);
      Current_overall_time = End_times(i) - Past_overall_time;
      Current_step_number = Current_overall_time/Current_time_step;

      for j=1:Current_step_number
          % Updating variables
          Sigma_uppercase = (1 + Ct_G./(C_A + Kd));
          Sigma_lowercase = 1 + Sigma_numerator_matrix./((C_A + Kd).^2);
          Ct_A = C_A.*Sigma_uppercase;
          
          % Updating reaction rate for the components
          R_S_total = -k1*C_S.*Ct_A - k4*C_S;
          R_A_total = (k1*C_S - k2*C_I - k3).*Ct_A + k4*C_S;
          R_I_total = -k2*C_I.*Ct_A;
          
          % Updating concentration of components
          C_S = C_S + Current_time_step*(diff_change(C_S, IC_S) + R_S_total + k0*(C0_S - C_S));
          C_A = C_A + Current_time_step*(diff_change(C_A, IC_A) + R_A_total - k0*C_A)./Sigma_lowercase;
          C_I = C_I + Current_time_step*(diff_change(C_I, IC_I) + R_I_total +  k0*(C0_I - C_I));
          
          Timer = Timer + Current_time_step;
          
        % Updating graphs
        if mod(round(Timer,4), Update_delay) == 0
    
            S_graph.CData = C_S;
            A_graph.CData = C_A;
            I_graph.CData = C_I;
            FullA_graph.CData = Ct_A;
    
            time_str = strcat('Time passed:', 32, num2str(Timer), 32, 's');
            set(First_annot, 'String', time_str);

            FullA_graph_second.CData = Ct_A;
            set(g_annot, 'String', time_str);

            drawnow;
            writeVideo(Movie_rec,getframe(First_figure));

            if round(Timer,4) == Switch_time
                C0_I = C0_I_end; % Updating fast inhibitor concentration
            end
    
        end
        
      end
      Past_overall_time = Past_overall_time + Current_overall_time;
end
close(Movie_rec);

%% Function for calculating concentration gradients
function out = diff_change(C, incr_coeff) 
    RC = circshift(C,[ 0, 1]) - C;
    LC = circshift(C,[ 0,-1]) - C;
    UC = circshift(C,[ -1, [ 0, 1]]) - C;
    DC = circshift(C,[ 1, [ 0,-1]]) - C; 
    RC(:, 1) = 0;
    LC(:, end) = 0;
    UC(end, :) = 0;
    DC(1, :) = 0;
    out = (RC + LC + UC + DC).*incr_coeff;
end
