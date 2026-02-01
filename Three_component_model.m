%%
clear;
clc;
tic;

%% General parameters
Timer = 0; % s
Time_steps = [0.001, 0.01, 0.25]; % s - first, second, and third time step
End_times = [1, 10, 10000000]; % s - end time of first, second, and third period
Update_delay = 1000; % s
Num_pts = 150;
Side_length = 15;% mm
Point_dist = Side_length/Num_pts;

%% Adjustable parameters
%CSTR concentrations - mM
C0_S = 30; % Substrate
C0_A = 0; % Autocatalyst
C0_I = 1.5; % Fast inhibitor

%% Switch for Inhibitor concentration
Switch_time = 3000; % s
C0_I_end = 6.5; % New concentrition of Inhibitor

% Related to gel matrix
Ct_G = 300; % Binding matrix - constant total concentration - G + GA
Kd = 50; % [G]*[A]/[GA]

%Diffusion coefficients (mm^2/s)
DC_S = 4*10^(-4); 
DC_A = 4*10^(-4);
DC_I = 2*10^(-4);

%Reaction kinetics
k0 = 0.005; % 1/s 
k1 = 1*10^(-3); % 1/mM*s   S + A/GA -> A/GA + A
k2 = 0.3; % 1/mM*s A/GA + I -> x
k3 = 7.5*10^(-3); % 1/s   A/GA -> x
k4 = 3*10^(-5); % 1/s   S -> A/GA

%% Adjust initial concentrations
C_S = zeros(Num_pts);
C_A = zeros(Num_pts);
C_I = zeros(Num_pts);

% Tunable concentrations
C_S(:,:) = C0_S;
C_A(:,:) = C0_A;
C_I(:,:) = C0_I;

% % Tunable concentrations
% C_S(:,:) = 10.0001;
% C_A(:,:) = 0.8655;
% C_I(:,:) = 0.0154;

% Gel related
Sigma_numerator_matrix = zeros(Num_pts) + Ct_G*Kd;
Sigma_uppercase = (1 + Ct_G./(C_A + Kd));
Ct_A = C_A.*Sigma_uppercase;

%% Actication
% %% Manual activation
% C_A(end-20:end-14,27:37) = 1;
% C_I(end-20:end-14,27:37) = 0;

% %% Random activation
% RA_matrix = rand (Num_pts, Num_pts) > 0.90;
% C_A(RA_matrix) = 1;
% C_I(RA_matrix) = 0;

%% Activation template
Template_pic = imread("Template.jpg");
Template_matrix = imresize(Template_pic(:,:,1), [Num_pts, Num_pts]);
Template_matrix(Template_matrix<125) = 0;
Template_matrix = ~logical (Template_matrix);
% C_A(Template_matrix) = 0.85;
C_A(Template_matrix) = 1;
C_I(Template_matrix) = 0;

%% Calculting diffusion increment
IC_S = DC_S/((Point_dist)^2);
IC_A = DC_A/((Point_dist)^2);
IC_I = DC_I/((Point_dist)^2);
IC_check = Time_steps(end)*max([IC_S, IC_A, IC_I]); % - not more than 0.1!!!!!

%% Loading concentrations
% load('PatternB_300');
Ct_A = C_A.*Sigma_uppercase;

%% Plotting
h = figure('Name','Concentrations', 'Position', [0 0 1080 1080]);

subplot(2,2,1);
ha = image(C_S,'CDataMapping','scaled');
set(gca, 'YTickLabel', {}, 'Position', [0.02 0.58 0.47 0.4]);
colorbar
% clim([0 C0_S]);
xlabel('Substrate, mM');

subplot(2,2,2);
hb = image(C_A,'CDataMapping','scaled');
set(gca, 'YTickLabel', {}, 'Position', [0.52 0.58 0.47 0.4]);
colorbar
% clim([0 C0_S]);
xlabel('Free activator, mM');

subplot(2,2,3);
hc = image(C_I, 'CDataMapping','scaled');
set(gca, 'YTickLabel', {},'Position', [0.02 0.10 0.47 0.4]);
colorbar
xlabel('Fast Inhibitor, mM');

subplot(2,2,4);
hd = image(Ct_A, 'CDataMapping','scaled');
set(gca, 'YTickLabel', {},'Position', [0.52 0.10 0.47 0.4]);
colorbar
xlabel('Activator - Free + Bound, mM');
% 
time_str = strcat('Time passed:', 32, num2str(Timer), 32, 's');
h_annot = annotation('textbox', [0.85, 0.45, 0.1, 0.1], 'string', time_str);
% 
% g = figure('Name','Concentrations', 'Position', [0 0 1080 1080], 'color', 'white');
% 
% ga = image(Ct_A, 'CDataMapping','scaled');
% set(gca,'YTickLabel', {}, 'XTick', [0:(Num_pts/5):Num_pts], 'XTickLabel', {'0', '3', '6', '9', '12', '15'}, 'Position', [0.15 0.15 0.8 0.8]);
% colorbar
% % clim([0 25]);
% xlabel('Autocatalyst concentration, mM');
% g_annot = annotation('textbox', [0.8, 0.02, 0.1, 0.1], 'string', time_str);

%% Recording
Movie_rec = VideoWriter("Movies", 'MPEG-4');
Movie_rec.FrameRate = 24;
open(Movie_rec);
writeVideo(Movie_rec,getframe(h));

%% Calculation + Animation
Past_overall_time = 0; % s

for i=1:length(Time_steps)
      Current_time_step = Time_steps(i);
      Current_overall_time = End_times(i) - Past_overall_time;
      Current_step_number = Current_overall_time/Current_time_step;

      for j=1:Current_step_number
          
          Sigma_uppercase = (1 + Ct_G./(C_A + Kd));
          Sigma_lowercase = 1 + Sigma_numerator_matrix./((C_A + Kd).^2);
          Ct_A = C_A.*Sigma_uppercase;
          
          R_S_total = -k1*C_S.*Ct_A - k4*C_S;
          R_A_total = (k1*C_S - k2*C_I - k3).*Ct_A + k4*C_S;
          R_I_total = -k2*C_I.*Ct_A;
          
          C_S = C_S + Current_time_step*(diff_change(C_S, IC_S) + R_S_total + k0*(C0_S - C_S));
          C_A = C_A + Current_time_step*(diff_change(C_A, IC_A) + R_A_total - k0*C_A)./Sigma_lowercase;
          C_I = C_I + Current_time_step*(diff_change(C_I, IC_I) + R_I_total +  k0*(C0_I - C_I));
          
          Timer = Timer + Current_time_step;
          
        if mod(round(Timer,4), Update_delay) == 0
    
            ha.CData = C_S;
            hb.CData = C_A;
            hc.CData = C_I;
            hd.CData = Ct_A;
    
            time_str = strcat('Time passed:', 32, num2str(Timer), 32, 's');
            set(h_annot, 'String', time_str);

%             ga.CData = Ct_A;
%             set(g_annot, 'String', time_str);

            drawnow;
            writeVideo(Movie_rec,getframe(h));

            if round(Timer,4) == Switch_time
                C0_I = C0_I_end;
            end
    
        end
        
      end
      Past_overall_time = Past_overall_time + Current_overall_time;
end
close(Movie_rec);
saveas (h, "Last_frame.svg");
save('Concentrations', 'C_S', 'C_A', 'C_I');
toc;

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
