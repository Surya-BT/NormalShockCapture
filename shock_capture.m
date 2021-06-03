% solution of nozzle flow with shock
clear all
close all
clc

%% read analytical Soln
mach_analytical = readmatrix('mach_analyticalSoln.csv');
pressure_analytical = readmatrix('pressure_analyticalSoln.csv');

%% grid
nx = 61;
L = 3;
x = linspace(0,L,nx);
dx = x(2) - x(1);
% number of time-steps
nt = 1600;

gamma = 1.4;
% coefficient for artificial viscostiy
cx = 0.5;

%% initialisation
% variation of area
a = ones(size(x)) + (2.2).*(x - 1.5).^2;

% initial conditions for other variables
for i=1:length(x)
    
   if x(i)>= 0 && x(i)<=0.5 
       rho(i) = 1.0;
       t(i) = 1.0;
   elseif x(i)>= 0.5 && x(i)<=1.5
       rho(i) = 1.0 - 0.366*(x(i)-0.5);
       t(i) = 1.0 - 0.167*(x(i) - 0.5);
   elseif x(i)>= 1.5 && x(i)<=2.1 
       rho(i) = 0.634 - 0.702*(x(i)-1.5);
       t(i) = 0.833 - 0.4908*(x(i) - 1.5);
   else
       rho(i) = 0.5892 + 0.10228*(x(i)-2.1);
       t(i) = 0.93968 - 0.0622*(x(i) - 2.1);       
   end
   
   v(i) = 0.59/(rho(i)*a(i));
   p(i) = rho(i)*t(i);
   
   u1(i) = rho(i) * a(i);
   u2(i) = rho(i) * a(i) * v(i);
   u3(i) = rho(i) * (t(i)/(gamma-1) + gamma/2*(v(i)^2)) * a(i);
end
p(nx) = 0.6784; % forcing the exit pressure to a constant value

dt = 0.01075;

%% solution - MacCormack method
% time loop
for k=1:nt
    
    u1_old = u1;
    u2_old = u2;
    u3_old = u3;
  
    for i = 1:nx
        f1(i) = u2(i);
        f2(i) = (u2(i)^2)/u1(i) + (gamma -1)/gamma*(u3(i) - gamma/2*(u2(i)^2)/u1(i));
        f3(i) = gamma*u2(i)*u3(i)/(u1(i)) - gamma*(gamma-1)/2*u2(i)^3/(u1(i)^2);
    end
    % predictor loop
    for j = 2:nx-1
        
       df1dx = (f1(j+1) - f1(j))/dx;
       df2dx = (f2(j+1) - f2(j))/dx;
       df3dx = (f3(j+1) - f3(j))/dx;
       dadx(j) = (a(j+1) - a(j))/dx; 
       j2(j) = 1/gamma*rho(j)*t(j)*dadx(j); 
       
       du1dt_p(j) = -1*df1dx;
       du2dt_p(j) = -1*df2dx + j2(j);
       du3dt_p(j) = -1*df3dx;
       
       % calculation of artificial viscosity
       d2pdx2(j) = p(j+1) - 2*p(j) + p(j-1);
       d2pdx2_d(j) = p(j+1) + 2*p(j) + p(j-1);
       d2u1dx2(j) = u1(j+1) - 2*u1(j) + u1(j-1);
       d2u2dx2(j) = u2(j+1) - 2*u2(j) + u2(j-1);
       d2u3dx2(j) = u3(j+1) - 2*u3(j) + u3(j-1);
       
       s1(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u1dx2(j);
       s2(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u2dx2(j);
       s3(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u3dx2(j);
       
       % solution update - primary variables
       u1(j) = u1(j) + (du1dt_p(j) * dt) + s1(j);
       u2(j) = u2(j) + (du2dt_p(j) * dt) + s2(j);
       u3(j) = u3(j) + (du3dt_p(j) * dt) + s3(j);
       
       % solution update - primitive variables
       rho(j) = u1(j)/a(j);
       t(j) = (gamma-1)*(u3(j)/u1(j) - (gamma/2)*(u2(j)/u1(j))^2);
       p(j) = rho(j)*t(j);
       
              
    end
    
    % boundary values
    u2(1) = 2*u2(2) - u2(3);
    v(1) = u2(1)/u1(1);
    u3(1) = u1(1) * ( t(1)/(gamma-1) + gamma/2*(v(1)^2));
    
    u1(nx) = 2*u1(nx-1) - u1(nx-2);
    u2(nx) = 2*u2(nx-1) - u2(nx-2);
    v(nx) = u2(nx)/u1(nx);
    p(nx) = 0.6784;
    u3(nx) = (0.6784*a(nx))/(gamma-1) + (gamma/2)*u2(nx)*v(nx);
    
    
    % solution update
    for i = 1:nx
        f1(i) = u2(i);
        f2(i) = (u2(i)^2)/u1(i) + (gamma -1)/gamma*(u3(i) - gamma/2*(u2(i)^2)/u1(i));
        f3(i) = gamma*u2(i)*u3(i)/(u1(i)) - gamma*(gamma-1)/2*u2(i)^3/(u1(i)^2);
    end
    
    % correction step
    for j = 2:nx-1
        
       df1dx = (f1(j) - f1(j-1))/dx;
       df2dx = (f2(j) - f2(j-1))/dx;
       df3dx = (f3(j) - f3(j-1))/dx;
       dadx(j) = (a(j) - a(j-1))/dx; 
       j2(j) = 1/gamma*rho(j)*t(j)*dadx(j); 
       
       du1dt_c(j) = -1*df1dx;
       du2dt_c(j) = -1*df2dx + j2(j);
       du3dt_c(j) = -1*df3dx;
       
    end
    % calculating the average value
    du1dt_avg = 0.5*(du1dt_p + du1dt_c);
    du2dt_avg = 0.5*(du2dt_p + du2dt_c);
    du3dt_avg = 0.5*(du3dt_p + du3dt_c);
    
    for j=2:nx-1
        
       % calculation of artificial viscosity
       d2pdx2(j) = p(j+1) - 2*p(j) + p(j-1);
       d2pdx2_d(j) = p(j+1) + 2*p(j) + p(j-1);
       d2u1dx2(j) = u1(j+1) - 2*u1(j) + u1(j-1);
       d2u2dx2(j) = u2(j+1) - 2*u2(j) + u2(j-1);
       d2u3dx2(j) = u3(j+1) - 2*u3(j) + u3(j-1);
       
       s1(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u1dx2(j);
       s2(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u2dx2(j);
       s3(j) = cx * (abs(d2pdx2(j))/d2pdx2_d(j)) * d2u3dx2(j);
               
       u1(j) = u1_old(j) + (du1dt_avg(j) * dt) + s1(j);
       u2(j) = u2_old(j) + (du2dt_avg(j) * dt) + s2(j);
       u3(j) = u3_old(j) + (du3dt_avg(j) * dt) + s3(j);
       
    end

    % boundary values
    u2(1) = 2*u2(2) - u2(3);
    v(1) = u2(1)/u1(1);
    u3(1) = u1(1) * ( t(1)/(gamma-1) + gamma/2*(v(1)^2));
    
    u1(nx) = 2*u1(nx-1) - u1(nx-2);
    u2(nx) = 2*u2(nx-1) - u2(nx-2);
    v(nx) = u2(nx)/u1(nx);
    p(nx) = 0.6784;
    u3(nx) = (0.6784*a(nx))/(gamma-1) + (gamma/2)*u2(nx)*v(nx);
    
    % solution variables
    for j=1:nx
       rho(j) = u1(j)/a(j);
       t(j) = (gamma-1)*(u3(j)/u1(j) - (gamma/2)*(u2(j)/u1(j))^2);
       v(j) = u2(j)/u1(j);
       p(j) = rho(j)*t(j);
       M(j) = v(j)/sqrt(t(j));
        
       % storing values for plotting
       rho_plot(j,k) = rho(j);
       v_plot(j,k) = v(j);
       t_plot(j,k) = t(j);
       p_plot(j,k) = p(j);
       m_plot(j,k) = u2(j);
       M_plot(j,k) = M(j);
    end
    % plot to check the progression of the solution
     plot(x,p)
     ylim([0,1]);
     title(sprintf("time-step=%d",k));
     xlabel("nozzle length in [m]");
     ylabel("pressure in [Pa]");
     pause(0.03)
    
end

%% variables for quantitative validation
throat_loc = find(a==1)
disp("throat")
disp("_____")
disp("")
rho(throat_loc)
v(throat_loc)
t(throat_loc)
p(throat_loc)
M(throat_loc)
m_plot(throat_loc,nt)

disp("exit")
disp("_____")
disp("")
rho(end)
v(end)
t(end)
p(end)
M(end)
m_plot(end,nt)

%% analytical solution for shock wave
% shock location - from the manual calculation
shock_loc = 2.1;

for i = 1:length(x)
    if x(i) < shock_loc
        
        if x(i) < 1.5
            flow_type = 'subsonic';
            mach_num_calc(i) = calcMachNum(a(i),flow_type);
            %disp('sub')
        elseif x(i) > 1.5
            flow_type = 'supersonic';
            mach_num_calc(i) = calcMachNum(a(i),flow_type);
            %disp('sup');
        else
            mach_num_calc(i) = 1;
        end
    else
        flow_type = 'subsonic';
        mach_num_calc(i) = calcMachNum(a(i),flow_type);
    end
    
end

plot(x,mach_num_calc)

%% numerical vs analytical plot
% Uncomment this section, if you do not want to save a video file
% for k = 1:nt
%     subplot(2,1,1)
%     %plot(x,mach_num_calc,'--r');
%     plot(pressure_analytical(:,1),pressure_analytical(:,2),'--r')
%     hold on;
%     plot(x,p_plot_noArt(:,k),'-b');
%     title(sprintf("time-step=%d\n Artificial Viscosity Coefficient = 0.0",k));
%     legend("analytical Solution","Numerical Solution",'Location','southwest');
%     hold off;
%     pause(0.03)
%     
%     subplot(2,1,2)
%     plot(pressure_analytical(:,1),pressure_analytical(:,2),'--r')
%     hold on;
%     plot(x,p_plot(:,k),'-b');
%     title(sprintf("Artificial Viscosity Coefficient = 0.2"));
%     hold off
%     pause(0.03)
% end

%% write video file
filename = "shockCapture.avi";
videoFile = VideoWriter(filename);
videoFile.FrameRate = 30;
open(videoFile)
for k = 1:nt
    subplot(2,1,1)
    %plot(x,mach_num_calc,'--r');
    plot(pressure_analytical(:,1),pressure_analytical(:,2),'--r')
    hold on;
    plot(x,p_plot_noArt(:,k),'-b');
    title(sprintf("time-step=%d\nDistribution of Pressure along the nozzle length\nArtificial Viscosity Coefficient = 0.0",k));
    legend("analytical Solution","Numerical Solution",'Location','southwest');
    xlabel("Non-dimensional length of the nozzle (x/L) [-]");
    ylabel("p/p_0 [-]");
    xlim([0,3])
    hold off;
    pause(0.003)
    
    subplot(2,1,2)
    plot(pressure_analytical(:,1),pressure_analytical(:,2),'--r')
    hold on;
    plot(x,p_plot(:,k),'-b');
    title(sprintf("Distribution of Pressure along the nozzle length\nArtificial Viscosity Coefficient = 0.2"));
    legend("analytical Solution","Numerical Solution",'Location','southwest');
    xlabel("Non-dimensional length of the nozzle (x/L) [-]");
    ylabel("p/p_0 [-]");
    xlim([0,3])
    hold off
    pause(0.003)
    frame = getframe(gcf);
    writeVideo(videoFile,frame);
 end
close(videoFile)

%% functions
% to calculate the area of the nozzle from the mach number
function mach_num = calcMachNum(AR,flow_type)

    g = 1.4;
    gp1 = g + 1;
    gm1 = g -1;
    %ARatio = 5;
    prob.objective = @(M) (1/M^2)*(((2 + gm1*M^2)/gp1)^(gp1/gm1)) - AR^2;
    prob.solver = 'fzero';
    prob.options = optimset(@fzero);
    prob.x0 = [1e-6 1];
    Msub = fzero(prob); 
    prob.x0 = [1+1e-6 50];
    Msup = fzero(prob);
    switch flow_type
        case 'subsonic'
            mach_num = Msub;
            %disp('sub');
        case 'supersonic'
            mach_num = Msup;
            %disp('sup');
        otherwise
            fprintf('invalid option');
    end
end