
function [CL_data, CD_data] = regularize_tables()
sd7032;

%===========================================================================
% for each Re. number, interpolate the tables to use regularly spaced AoA
%===========================================================================

global CL_dat CD_dat 
global alpha_min dalpha alpha_max
global Re_step Re_min Re_vec Re_max
global alpha_zero  
global GE_coeffs
alpha_min   = max([min(CL_data(:,1)),min(CD_data(:,1))]);
alpha_max   = min([max(CL_data(:,1)),max(CD_data(:,1))]);

dalpha      = 1;                       % resolution in alpha

Re_max      = 1000000;                  % upper Reynolds number limit for tables
Re_min      = 10000;                   % lower Reynolds number limit for tables
Re_step     = 10000;                   % step size in Reynolds number
num_Re      = floor(Re_max/Re_step);   % # intermediate Reynolds numbers
Re_vec      = Re_min:Re_step:Re_max;   % Reynolds number array

%===========================================================================
% make sure tables don't extend beyond known values
%===========================================================================

if(alpha_min <0)
    alpha_min = dalpha*(1+floor(alpha_min/dalpha));
else
    alpha_min = dalpha*(0+floor(alpha_min/dalpha));
end

if(alpha_max >0)
    alpha_max = dalpha*(0+floor(alpha_max/dalpha));
else
    alpha_max = dalpha*(1+floor(alpha_max/dalpha));
end

%===========================================================================
% construct angle of attack vector
%===========================================================================

alpha_vector = alpha_min:dalpha:alpha_max;

[a_lift,b_lift]=size(CL_data);
[a_drag,b_drag]=size(CD_data);

%===========================================================================
%Interpolating for CL
%===========================================================================

CL_matrix(:,1) = alpha_vector;
for i=1:length(alpha_vector)
    break_switch = 0;
    for j=1:a_lift-1                  % sweep over table aoa
        if(break_switch == 1)
            break
        end
        x_left = CL_data(j,1);
        x_right= CL_data(j+1,1);
        
% look for the current alpha in the data table
        if(alpha_vector(i) >= x_left && alpha_vector(i) <=x_right)
            break_switch = 1;
% if interval is found, no need to continue the loop
            for k=2:b_lift           
                y_left = CL_data(j,k);
                y_right= CL_data(j+1,k);
                CL_matrix(i,k) = y_left + (y_right-y_left)/(x_right-x_left)*(alpha_vector(i)-x_left);
            end
        end
    end
end

%===========================================================================
%Done interpolation for CL-alpha, preserving the original Reynolds number
%===========================================================================

%===========================================================================
%Interpolating for CD
%===========================================================================

CD_matrix(:,1) = alpha_vector;
for i=1:length(alpha_vector)
    break_switch = 0;
    for j=1:a_drag-1                  %sweep over table aoa
        if(break_switch == 1)
            break
        end
        x_left = CD_data(j,1);
        x_right= CD_data(j+1,1);
        
        %look for the current alpha in the data table
        if(alpha_vector(i) >= x_left && alpha_vector(i) <=x_right)
            break_switch = 1;
            %if interval is found, no need to continue the loop
            for k=2:b_drag           
                y_left = CD_data(j,k);
                y_right= CD_data(j+1,k);
                CD_matrix(i,k) = y_left + (y_right-y_left)/(x_right-x_left)*(alpha_vector(i)-x_left);
            end
        end
    end
end

%===========================================================================
%Now interpolate for Reynolds numbers
%===========================================================================

CL_dat = zeros(length(CL_matrix),num_Re+1);%initialize data array lift
CD_dat = zeros(length(CL_matrix),num_Re+1);%initialize data array drag

%===========================================================================
%Interpolating for CL for various Reynolds numbers
%===========================================================================

CL_dat(:,1) = alpha_vector;
for j=1:num_Re                  %sweep over table Re
    if(Re_vec(j) <= Re_data(1))
        CL_dat(:,j+1) = CL_matrix(:,2);
    elseif(Re_vec(j) > max(Re_data))
        CL_dat(:,j+1) = CL_matrix(:,end);
    else
%search for the interval where the required Reynolds number lies
        break_switch = 0;
        for i=1:length(Re_data)-1
            if(break_switch == 1)
                break
            end            
            left_Re    = Re_data(i);
            right_Re   = Re_data(i+1);
            left_array = CL_matrix(:,i+1);
            right_array= CL_matrix(:,i+2);
            if(Re_vec(j) >= left_Re && Re_vec(j) <= right_Re)
                break_switch = 1;
                CL_dat(:,j+1) = left_array + (right_array-left_array)/...
                                (right_Re-left_Re)*(Re_vec(j)-left_Re);
            end
        end
    end
end

%===========================================================================
%Interpolating for CD for various Reynolds numbers
%===========================================================================

CD_dat(:,1) = alpha_vector;
for j=1:num_Re                  %sweep over table Re
    if(Re_vec(j) <= Re_data(1))
        CD_dat(:,j+1) = CD_matrix(:,2);
    elseif(Re_vec(j) > max(Re_data))
        CD_dat(:,j+1) = CD_matrix(:,end);
    else
%search for the interval where the required Reynolds number lies
        break_switch = 0;
        for i=1:length(Re_data)-1
            if(break_switch == 1)
                break
            end
            left_Re    = Re_data(i);
            right_Re   = Re_data(i+1);
            left_array = CD_matrix(:,i+1);
            right_array= CD_matrix(:,i+2);
            if(Re_vec(j) >= left_Re && Re_vec(j) <= right_Re)
                break_switch = 1;
                CD_dat(:,j+1) = left_array + (right_array-left_array)./...
                                (right_Re-left_Re)*(Re_vec(j)-left_Re);
            end
        end
    end
end

%===========================================================================
%Converting aoa units from degrees to radians
%===========================================================================

alpha_vector = alpha_vector*pi/180;
alpha_min    =    alpha_min*pi/180;
alpha_max    =    alpha_max*pi/180;
dalpha       =       dalpha*pi/180;
alpha_zero   =   alpha_zero*pi/180;
CL_dat(:,1)  =  CL_dat(:,1)*pi/180;
CD_dat(:,1)  =  CD_dat(:,1)*pi/180;

return 