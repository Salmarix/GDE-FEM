function [ n_FEM,tt,comp_time ] = TimeEvolutionGDE( M,K,B,C,n_init,t_interval,h,C_in,dC_in,source,method )
%TimeEvolutionGDE Solves FEM matrix equation for whole GDE with
% Euler, Implicit Euler and with Crank-Nicolson.

%   INPUTS:
%   M : mass matrix (matrix in the left hand side of FEM equation)
%   K : Stiffness matrix (condensation and removal term)
%   B : Coagulation formation matrix (set input as empty if no
%   coagulation) CELL ARRAY
%   C : Coagulation removal matrix
%   n_init : Initial size distribution for the problem
%   t_interval : first and last time index in the vector [0,tmax]
%   h : time step for Crank-Nicolson
%   C_in : Border condition
%   dC_in : derivative of border condition
%   source : source term (independent of size distribution function)
%   method : decision for method (Defalut Crank-Nicolson method)
%
%   OUTPUTS
%   n_FEM : time evolution of size distribution
%   tt : time vector for the evolution
%
% Teemu Salminen
% Department of Technical Physics
% 2020
%% Intitializations for FEM and time evolution
tt = t_interval(1):h:t_interval(2);
tt = tt(:);
n_original = n_init;
n_FEM = zeros(length(n_original),length(tt));

% If boundary conditions exists
if ~isempty(C_in)
    
    % Applying boundary conditions
    M(1,:) = [];
    K(1,:) = [];
    M_b = M(:,1);
    K_b = K(:,1);
    M(:,1) = [];
    K(:,1) = [];
    K_init = K;
    n_init(1) = [];
    n_FEM(2:end,1) = n_init;
    
else
    
    % Without boundary condition
    n_FEM(:,1) = n_init;
    K_init = K;
    
end


% Setting Crank-Nicolson to be default
if nargin < 11
    
    method = 'crank';
    
end

% If Coagulation is present
if ~isempty(B) && ~isempty(C)
    
    for jj = 1:length(B)
        
        if ~isempty(C_in)
            
            % Applying boundary conditions for the FEM matrices
            BB = B{jj};
            CC = C{jj};
            BB(1,:) = [];
            BB(:,1) = [];
            CC(1,:) = [];
            CC(:,1) = [];
            
            % Formation of coagulation FEM matrix with boundary condition
            if jj > 1
                
                BC{jj-1,1} = BB-CC;
                
            end
            
        else
            
            % Formation of coagulation FEM matrix without boundary
            % conditions
            
            BC{jj,1} = B{jj}-C{jj};
            
        end
        
    end
    
else
    
    
    
    
end

% Calculating time evolution for FEM
while ~strcmp(method,'euler') && ~strcmp(method,'impl.euler')  && ~strcmp(method,'crank')
    prompt = 'Unknown method as input. Set method to be either euler, impl.euler or crank: \n';
    method = input(prompt,'s');
    
    
    
end
tic
for jj = 1:length(tt)-1
    
    if strcmp(method,'crank')
        % Check for if coagulation is present
        if exist ('BC','var')
            % Set K matrix to empty
            K = [];
            HH = [];
            
            % Linearizing cougulation FEM matrix with or without boundary
            % conditions
            if ~isempty(C_in)
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(2:end,jj)'*BC{ii};
                    
                end
            else
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(:,jj)'*BC{ii};
                    
                end
                
            end
            
            % Forming new stiffness matrix
            K = K_init+HH;
            
        end
        % Left Hand side for Crank-Nicolson method
        S_l = [eye(size(M))-(h/2).*(M\(K))];
        % Right hand side for Crank-Nicolson (This part alway exists)
        S_r = [eye(size(M))+(h/2).*(M\(K))];
        
        % If Source term is included
        if ~isempty(source) && ~isempty(C_in)
            S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))))+source(tt(jj+1));
            S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))))+source(tt(jj));
            
            % Without source term
        elseif ~isempty(C_in) && isempty(source)
            S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))));
            S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))));
            
            % Only source term but no boundary conditions
        elseif isempty(C_in) && ~isempty(source)
            S2= source(tt(jj+1));
            S1= source(tt(jj));
        end
        if exist('S1','var')
            n_FEM(2:end,jj+1) = S_l\S_r*n_FEM(2:end,jj) + (h/2)*S_l\(S1+S2);
            n_FEM(1,jj+1) = C_in(tt(jj+1));

            % Positivity constraint
            n_FEM(n_FEM < 0) = 0;
        else
            n_FEM(:,jj+1) = S_l\(S_r)*n_FEM(:,jj);
        end
        
        
        % If the Euler method is chosen
    elseif strcmp(method,'euler')
        
        % Check for if coagulation is present
        if exist ('BC','var')
            % Set K matrix to empty
            K = [];
            HH = [];
            
            % Linearizing cougulation FEM matrix with or without boundary
            % conditions
            if ~isempty(C_in)
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(2:end,jj)'*BC{ii};
                    
                end
            else
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(:,jj)'*BC{ii};
                    
                end
                
            end
            
            % Forming new stiffness matrix
            K = K_init+HH;
            
        end
        
        % If Source term is included
        if ~isempty(source) &&  ~isempty(C_in)
            S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))))+source(tt(jj));
            
            % Without source term
        elseif ~isempty(C_in) && isempty(source)
            S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))));
            % Only source term but no boundary conditions
        elseif isempty(C_in) && ~isempty(source)
            S1= source(tt(jj));
        end
        if exist('S1','var')
            n_FEM(2:end,jj+1) = n_FEM(2:end,jj) + (h)*(M\(K)*n_FEM(2:end,jj)+S1);
            n_FEM(1,jj+1) = C_in(tt(jj+1));
        else
            n_FEM(:,jj+1) = n_FEM(:,jj) + (h)*(M\(K*n_FEM(:,jj)));
        end
        
        % Implicit Euler time evolution calculation
    elseif strcmp(method,'impl.euler')
        
        % Check for if coagulation is present
        if exist ('BC','var')
            % Set K matrix to empty
            K = [];
            HH = [];
            
            % Linearizing cougulation FEM matrix with or without boundary
            % conditions
            if ~isempty(C_in)
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(2:end,jj)'*BC{ii};
                    
                end
            else
                
                for ii = 1:length(BC)
                    
                    HH(ii,:) = n_FEM(:,jj)'*BC{ii};
                    
                end
                
            end
            
            % Forming new stiffness matrix
            K = K_init+HH;
            
        end
        % Left Hand side for Crank-Nicolson method
        S_l = [eye(size(M))-(h).*(M\(K))];
        
        % If Source term is included
        if ~isempty(source) && ~isempty(C_in)
            S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))))+source(tt(jj+1));
            
            % Without source term
        elseif ~isempty(C_in) && isempty(source)
            S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))));
            
            % Only source term but no boundary conditions
        elseif isempty(C_in) && ~isempty(source)
            S2= source(tt(jj+1));
        end
        if exist('S2','var')
            n_FEM(2:end,jj+1) = S_l\(n_FEM(2:end,jj) + (h)*S2);
            n_FEM(1,jj+1) = C_in(tt(jj+1));
        else
            n_FEM(:,jj+1) = S_l\n_FEM(:,jj);
        end
        
    end
    
    
end
comp_time = toc;
end

