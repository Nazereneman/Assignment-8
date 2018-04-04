% 1D Diffusion Solver (du/dt=D*d^2u/dx^2+F(x,t))
% Using Crank-Nicholson method 
% Different spaced nodes in x and t
clc
clear;

%Problem 1
Nx = 10;        % N+1 node is at L, boundary
                % 1st node is on boundary at 0
Nt = Nx;        % This code assumes that the same number of nodes are used
                % in both x and t, but it can easily be changed
Lx = pi;        % x Domain for testing
Lt = 10;        % t Domain for testing, given as T
k = 1;			% k is given to be an integer, so choose 1
F = 0;          % F is given to be 0            
D = 0.1;        % D is given to be 0.1
g0 = 0;         % g0 is given as 0
gL=0;           % gL is given as 0
%f = sin(k*x);  % f is given
%Utnot = f;     % u(x,t=0) = f(x) is given and is defined below
Uxnot = g0;     % Dirichlet BC (u(x=0,t)=g0)
UxL = gL;       % Dirichlet BC (u(x=L,t)=gL)

era = 0.00001;% Used to compare itteration results
fprintf('The solution to the 1D using Crank-Nicolson:\n\n');
fprintf('Absolute error between itterations:   %1.1e\n\n', era);
fprintf('Number of equally spaced intervals in x:   %d\n', Nx);
fprintf('Number of equally spaced intervals in y:   %d\n', Nt);

% Define constants and RHS vector and apply BC to system
for m = 1:2
	Nx = m*Nx;              % Check to see difference between Nx and 2Nx
    Nt = m*Nt;              % Check to see difference between Ny and 2Ny
    max_er = 2*era;         % Initial error for while loop itteration stop condition
    it = 0;                 % iterations
    hx = Lx/Nx;             % Step size in x
    ht = Lt/Nt;             % Step size in y
   %Create constants for equations to reduce the # of operations in loop
   S_1 = -2*k;             % 
   S_2x = hx*k;             % f(x) is sin(k*x) so constant is k*x
   S_2y = ht*k;             % Constant in y
   %S_3x = hx*hx;            % multiply by .25 instead of divide by 4 (to better optomize) in x
   %S_3y = hy*hy;            % Same as above but in y  
   %S_4 = -2*S_3x-2*S_3y;    %Constant for if delta x and delta y are not the same
   %S_5=S_3x*S_3y;           %constant for f
    U = zeros(Nx+1,Nt+1);	% Create solution matix with initial guess of ZERO
    RHS = zeros(Nx+1,Nt+1);	% Create non homogeneous RHS vector
                            % Apply Dirichlet boundary conditions to solution matrix
    U(:,1) = Uxnot;         % LHS Boundary at x = 0 (LHS = left hand side)
  U(:,Nx+1) = UxL;          % RHS Boundary at x = L (RHS = right hand side)
    c=0;
  for j=2:Nt+1                % Bounday at t =0, given in problem
      c=c+hx;
      U(1,j) = sin(S_2x*c);    
  end  
end
%{
for i = 2:Nx             % Create RHS vector for all interior U, i,j = 1,N+1 on boundary		
        for j = 2:Nt 
       	RHS(i,j) = S_1*sin(S_2x*(j-1))*cosh(S_2y*(i-1));
        end
end

RHS = (S_5).*RHS;      %For Nx and Ny different than each other

%Do twice to compare N and 2N

while max_er > era          % RHS of Condition is percent error, no error in matrix can be larger than RHS

it = it + 1;                % Track the number of itterations with count
Uprev = U;                  % Must rewrite previous matrix so that can compare to new and determine error

for i = 2:Nx                % Gauss - Seidel itteration 
    for j = 2:Nt
		U(i,j) = (RHS(i,j)- S_3y.*U(i+1,j) -S_3y.* U(i-1,j) - S_3x.*U(i,j+1) - S_3x.*U(i,j-1) )/S_4;
		it_er(i,j) = abs( (U(i,j) - Uprev(i,j)));
    end
end
max_er = max(it_er(:));
end

for i=1:Nx+1
    for j = 2:Nt+1
    Uexact(i,j) = (Ly-(i-1)*ht)*sin(S_2x*(j-1))*sinh(S_2y*(i-1));
    end
end

for i = 2:Nx
    for j = 2:Nt
    L1_er(i,j) = abs( Uexact(i,j) - U(i,j));    
    end
end
if k == 2
	U2N = U;
    U2Ne= Uexact;
    L1err2N = (1/(Nx-1)^2) * sum(L1_er(:)); %Unsure how to do L1 error if using a different N for x and for y
    fprintf('Number of equally spaced intervals in x:   %d\n', Nx);
    fprintf('Number of equally spaced intervals in y:   %d\n', Nt);
    fprintf('Mean L1 error of all internal points: %1.5f\n\n', L1err2N);
else
    UN = U;
    UNe = Uexact;
    L1errN = (1/(Nx-1)^2) * sum(L1_er(:));
    fprintf('Mean L1 error of all internal points: %1.5f\n\n', L1errN);
end
end
% Plot and table data
x = linspace(0,pi,Nx+1);
figure;
plot(x,U2N(Nx/4+1,:),'g*',x,U2Ne(Nx/4+1,:),'g',x,U2N(Nx/2+1,:),'bo',x,U2Ne(Nx/2+1,:),'b',x,U2N(3*Nx/4+1,:),'rx',x,U2Ne(3*Nx/4+1,:),'r');
legend('U(x,\pi/4) approximate','U(x,\pi/4) exact','U(x,\pi/2) approximate','U(x,\pi/2) exact','U(x,3\pi/4) approximate','U(x,3\pi/4) exact');
figure(2);
y = linspace(0,pi,Nt+1);
plot(U2N(Nt/4+1,:),y,'g*',U2Ne(Nt/4+1,:),y,'g',U2N(Nt/2+1,:),y,'bo',U2Ne(Nt/2+1,:),y,'b',U2N(3*Nt/4+1,:),y,'rx',U2Ne(3*Nt/4+1,:),y,'r');
legend('U(\pi/4,y) approximate','U(\pi/4,y) exact','U(\pi/2,y) approximate','U(\pi/2,y) exact','U(3\pi/4,y) approximate','U(3\pi/4,y) exact');
%}
