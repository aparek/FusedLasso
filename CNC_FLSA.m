function [ x, cost] = CNC_FLSA( y, lam0, lam1, a0, a1, Nit, pen)
% function [ x, cost] = CNC_FLSA( y, lam0, lam1, a0, a1, Nit, pen)
% 
%     This function minimizes the CNC FLSA objective function
%     
%     F(x) = 0.5||y-x||_2^2 + lam0*phi(x,a0) + lam1*phi(Dx,a1)
%     
%     using the majorization-minimization technique, where phi 
%     is a non-convex penalty function. 
%     
%     Input: 
%         y           - Noisy input signal      
%         lam0,lam1   - Regularization parameters         
%         a0,a1       - non-convex penalty parameters 
%                     (0 <= a0*lam0 + 4*a1*lam1 <= 1)                     
%         Nit         - Number of iterations         
%         pen         - Non-convex penalty to be used
%                     (select from 'L1','log','atan')
%     Output:
%         x           - Denoised output signal         
%         cost        - Cost function history          
%             
%     Ankit Parekh, NYU School of Engineering.
%     Ref.:   Convex fused lasso denoising with non-convex regularization 
%             and its
%             use for pulse detection.
%             Ankit Parekh and Ivan W. Selesnick. 2015. 
%          

% Function definitions
Fl = @(x,lam1,lam2,N) soft(tvd(x,N,lam2),lam1);                             % L1 FLSA two-step solution

% Select the non-convex penalty to be used
switch pen
    case 'L1'        
        phi = @(x) abs(x);                                                  %phi = phi(x;a)
        ds = @(x) zeros(size(x));                                           %ds = s'(x;a)
        
    case 'log'        
        phi = @(x,a) 1/a * log(1 + a*abs(x));
        ds = @(x,a) - a * x ./ (1 + a*abs(x));
            
    case 'atan'        
        phi = @(x,a) 2./(a*sqrt(3)) .* (atan((2*a.*abs(x)+1)/sqrt(3)) - pi/6);
        dphi = @(x,a) 1 ./(1 + a*abs(x) + a.^2.*abs(x).^2) .* sign(x);
        ds = @(x,a) dphi(x,a) - sign(x);
    otherwise
        disp('Please select a non-convex penalty function from - L1,log,atan')
        return
end

%Initialize
cost = zeros(1, Nit);                                                       % Cost function history
y = y(:).';                   
N = length(y);
D = @(x) diff(x);                                                           % First order difference matrix
DT = @(x) [-x(1) -diff(x) x(end)];                                          % Transpose of the matrix D
x = Fl(y,lam0,lam1,N);                                                      % Initialize the algorithm with L1 FLSA solution
Dx = D(x);

% Run the MM algorithm for Nit number of iterations
for k = 1:Nit
    x = Fl(y - lam1*DT(ds(Dx,a1))-lam0*ds(x,a0),lam0,lam1,N);
    Dx = D(x);
    cost(k) = lam0 * sum(phi(x, a0)) +...
                lam1 * sum(phi(Dx, a1)) + 0.5 * sum(abs(x(:)-y(:)).^2);
end

end

