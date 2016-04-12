function [ X_t ] = ParticleFilter( X_t_1, z_t)
%Particle Filter
M = 25; % Number of particles
x_r = 30; % Measurement Noise covariance
x_n = 30; % state noise covariance
X_t = ones(M, 1);
P_w = zeros(M, 1);
X_t_u = ones(M, 1);

for m = 1:M

   if rand < 0.5
       X_t_u(m) = X_t_1(m) + x_n * rand; 
   else
       X_t_u(m) = X_t_1(m) - x_n * rand; 
   end

   P_w(m) = normpdf(X_t_u(m), z_t, x_r); % calculate weights
   
end
  for m = 1:M
    P_w(m) = P_w(m) + 50;
  end
  
  A = find(P_w > 50);
  S = size(A);
  i = 1;
  while i < S(1)
      P_w(A(i)) = P_w(A(i)) + 100*rand;
      i = i + 1;
  end
  P_w = P_w ./ sum(P_w);

  %Resampling
  for m = 1:M
     X_t(m) = X_t_u(find(rand < cumsum(P_w), 1)); 
  end
end