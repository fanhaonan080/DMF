function [alphaCRIT,lambdaCRIT] = intersections2(a,theta1,theta2,lambda1,lambda2,A,B,tauinh)
% Increasing accuracy of the zero-crossings
n = size(A,1);
iteration = 50;
tol = a;
cnt = 1;
flag1 = 0;
flag2 = 0;

% Zero-crossings
%--------------------------------------------------------------------------
    for i = 1:1
        if objectiveFun(lambda1(i,1),tauinh)*objectiveFun(lambda2(i,1),tauinh) < 0
            cnt2 = i;
            break
        end
    end
    f2s = lambda1(cnt2,1);
    f2e = lambda2(cnt2,1);

% Iterations
%--------------------------------------------------------------------------
    while cnt<=iteration
          p = theta1+(theta2-theta1)/2;
          tempF2 = eig( A+B*exp(-sqrt(-1)*p) ); % using the upper semi-circle: alpha + j*sqrt(1-alpha^2)
          [~,iii] = min(abs(real(tempF2)-(real(f2e)+real(f2s))/2));
          lambdaCRIT = tempF2(iii);

          if abs(theta2-theta1)/2 <= tol
             flag1 = 1;
          end
          if cnt==iteration
             flag2 = 1;
          end
          if (flag1==1) || (flag2==1)
              alphaCRIT = p;
              break
          end
          cnt = cnt+1;
          if f2s*real(lambdaCRIT) < 0
             theta2 = p;
             f2e = lambdaCRIT;
          end
          if real(lambdaCRIT)*f2e < 0
             theta1 = p;
             f2s = lambdaCRIT;
          end
    end
end

function obj = objectiveFun(FamLambda,tauinh)
ReFamLambda = real(FamLambda);
ImFamLambda = imag(FamLambda);
omega = sqrt(ReFamLambda.^2+ImFamLambda.^2);

% obj1 = (ReFamLambda./ImFamLambda)+tan(omega*tauinh);
obj = ReFamLambda.*cos(omega*tauinh) + ImFamLambda.*sin(omega*tauinh);
% obj3 = ReFamLambda + omega.*sin(omega*tauinh);
% obj4 = ImFamLambda - omega.*cos(omega*tauinh);
end