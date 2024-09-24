function[A,B]=matrices(n)
% Generation of Hurwitz stable (I-C)^-1*(A+B) matrix
  % 0<=S<=1, density of sparse matrices 
  % S=1   =>  A, B all populated
  % S=0   =>  A, B are zero 
  % S=0.1 =>  A, B are strictly diagonally dominant

 S=zeros(1,2);
 S(1)=0.3+(0.5-0.3)*rand(1);
 S(2)=0.2+(0.3-0.2)*rand(1);
 
 A1 = full(sprandn(n,n,S(1)))/n; % Sparse non-shifted matrices
 B1 = full(sprandn(n,n,S(2)))/n;

 
 

 [max_radiusA] = gershdisc(A1);
 [max_radiusB] = gershdisc(B1);
 
                                  
 % Shift spectrum of A1 and B1
 % -------------------------------
          
 A = A1-max_radiusA*eye(n); % Sparse shifted matrices
 B = B1-max_radiusB*eye(n); %-(B1-max_radius*eye(n));

 k=2;    % larger k (>3) smaller DM, smaller k (<3)larger DM
 B=B*k;
 
 if 0==exist('ABmatrices','dir')       % check if Data folder already exists
    mkdir ABmatrices                   % create new directory to save the data
 end 
 

fname=sprintf('ABex%d.mat',n);
savdir = '.\ABmatrices';
save(fullfile(savdir,fname),'A','B');
 
 M1=A+B;
 ff=figure(1);
 subplot(2,1,1); box on
 
 
 plot(real(eig(M1)),imag(eig(M1)),'x','MarkerSize',5,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
 mM1=max(imag(eig(M1))); line([0 0],2*[-mM1-0.1*mM1-0.01 mM1+0.1*mM1+0.01],'LineStyle','--','Color','w');
 nMA=min(real(eig(M1))); xlim([nMA+0.1*nMA-0.01 -nMA-0.1*nMA+0.01]); ylim([-mM1-0.1*mM1-0.01 mM1+0.1*mM1+0.01]); 
 xlabel('R\lambda^{M(1)}'), ylabel('I\lambda^{M(1)}'), title('Eigenvalues of M(1)');
 set(gca,'Color',[0.1 0.1 0.1],'FontName','Times New Roman','FontSize',10);

 subplot(2,2,3)
 spy(A); xlabel('columns'); ylabel('rows'); title('Sparsity pattern of A');

 subplot(2,2,4)
 spy(B); xlabel('columns'); ylabel('rows'); title('Sparsity pattern of B');         
            
 
 
 set(ff,'units','normalized','outerpos',[0.6 0.07 0.4 0.5]);


 
 
 

function[max_radius]=gershdisc(A)
 error(nargchk(nargin,1,1));
 if size(A,1) ~= size(A,2)
    error('Matrix should be square');
    return;
 end

 radius=zeros(1,size(A,1));
 for i=1:size(A,1)
    % The circle has center in (h,k) where h is the real part of A(i,i) and
    % k is the imaginary part of A(i,i)   :
    % h=real(A(i,i)); k=imag(A(i,i)); 
 
    % Now we try to compute the radius of the circle, which is nothing more
    % than the sum of norm of the elements in the row where i != j
    r=0;
    for j=1:size(A,1)
       if i ~= j 
           r=r+(norm(A(i,j)));
       end    
    end 
    
    radius(1,i)=r;

    N=256;
    t=(0:N)*2*pi/N;
    
    % Now we're able to map each of the elements of this vector into a circle:
    % plot( r*cos(t)+h, r*sin(t)+k ,'-');

    % We also plot the center of the circle for better undesrtanding:
    % hold on;
    % plot( h, k,'+');
end

[max_radius]=max(radius);
%axis equal;
% hold on
%
% Now we plot the actual eigenvalues of the matrix:
% ev=eig(A);
% for i=1:size(ev)
%     rev=plot(real(ev(i)),imag(ev(i)),'.','MarkerSize',14,'MarkerFaceColor',[0.0 0.45 0.74],'MarkerEdgeColor',[0.0 0.45 0.74]);
%     grid on
% end
% 
% hold off;
%legend(rev,'Actual Eigenvalues');

