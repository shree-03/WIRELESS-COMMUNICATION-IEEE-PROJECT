% WMMSE

clc;
clear all      
close all

% defining variables 
R = 2; % Number of antennas at each user
T = 3; % Number of antennas at each BS
K = 10; % Number of cells and Number of users 
epsilon=0.01;
max_iter_1=100;
iter_1=0;

SNRvec  = 5:5:30;
store_avg_sr=zeros(max_iter_1,length(SNRvec));

% calculating wsr for each SNR with a step of 5
  while iter_1<max_iter_1
      iter_1=iter_1+1;
    
    %channel is initialized as H(bs,user)--> H11,H12,H13,H14,H21,H22,H23...etc
    for ii=1:length(SNRvec)
          SNR = SNRvec(ii);
      % Generate channel between each BS and each user (K.K channels)
      H=[];
      for k= 1:K
         for i = 1:K
             H=[H channel(R,T)];
         end
      end

      alpha=ones(K,1);% user priority
      noise_power=1;%sigma squared 
      P=10^(SNR/10);%SNR in dB scale
      d=1; % Number of data streams for each user
      Ik=1;% single user per cell

  
     hold_wsr=[];
     mu_store=[];

     %initialization of matrices
     %assume Ik=1
     U=zeros(R,K*d)+1j*zeros(R,K*d);
     %initialize V with power constraint 
     for k=1:K
         V(:,k)=(H(1,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)');
     end
     for k = 1:K
         V(:,k) = sqrt(P/T)*V(:,k)/norm(V(:,k));
     end

     W=ones(d,K*d); 

     new_U=zeros(R,K*d)+1j*zeros(R,K*d);
     new_V=zeros(T,K*d)+1j*zeros(T,K*d);
     new_W=ones(d,K*d); 


     %initialization 
     iterations=0;
     BC=10; % break condition initialised to enter the loop
     max_iterations=500;
     iter_hold=[];

     % ---------ALGORITHM--------------------


     %  while loop begins
     while BC>epsilon && iterations<max_iterations
      iterations=iterations+1;

      iter_hold(end+1) = iterations; % storing the iterations in an array to plot
     
   
      % 1. receive beamformer updation
      s1=zeros(R,R*K);  % assume Ik=1
      for k=1:K % find U for each cell
            for j=1:K
                        s1(:,(k-1)*R+1:k*R)=s1(:,(k-1)*R+1:k*R) + H(:,(j-1)*T*K+(k-1)*T+1:(j-1)*T*K+k*T)*V(:,j)*(V(:,j)')*((H(:,(j-1)*T*K+(k-1)*T+1:(j-1)*T*K+k*T))');
            end         
            new_U(:,k)=inv(s1(:,(k-1)*R+1:k*R)+noise_power*eye(R))*H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)*V(:,k);
      end


     % 2. weights updation
      for k=1:K
             new_W(:,k)= real(inv(eye(1)-(new_U(:,k)')*H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)*V(:,k)));
      end


     % 3.  finding mu using bisection search
      mu_array=findmuh(H,new_U,new_W,K,P,T);
      mu_store=[mu_store mu_array];%for checking mu over all iterations


     %4. transmit beamformer updation 
     s2=zeros(T,T*K);
     for k=1:K
          for j=1:K
                     s2(:,(k-1)*T+1:k*T)=s2(:,(k-1)*T+1:k*T)+(alpha(j,:)*(H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T)')*new_U(:,j)*new_W(:,j)*(new_U(:,j)')*H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T));
           
          end       %---------------------------------
          new_V(:,k)=alpha(k,:)*inv(s2(:,(k-1)*T+1:k*T)+(mu_array(:,k)*eye(T)))*(H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)')*new_U(:,k)*new_W(:,k);
     end
         

     %5.Break condition calculation for termination of the loop 
     BC1=0;
     BC2=0;
     for j=1:K
                  BC1=BC1+log2(det(new_W(:,j)));
                  BC2=BC2+log2(det(W(:,j)));
     end
     BC=abs(BC1-BC2);


     %6.copying new values to old ones 
     U=new_U;
     V=new_V;
     W=new_W;   


     %7.Sum rate calculation
     wsr=compute_wsr(alpha,H,V,noise_power,K,T,R);
     hold_wsr=[hold_wsr wsr];% to put sum rate values in an array for plotting 


     end% end of while loop

   %snr and sum rate calculation
   avg_sum_rate=mean(hold_wsr);
   store_avg_sr(iter_1,ii)=avg_sum_rate;dlmwrite('YourOutputFile.txt', xy, 'delimiter', ',');
   %store_snr=[store_snr SNR]


   

   %----------PLOT-----------
   %PLOT :iterations vs Sumrate
   % figure
   % hold on 
   % grid on
   % plot(iter_hold,hold_wsr,'b');
   % plot(iter_hold,hold_wsr,'r*');
   % xlabel('number of iterations');
   % ylabel('sum rate');

 end  %End 


end
%----------PLOT-----------
%PLOT: SNR vs Average sum rate
figure
grid on;
plot(SNRvec,mean(store_avg_sr),'b');
%plot(store_snr,store_avg_sr,'r*');
xlabel('SNR');
ylabel('Average sum rate');


% -----FUNCTIONS-----------

%18 equation in the paper 
function result=muh_function(phi,lambda,muh,M,power)
        
         sum=0;
         for m=1:M
             sum=sum+(real(phi(m,m))/real((lambda(m,m))+muh)^2) ;
         end       
         result=sum-power;    
end 


% Sum Rate calculation function 
function result=compute_wsr(alpha,H,V,noise_power,K,T,R)
         
         s3=zeros(R,R*K);
         Rik=zeros(1,K);
         for k=1:K
                 for j=1:K           
                         if k~=j 
                            s3(:,(k-1)*R+1:k*R)=s3(:,(k-1)*R+1:k*R)  + H(:,(j-1)*T*K+(k-1)*T+1:(j-1)*T*K+k*T)*V(:,j)*(V(:,j)')*(H(:,(j-1)*T*K+(k-1)*T+1:(j-1)*T*K+k*T)');
                         end
                 end
                 Rik(:,k)=log2(real(det(eye(R)+(H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)*V(:,k)*(V(:,k)')*(H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)')*inv( s3(:,(k-1)*R+1:k*R)+ (noise_power*eye(R)) ) ) ) ));   
         end
         result=0;
         for k=1:K
                 result=result+alpha(k,:)*Rik(:,k);
         end

end 


% for finding mu using bisection search function
function mu_array=findmuh(H,U,W,K,power,T)

        mat=zeros(T,T*K); 
        phi=zeros(T,T*K);
        Vector=zeros(T,T*K);
        s4=zeros(T,T*K);% temporary variable to hold sum for the middle value of Phi calc
        mu_array=zeros(1,K);

        for k=1:K
                for j=1:K
                         % matrix that needs decomposition
                        mat(:,(k-1)*T+1:k*T)=mat(:,(k-1)*T+1:k*T)+(H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T)')*U(:,j)*W(:,j)*(U(:,j)')*H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T);
                end
        end
        for k=1:K  %---------------------------------------
               %eigen decomposition

               [Vector(:,(k-1)*T+1:k*T),lambda(:,(k-1)*T+1:k*T)]=eig(mat(:,(k-1)*T+1:k*T));
        end

        for k=1:K  
              %phi calculation intermediate steps 

              s4(:,(k-1)*T+1:k*T)=s4(:,(k-1)*T+1:k*T) + (H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)')*U(:,k)*(W(:,k)^2)*(U(:,k)')*H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T);
        end
        for k=1:K

              % phi calculation
              phi(:,(k-1)*T+1:k*T)=(Vector(:,(k-1)*T+1:k*T)')*s4(:,(k-1)*T+1:k*T)*(Vector(:,(k-1)*T+1:k*T));
        end
       
        for k=1:K
             
              %BISECTION SEARCH STARTS
              low=1e-10;
              high=10;
              iter=0;
              maxitr=5000;
              TOL=1e-10;

              A=trace(voptimum(0,H,U,W,T,K,k)*(voptimum(0,H,U,W,T,K,k)'));
              invertible_cond=vinvb(H,U,W,T,K,k);
              l=det(invertible_cond(:,(k-1)*T+1:k*T));

             if (A<=power) && (l)>1e-10
                  p=0;
                  mu_array(:,k)=p;
             else

               while iter<=maxitr
                 mid=(low+high)/2;
                 iter=iter+1;

                 guess=muh_function(phi(:,(k-1)*T+1:k*T),lambda(:,(k-1)*T+1:k*T),mid,T,power);
                 result1=muh_function(phi(:,(k-1)*T+1:k*T),lambda(:,(k-1)*T+1:k*T),low,T,power);
                    
               if guess==0 || abs(guess) < TOL
                  mu_array(:,k)=mid;
                  break

               elseif sign(guess)==sign(result1)
                          low=mid;
                          high=high;
               else
                          high=mid;
                          low=low;
          
               end
               end
             end
        end  
end 


%function used in bisection search
function V_opt=voptimum(mu,H,U,W,T,K,k)
          alpha=ones(K,1);
          V_opt=zeros(T,T*K);
          s5=zeros(T,T*K);
         
             for j=1:K
                     s5(:,(k-1)*T+1:k*T)=s5(:,(k-1)*T+1:k*T)+(alpha(j,:)*(H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T)')*U(:,j)*W(:,j)*(U(:,j)')*H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T));
           
             end       %---------------------------------
          V_opt(:,k)=alpha(k,:)*inv(s5(:,(k-1)*T+1:k*T)+(mu*eye(T)))*(H(:,(k-1)*T*K+(k-1)*T+1:(k-1)*T*K+k*T)')*U(:,k)*W(:,k);
 
      
end


%function used in bisection search
function invbmat=vinvb(H,U,W,T,K,k)
          alpha=ones(K,1);
          invbmat=zeros(T,T*K);
         
             for j=1:K
                     invbmat(:,(k-1)*T+1:k*T)=invbmat(:,(k-1)*T+1:k*T)+(alpha(j,:)*(H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T)')*U(:,j)*W(:,j)*(U(:,j)')*H(:,(k-1)*T*K+(j-1)*T+1:(k-1)*T*K+j*T));
             end      
       
end


%CHANNEL GENERATION FUNCTION
function result =channel(R,T)

         rp=(1/sqrt(2))*randn(R,T); % or try mean+sigma*randn(dim)
         imgp=(1/sqrt(2))*randn(R,T);
         result=rp+1j*imgp;
end













 
