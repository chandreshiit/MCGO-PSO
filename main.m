final_mean=0;
for fm=1:10
Data=csvread('dataset/chess.csv',1,0);
% n1=ceil(.6*8774);
% ntotal=8774;
% k=7;
% d=524;
% c1=1.5;
% c2=.5;
k=3;
n1=433;
ntotal=533;
d=6;
c1=1.2;
c2=2;

Classes=unique(Data(:,d+2));
% k=7;
% n=50;
% d=524;
%% Parameters Initialization
pop=70;
gen=30;
tic;
prtcl_indx=zeros(1,pop);
psomatrix_x=zeros(k,d,pop);
particle_store_x=zeros(k,d,pop);
psomatrix_v=zeros(k,d,pop);
particle_store_v=zeros(k,d,pop);
 particle_store_Fitness=zeros(1,pop);
 pbest=particle_store_x;
  for j=1:pop  
  psomatrix_x(:,:,j)=rand(k,d); 
  psomatrix_v(:,:,j)=((2*rand(1,1))-1); %NEED TO CORRECT
 
 [fitnessVal] = Fitness( Data, psomatrix_x(:,:,j),k,d,n1); % Fitness
%  if (f_val<= fitnessVal)
%        f_val= fitnessVal;
%  end
 particle_store_Fitness(1,j)=fitnessVal; %store fitness
 particle_store_x(:,:,j)=psomatrix_x(:,:,j);
  particle_store_v(:,:,j)=psomatrix_v(:,:,j);
 pbest(:,:,j)= particle_store_x(:,:,j); %initialize Pbest
 end
%INITIALIZATION DONE..

%% Best Fitness Initialization
  [best_fitness pind]=min(particle_store_Fitness); 
  %FINAL_MWT1=best_fitness; 
  gbest=particle_store_x(:,:,pind); 
  gbest_node_sequence=particle_store_x(1,:,pind);
  gbest_prtcl_end_indx=prtcl_indx(pind);
   best_fitness1= best_fitness
%% Iterations

   for k1=1:gen
    display(k1)
     pbest1=pbest;
     gbest_prtcl_end_indx1=gbest_prtcl_end_indx;
     gbest1=gbest;
     %gbest_node_sequence1=gbest_node_sequence;
    particle_store_Fitness1= particle_store_Fitness;
     psomatrix1_x=particle_store_x;
     psomatrix1_v=particle_store_v;
     w=0.9*rand(1,1);
    
 
     r1=rand(1,1);
     r2=rand(1,1);
    prtcl_indx=zeros(1,pop);
    for i=1:pop
        %    w=w*rand(1,1);  
  psomatrix_v(:,:,i)=w*psomatrix1_v(:,:,i)+((c1*r1)*( pbest1(:,:,i)-psomatrix1_x(:,:,i)))+((c2*r2)*( gbest1(:,:,1)-psomatrix1_x(:,:,i))); 
  A = abs(psomatrix_v(:,:,i));
  [max_value,index] = max(A(:));
  psomatrix_v(:,:,i) = psomatrix_v(:,:,i)/max_value;
  [m,n] = size(psomatrix_v(:,:,i));

  psomatrix_x(:,:,i)=psomatrix1_x(:,:,i)+psomatrix_v(:,:,i);
  
  psomatrix_x(:,:,i)=mod(psomatrix_x(:,:,i),1);
  
%   [max_value,index1] = max(A(:));
%   [min_value,index2] = min(A(:));
%   psomatrix_x(:,:,i) = psomatrix_x(:,:,i)-min_value/(max_value-min_value);
%  psomatrix_x(:,:,i)

[fitnessVal] = Fitness(Data, psomatrix_x(:,:,i),k,d,n1 ); % Fitness
 %if (particle_store_Fitness(1,i)<= fitnessVal)
  particle_store_Fitness(1,i)= fitnessVal;
  particle_store_x(:,:,i)=psomatrix_x(:,:,i);
  particle_store_v(:,:,i)=psomatrix_v(:,:,i);
% end
 
 
  %particle_store_Fitness(1,i)=fitnessVal;
 if( particle_store_Fitness(1,i)< particle_store_Fitness1(1,i))
     pbest(:,:,i)= particle_store_x(:,:,i);
  end
    end
    
    [best_fit, pind]=min(particle_store_Fitness);
  if(best_fit < best_fitness1)
    gbest(:,:,1)=particle_store_x(:,:,pind);
     gbest_prtcl_end_indx=prtcl_indx(pind);
     best_fitness1 = best_fit;
     
  end

%   s(k1)=best_fitness1;
%   r(k1)=best_fit;
   %% Plotting the swarm
%    { clf    
%     plot(k1,best_fitness1,'x')   % drawing swarm movements
%     drawnow
%     title('Plot of Global Fitness Vs No. of Iterations')
%     axis([-2 gen -2 1000]);
% pause(.2)
   end
t=0;
r=zeros(1,k);
zu=zeros(1,k);
i=n1;
jp=ntotal-n1+1;
% jp=n1;
for t=1:jp
    
    x=Data(i,2:d+1);
    y=Data(i,d+2);
    p(i)=Data(i,d+2);
    R=gbest(:,1:d)*transpose(x);
    [value,idx]=max(R);
    index=find(Classes==y);
    zu(index) = zu(index)+ 1;
    if(index==idx )
        r(idx)= r(idx) + 1;
   
    end
    
    end
   i=i+1;

mul=1;

for i=1:k
    if(zu(i)~=0)
        if(r(i)==0)
            r(i)=r(i)+1;
        end
        vas=r(i)/zu(i);
    
        vas
        mul=mul*vas;
    end
    
   
end
acc=sum(r)/sum(zu)
%Gmean value
gmean=nthroot(mul,k)
final_mean=final_mean+gmean;
end
% Avg. Gmean value after 10 iteration
gm=final_mean/10;
%     x=1:gen;
%     plot(x,s(x))
%     drawnow
%     text(gen,s(gen),num2str((s(gen))))
%     axis([1 gen -2 1000]);
%     title('Plot of Global Fitness Vs No. of Iterations');
%     xlabel('No. of Iterations');
%     ylabel('Global Fitness');
%    } 
toc;   
%%
% PSO ends here-------------------------------------



%[Out] = PSO(Data, pop,gen,k,d,n);
% best fitness to best fit ;changed values between -1 and 1