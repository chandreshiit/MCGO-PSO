% dataset={'connect-4','covtype','dna','glass','iris','letter','shuttle','usps','vehicle','vowel2csv','wine2csv'};
dataset={'connect-4','dna','glass','iris','letter','shuttle','usps','vehicle','vowel2csv','wine2csv','combined','acoustic','pendigits','poker','satimage','segment','seismic','svmguide4','sensorless'};
nn={67556,1999,213,149,14999,43499,7290,845,527,177,78822,78822,7493,25009,4434,1210,2309,299,58508};
dd={126,180,9,4,16,9,256,18,10,13,100,50,16,10,36,19,50,10,48};
kk={3,3,6,3,26,7,10,4,11,3,3,3,10,10,6,7,3,6,11};
for w=1:length(dataset)

    str=['dataset/',char(dataset(w)),'.csv'];
    Data=csvread(str,1,0);

    k=cell2mat(kk(w));
    % n1=ceil(.6*108000);
    % ntotal=108000;
    n1=cell2mat(nn(w));
    d=cell2mat(dd(w));
    c1=1.2;
    c2=2;
    Classes=unique(Data(:,1));
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

     [fitnessVal] = Fitness2( Data, psomatrix_x(:,:,j),k,d,n1); % Fitness
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

    [fitnessVal] = Fitness2(Data, psomatrix_x(:,:,i),k,d,n1 ); % Fitness
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
    q=zeros(1,k);
    zu=zeros(1,k);
    %i=n1;
    %jp=ntotal-n1+1;
    jp=n1;
    for i=1:jp

        x=Data(i,2:d+1);
        y=Data(i,1);
        %p(i)=Data(i,1);
        R=gbest(:,1:d)*transpose(x);
        [value,idx]=max(R);
        index=find(Classes==y);
        zu(index) = zu(index)+ 1;
        if(index==idx)
            r(idx)= r(idx) + 1;
        else
            q(idx)= q(idx) + 1;
        end
        if(max(x)~=0)
        [value1,idx1] = max(R(R<value));
        idx1=find(R==value1);
        [value2,idx2] = max(R(R<value1));
        idx2=find(R==value2);
        if(index==idx1)
        r(idx1)= r(idx1) + 1;
        q(idx)= q(idx) - 1; 
    else
         q(idx)= q(idx) - 1; 
      q(idx1)= q(idx1) + 1;
    end
 end
   
    end 
    mul=1;
    f_sum=0;
    for i=1:k
        if(zu(i)~=0)
         if(r(i)==0)
        r(i)=r(i)+1;
        end
            vas=r(i)/zu(i);

            vas
            mul=mul*vas;
        else
            vas=0;
        end

        if((r(i)+q(i))~=0)
            prec=r(i)/(r(i)+q(i));
        else
            prec=0;
        end
        if((prec+vas)~=0)
          f_scoe=2*(prec*vas)/(prec+vas);
        else
            f_scoe=0;
        end
    f_sum=f_sum+f_scoe;
    end
    %Accuracy value
    acc=sum(r)/sum(zu)
    %Gmean value
    gmean=nthroot(mul,k)
    %F-score value
     f_score=f_sum/k

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
end
%%
% PSO ends here-------------------------------------



%[Out] = PSO(Data, pop,gen,k,d,n);
% best fitness to best fit ;changed values between -1 and 1