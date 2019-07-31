function [Value] = Fitness2( Data, psomatrix_x,k,d,n1 )
%Data=csvread('Data.csv',1,0);
Classes=unique(Data(:,1));
%k=7;
% n=8774;
% d=524;
Clas_count=zeros(1,k);
CYI=zeros(1,k);
W=psomatrix_x(:,:);
z=zeros(1,n1);
length(W);
length(z);
for i=1:k
  for j=1:n1
      if (Data(j,1)==Classes(i))
          Clas_count(i)= Clas_count(i)+1;
      end
  end
  CYI(i)=10000/Clas_count(i);
end
summ=0;
t=(1/2)*norm(W)*transpose(norm(W));
for i=1:n1
    x=Data(i,2:d+1);
    y=Data(i,1);
    index=find(Classes==y);
    
    a=CYI(index);
    mas=-999;
    for j=1:k
        if(j~=index)
                mul1=x*transpose(W(j,:));
                mul2=x*transpose(W(index,:));
                num2=mul1-mul2;
%                 num3=mul2-mul1;
%                 if(num3<(1-z(index)))
%                     add=1-z(index)-num3;
% 
%                     z(index)=z(index)+add;
% 
%                 end
                if(num2>mas)
                    mas=num2;
                end
            
        end
  
    end
    if(mas<=0)
        s=0;
    else
        s=mas+1;
    end
    z(i)=a*s*s;
    summ=summ+z(i);
    
end
value=summ+t;

Value =min(value);
Value=min(Value);
end