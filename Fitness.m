function [Value] = Fitness( Data, psomatrix_x,k,d,n )
%Data=csvread('Data.csv',1,0);
Classes=unique(Data(:,d+2));
%k=7;
% n=8774;
% d=524;
Clas_count=zeros(1,k);
CYI=zeros(1,k);
W=psomatrix_x(:,1:d);
z=psomatrix_x(:,d);
length(W);
length(z);
for i=1:k
  for j=1:n
      if (Data(j,d+2)==Classes(i))
          Clas_count(i)= Clas_count(i)+1;
      end
  end
  CYI(i)=1/Clas_count(i);
end

summ=0;
t=(1/2)*norm(W*transpose(W));
for i=1:n
    x=Data(i,2:d+1);
    y=Data(i,d+2);
    index=find(Classes==y);
    
    a=CYI(index);
    mas=0;
    for j=1:k
        if(j~=index)
                mul1=x*transpose(W(j,:));
                mul2=x*transpose(W(index,:));
                num2=mul1-mul2;
                num3=mul2-mul1;
                if(num3<(1-z(index)))
                    add=1-z(index)-num3;

                    z(index)=z(index)+add;

                end
                if(num2>mas)
                    mas=num2;
                end
            
        end
  
    end
    if(mas+1<0)
        s=0;
    else
        s=mas+1;
    end
    g=a*s;
    summ=summ+g;
    
end
value=summ+t;

Value =min(value);
Value=min(Value);
end