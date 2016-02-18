I=imread('f.jpg');
Is=im2double(I);

temp=size(Is);
image_rows=temp(1,1);
image_columns=temp(1,2);
image_colors=temp(1,3);
w_low=double(0.7);
w_up=double(0.3);
ahh=0;
l=1;
for k=1:image_colors
    for i=1:image_rows
        for j=1:image_columns
            Im(l,k)=Is(i,j,k);
            l=l+1;
        end
    end
    l=1;
end

testing=0.001:0.001:0.01;
dimensions=size(testing);
limit=dimensions(1,2);
d_values=double(zeros(1,limit)+0.0);
for indexing=1:limit

temp=size(Im);
datapoints=temp(1,1);
attributes=temp(1,2);
C=4;
c_means=double(zeros(attributes,C)+0.0);

c_means(:,1)=[0.22;0.13;0.56];
c_means(:,2)=[0.12;0.33;0.86];
c_means(:,3)=[0.256;0.653;0.56];
c_means(:,4)=[0.211;0.82;0.26];


c_approx=-1+zeros(2,C,datapoints);


distance=0.0+zeros(datapoints,C);


U=zeros(datapoints,C);
for i=1:C
    c_approx(1,i,1)=i;
    %c_approx(2,i,1)=i;
    U(i,i)=1.0;
end



X_bar=zeros(1,attributes);
for i=1:datapoints
    X_bar=X_bar+Im(i,:);
end
X_bar=X_bar./datapoints;



sig_square=0.0;
for i=1:datapoints
    sig_square=sig_square+sum((Im(i,:)-X_bar).^2);
end
sig_square=sig_square/datapoints;
sig_square=2*sig_square;

for index=(C+1):datapoints
    index
    for j=1:C
        Im(index,:);
        c_means(:,j);
        distance(index,j)=sum(((Im(index,:)).^2-(transpose(c_means(:,j))).^2).^2);
        distance(index,j)=exp((-1.0*distance(index,j))/sig_square);
        distance(index,j)=2.0*(1.0-distance(index,j));
    end
    
    
    
    for j=1:C
        if distance(index,j)==0
            U(index,j)=1;
        else
            summation=0;
            for k=1:C
                summation=summation+(distance(index,j)/distance(index,k))^2;
            end
            U(index,j)=1/summation;
        end
        pi=1.0-U(index,j)-((1.0-U(index,j))/(1.0+2.0*U(index,j)));
        U(index,j)=U(index,j)+pi;
    end
    
    
    [A,pos]= sort(U(index,:),'descend');
    max1=pos(1);
    max2=pos(2);
    
    diff=U(index,max1)-U(index,max2);
    max1;
    max2;
    
    
    if (U(index,max1)-U(index,max2))<=testing(indexing)
        flag=1;
        while(c_approx(2,max1,flag)~=-1) 
            flag=flag+1;
        end
            
        c_approx(2,max1,flag)=index;
        
        count2=1;
        dist_sum=double(0.0);
        while(c_approx(2,max1,count2)~=-1 && count2<=datapoints)
            dist_sum=dist_sum+(U(c_approx(2,max1,count2),max1))^2;
            count2=count2+1;
        end
        gamma=dist_sum/double(count2-1);
      
        
        count1=1;
        count2=1;
        sum1=double(zeros(1,attributes)+0.0);
        sum2=double(zeros(1,attributes)+0.0);
        sum2a=double(0.0);
    
        while((c_approx(1,max1,count1)~=-1 && count1<=datapoints))  
            sum1=sum1+double(Im(c_approx(1,max1,count1),:));
            count1=count1+1;
        end
        sum1;
        while(c_approx(2,max1,count2)~=-1 && count2<=datapoints)   
            sum2=sum2+double(Im(c_approx(2,max1,count2),:).*double((1.0/(1.0+((U(c_approx(2,max1,count2),max1))^2)/gamma))^2));
            sum2a=sum2a+double((1.0/(1.0+((U(c_approx(2,max1,count2),max1))^2)/gamma))^2);
            count2=count2+1;
        end
        
        if ((count1-1)~=0 && (count2-1)~=0)
            c_means(:,max1)=((sum1.*w_low)./(double(count1-1))+((sum2.*w_up)./sum2a)); 
        end               
        if ((count1-1)~=0 && (count2-1)==0)
            c_means(:,max1)=(sum1./double(count1-1));
        end
        lower='no';
        c_means(:,max1);
        %for second cluster
        flag=1;
        while(c_approx(2,max2,flag)~=-1) 
            flag=flag+1;
        end
            
        c_approx(2,max2,flag)=index;
        
        count2=1;
        dist_sum=double(0.0);
        while(c_approx(2,max2,count2)~=-1 && count2<=datapoints)
            dist_sum=dist_sum+(U(c_approx(2,max2,count2),max2))^2;
            count2=count2+1;
        end
        gamma=dist_sum/double(count2-1);
      
        
        count1=1;
        count2=1;
        sum1=double(zeros(1,attributes)+0.0);
        sum2=double(zeros(1,attributes)+0.0);
        sum2a=double(0.0);
    
        while((c_approx(1,max2,count1)~=-1 && count1<=datapoints))  
            sum1=sum1+double(Im(c_approx(1,max2,count1),:));
            count1=count1+1;
        end
        sum1;
        while(c_approx(2,max2,count2)~=-1 && count2<=datapoints)   
            sum2=sum2+double(Im(c_approx(2,max2,count2),:).*double((1.0/(1.0+((U(c_approx(2,max2,count2),max2))^2)/gamma))^2));
            sum2a=sum2a+double((1.0/(1.0+((U(c_approx(2,max2,count2),max2))^2)/gamma))^2);
            count2=count2+1;
        end
       
        if ((count1-1)~=0 && (count2-1)~=0)
            c_means(:,max2)=((sum1.*w_low)./(double(count1-1))+((sum2.*w_up)./sum2a)); 
        end               
        if ((count1-1)~=0 && (count2-1)==0)
            c_means(:,max2)=(sum1./double(count1-1));
        end
        lower='no';
        c_means(:,max2);
        
    else
        flag=1;
        while(c_approx(1,max1,flag)~=-1) 
            flag=flag+1;
        end
        c_approx(1,max1,flag)=index;
        
        count2=1;
        dist_sum=double(0.0);
        while(c_approx(2,max1,count2)~=-1 && count2<=datapoints)
            dist_sum=dist_sum+(U(c_approx(2,max1,count2),max1))^2;
            count2=count2+1;
        end
        gamma=dist_sum/double(count2-1);
        
        
        count1=1;
        count2=1;
        sum1=double(zeros(1,attributes)+0.0);
        sum2=double(zeros(1,attributes)+0.0);
        sum2a=double(0.0);
    
        while((c_approx(1,max1,count1)~=-1 && count1<=datapoints))  
            sum1=sum1+double(Im(c_approx(1,max1,count1),:));
            count1=count1+1;
        end
        
        
        
        while(c_approx(2,max1,count2)~=-1 && count2<=datapoints)   
            sum2=sum2+double(Im(c_approx(2,max1,count2),:).*double((1.0/(1.0+((U(c_approx(2,max1,count2),max1))^2)/gamma))^2));
            sum2a=sum2a+double((1.0/(1.0+((U(c_approx(2,max1,count2),max1))^2)/gamma))^2);
            count2=count2+1;
        end
            
        
        if ((count1-1)~=0 && (count2-1)~=0)
            c_means(:,max1)=((sum1.*w_low)./(double(count1-1))+((sum2.*w_up)./sum2a)); 
        end               
        if ((count1-1)~=0 && (count2-1)==0)
            c_means(:,max1)=(sum1./double(count1-1));
        end
        lower='yes';
        c_means(:,max1);
    end
    
   
end

c_means;

C1=c_approx(1,1,:)~=-1;
C2=c_approx(1,2,:)~=-1;
C3=c_approx(1,3,:)~=-1;
C4=c_approx(1,4,:)~=-1;


C1=c_approx(1,1,C1);
C2=c_approx(1,2,C2);
C3=c_approx(1,3,C3);
C4=c_approx(1,4,C4);

C1(:)
C2(:)
C3(:)
C4(:)





size_approx=zeros(C,2);
q=zeros(1,datapoints);
r=zeros(1,datapoints);
for j=1:C
    for i=1:datapoints
        q(i)=c_approx(1,j,i);
        r(i)=c_approx(2,j,i);
    end
    a=q~=-1;
    answer=q(a);
    a=size(answer);
    size_approx(j,1)=a(2);
    a=r~=-1;
    answer=r(a);
    a=size(answer);
    size_approx(j,2)=a(2);
end
size_approx
flag=1;
while(c_approx(1,1,flag)~=-1 && flag<=datapoints) 
    Im(c_approx(1,1,flag),:)=[0.1111,0.0,0.0];
    flag=flag+1;
end

flag=1;
while(c_approx(1,2,flag)~=-1 && flag<=datapoints) 
    Im(c_approx(1,2,flag),:)=[0.0,0.3333,0.0];
    flag=flag+1;
end

flag=1;
while(c_approx(1,3,flag)~=-1 && flag<=datapoints) 
    Im(c_approx(1,3,flag),:)=[0.0,0.0,0.6666];
    flag=flag+1;
end

flag=1;
while(c_approx(1,4,flag)~=-1 && flag<=datapoints) 
    Im(c_approx(1,4,flag),:)=[0.9999,0.0,0.9999];
    flag=flag+1;
end

l=1;
for k=1:image_colors
    for i=1:image_rows
        for j=1:image_columns
            Is(i,j,k)=Im(l,k);
            l=l+1;
        end
    end
    l=1;
end

Iback=im2uint8(Is);
imshow(Iback);




%ACCURACY MEASURES
S=double(zeros(1,C)+0.0);
for i=1:C
    count2=1;
    dist_sum=double(0.0);
    while(c_approx(2,i,count2)~=-1 && count2<=datapoints)
        dist_sum=dist_sum+(U(c_approx(2,i,count2),i))^2;
        count2=count2+1;
    end
    gamma=dist_sum/double(count2-1);
    count1=1;
    count2=1;
    sum1=double(0.0);
    sum2=double(0.0);
    sum2a=double(0.0);
    
    while((c_approx(1,i,count1)~=-1 && count1<=datapoints))  
        sum1=sum1+double((distance(c_approx(1,i,count1),i))^2);
        count1=count1+1;
    end
    
    while(c_approx(2,i,count2)~=-1 && count2<=datapoints)   
        sum2=sum2+double((distance(c_approx(2,i,count2),i))^2)*double((1.0/(1.0+((U(c_approx(2,i,count2),i))^2)/gamma))^2);
        sum2a=sum2a+double((1.0/(1.0+((U(c_approx(2,i,count2),i))^2)/gamma))^2);
        count2=count2+1;
    end
    
    if ((count1-1)~=0 && (count2-1)~=0)
        S(i)=((sum1*w_low)/(double(count1-1))+((sum2*w_up)/sum2a)); 
    end               
    if ((count1-1)~=0 && (count2-1)==0)
        S(i)=(sum1/double(count1-1));
    end
end
%till here its correct
c_distance=double(zeros(C)+0.0);

for i=1:C
    for j=1:C
        c_distance(i,j)=sqrt(sum((c_means(:,i)-c_means(:,j)).^2));
    end
end

t1=double(zeros(1,C-1)+0.0);
t2=double(zeros(1,C-1)+0.0);
t3=double(zeros(1,C)+0.0);

db_index=double(0.0);
d_index=double(0.0);
for i=1:C
    k=1;
    for j=1:C
        if i~=j
            t1(k)=((S(i)+S(j))/c_distance(i,j));
            k=k+1;
        end
    end
    db_index=db_index+max(t1);
end
db_index=db_index/double(C)
max_S=double(max(S));
f=1;

for i=1:C
    k=1;
    for j=1:C
        if i~=j
            t2(k)=c_distance(i,j)/max_S;
            k=k+1;
        end
    end
    t3(f)=min(t2);
    f=f+1;
end
d_index=min(t3)
d_values(indexing)=d_index;

clearvars -except Im testing dimensions limit d_values indexing temp image_rows image_colors image_columns w_low w_up ahh

end

[d_values,posit]= sort(d_values,'descend');
max1=posit(1);
testing(max1)

























    
    




