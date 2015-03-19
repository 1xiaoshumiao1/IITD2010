function [mat1,c,im3,result]=identify1(im,mat)
%im=imread(fname);
%im=rgb2gray(im);
%im=imresize(im,1/4);

%mat=roipoly(im);


x=size(mat,1);
y=size(mat,2);
mat1=zeros(x,y);
result=zeros(x,y);
result=mat;
% iteration will start here
% Calculating boundry (initial)
% mat contains tge region to be inpainted
for(i=1:x)
    for(j=1:y)
            if(mat(i,j)==1)
            if(mat(i+1,j)==0||mat(i-1,j)==0||mat(i,j+1)==0||mat(i,j-1)==0)
                mat1(i,j)=1;
            end;
        end;
    end;
end;
mat2=zeros(x,y);
mat2=mat1;
% mat1 contains the boundry pixels
c=im;
im3=c;
confidence=ones(x,y);
confidence=confidence-mat;
c=(im).*uint8(confidence);
im3=c;

% confidence contains the confidence of each pixel
% c contains the updated image
for(i=1:x)
    for(j=1:y)
        if(mat1(i,j)==1)
           confidence(i,j)=confidence(i,j-1)+confidence(i,j+1)+confidence(i,j)+confidence(i-1,j-1)+confidence(i-1,j+1)+confidence(i-1,j)+confidence(i+1,j-1)+confidence(i+1,j+1)+confidence(i+1,j);
           confidence(i,j)=confidence(i,j)/9;
        end;
    end;
end;
%imshow(confidence);


           h=-1;
  while(h<100)
    
      max=0;
normal=zeros(x,y);

dotp=zeros(x,y);



iy=edge(c,'sobel','v');
ix=edge(c,'sobel','h');
normal=zeros(x,y);
for(i=1:x)
    for(j=1:y)
        if(mat1(i,j)==1)
        normal(i,j)=-iy(i,j)+(ix(i,j)*1i);
    end;
end;
end;
a=[-1 -1 -1;-1 8 -1; -1 -1 -1];
l=conv2(double(im3),double(a));
dl=zeros(x,y);
for(i=1:x)
    for(j=1:y)
        if(mat1(i,j)==1)
            k=l(i+1,j)-l(i-1,j);
            m=l(i,j+1)-l(i,j-1);
            dl(i,j)=k+(m*1i);
        end;
    end;
end;


beta=zeros(x,y);
for(i=1:x)
    for(j=1:y)
        if(mat1(i,j)==1)
            beta(i,j)=(real(dl(i,j))*real(normal(i,j)))-(imag(dl(i,j))*imag(normal(i,j)));
            norm=real(normal(i,j))*real(normal(i,j));
            norm=norm+(imag(normal(i,j))*imag(normal(i,j)));
            norm=sqrt(norm);
            if(norm~=0)
           beta(i,j)=beta(i,j)/norm;
            end;
        end;
    end;
end;
ix=edge(c,'sobel','h');
iy=edge(c,'sobel','v');
di=zeros(x,y);
for(i=2:x-1)
    for(j=2:y-1)
        if(mat1(i,j)==1)
            cx=ix(i,j)-ix(i,j-1);
            dx=ix(i,j+1)-ix(i,j);
            di(i,j)=cx*cx+dx*dx;
            cx=iy(i,j)-iy(i-1,j);
            dx=iy(i+1,j)-iy(i,j);
            di(i,j)=di(i,j)+cx*cx+dx*dx;
            di(i,j)=sqrt(di(i,j));
        end;
    end;
end;
di=di.*beta;
priority=zeros(x,y);
%priority=confidence.*dotp;
priority=confidence.*di;
%priority=confidence.*G;
maxp=-999;
u=2;
v=2;
for(i=2:x-1)
    for(j=2:y-1)
        if(mat1(i,j)==1)
            if(priority(i,j)>maxp)
                u=i;
                v=j;
                maxp=priority(i,j);
            end;
        end;
    end;
end;


% finding a similar source fragment for a target fragment
max1=zeros(x,y);
xcord=zeros(x,y);
ycord=zeros(x,y);
im3=double(im3);
% converting image type to uint32
% max1 contains value for a target-source fragment similarity
                          i=double(u);
                          j=double(v);
                          e=u;
                          f=v;

            max=double(99999999999999999999999999999);
            u=0;
            v=0;
            for(p=4:x-3)
                for(q=4:y-3)
                     a=0;
                   for(k=p-3:p+3)
                       for(l=q-3:q+3)
                           if(mat(k,l)==1)
                              a=1;
                           end;
                       end;
                   end;
                   value=0;
                   if(a==0)
                       if(e>0&&f>0)
                       for(m=-3:1:3)
                           for(n=-3:1:3)
                               if(mat(e+m,f+n)==0)
                                   k=double(im3(p+m,q+n));
                                   
                                   
                                   l=double(im3(i+m,j+n));
                                   
                                   value=double(value+double(abs(double(k)*double(k)-(double(l)*double(l)))));
                                 
                               end;
                           end;
                       end;
                       end;
                      
                        
                      if(value<max)
                          u=p;
                          v=q;
                          max=double(value);
                      end;
                   end;
                end;
            end;
            
            
           
            
          
            
            max1(i,j)=max;
            xcord(i,j)=u;                  % xcord contains x coordinate of the source patch
                                            
            ycord(i,j)=v;                  % ycord contains y coordinate of source patch
          
            
           for(i=-3:3)
               for(j=-3:3)
                   if(mat(e+i,f+j)==1)
                       im3(e+i,f+j)=im3(u+i,v+j);
                       mat(e+i,f+j)=0;
                       confidence(e+i,f+j)=confidence(e,f);
                   end;
               end;
           end;
           
           mat1=zeros(x,y);
           for(i=1:x)
              for(j=1:y)
                if(mat(i,j)==1)
                   if(mat(i+1,j)==0||mat(i-1,j)==0||mat(i,j+1)==0||mat(i,j-1)==0)
                     mat1(i,j)=1;
            end;
        end;
    end;
end;

           
           
  for(i=1:x)
    for(j=1:y)
        if(mat1(i,j)==1)
           confidence(i,j)=confidence(i,j-1)+confidence(i,j+1)+confidence(i,j)+confidence(i-1,j-1)+confidence(i-1,j+1)+confidence(i-1,j)+confidence(i+1,j-1)+confidence(i+1,j+1)+confidence(i+1,j);
           confidence(i,j)=confidence(i,j)/9;
        end;
    end;
end;
h=1000;
for(i=1:x)
   for(j=1:y)
        if(mat(i,j)==1)
            h=-1;
        end;
    
  end;
end;
  end;
  im3=uint8(im3);
  im3=double(im3);
  avg=0;
%  for(i=1:x)
 %     for(j=1:y)
  %        if(mat2(i,j)==1)
   %           avg=0;
    %          for(p=-1:1:1)
     %            for(q=-1:1:1)
      %              avg=double(double(avg)+double(im3(i+p,j+q)));
      %            end;
       %       end;
       %       avg=double(avg/9);
  %            im3(i,j)=avg;
   %       end;
    %  end;
  %end;
 
      im3=uint8(im3);
      %a=[-1 -1 -1;0 0 0;1 1 1];
      %im4=conv2(double(im3),double(a));
      %im3=double(im3)+double(a);
      %im3=im3+im4;
result=uint8(result).*im3;