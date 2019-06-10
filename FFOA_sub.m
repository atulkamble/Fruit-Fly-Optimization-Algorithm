function [bestmem,bestval,nfeval,Inf_matrix] = FFOA_sub (fhd,D,PS,itermax,XVmin,XVmax,Record_MaxFES,fopt,j2,Inf_matrix,tipoF,tipoCR,opposition,x__axis,y__axis,varargin)

rand('state',sum(100*clock));

without_improve =  0;  

if length(XVmin)==1
  XVmin=repmat(XVmin,1,D);
  XVmax=repmat(XVmax,1,D);
end

XVmin=repmat(XVmin,PS,1);
XVmax=repmat(XVmax,PS,1);
pop=XVmin+(XVmax-XVmin).*rand(PS,D);
val=feval(fhd,pop',varargin{:});
nfeval=PS;
[bestval,indice] = min(val);
bestval_ant = bestval;

bestmem=pop(indice,:);

rot=(0:1:PS-1);               
rotd=(0:1:D-1);               

iter=1;
popold=pop;
CR_ant=rand;
VTR=1e-8;
cont= 1;

while ((iter < itermax) && ((bestval-fopt) >= VTR))
  previous_nfeval= nfeval;
  
  variance_pop(iter,:)=var(pop);    
  max_variance=max(variance_pop);
  diversity_norm=(D/20)*(variance_pop(iter,:)./max_variance);
  L= max(  min(0.6,1-diversity_norm), 0.3);
  
  % F
  Fdec= 0.7*(iter/itermax) + 0.2;   
  for i=1:PS
    for j=1:D  
      if rand>0.2  
        F(i,j)=L(1,j);
      else
        F(i,j)=Fdec;                   
      end
    end
  end  
       
  % CR
  if without_improve > 90  % without_improve >>10 para F14
    CR = 0.1*max(1-diversity_norm);            
  else
    CR = max(1-diversity_norm);
  end 
      
  
  popold=pop;                   
  ind=randperm(2);              

  a1=randperm(PS);             
  rt=rem(rot+ind(1),PS);        
  a2=a1(rt+1);                 
  rt=rem(rot+ind(2),PS);
  a3=a2(rt+1);                
  
  pm1=popold(a1,:);             
  pm2=popold(a2,:);             
  pm3=popold(a3,:);             
  
  
  bm=repmat(bestmem,PS,1);
  
  if length(CR) > 1
      mui=rand(PS,D)< repmat(CR,PS,1);
  else
    mui=rand(PS,D)<CR;          
  end
  
  mui=sort(mui');	          
  for i=1:PS
    n=floor(rand*D);
    rtd=rem(rotd+n,D);
    mui(:,i)=mui(rtd+1,i); 
  end
  mui=mui';			  
  mpo=mui<0.5;        

  
  if rand > 0.5  
    
    ui = popold + F.*(bm-popold) + F.*(pm1 - pm2);         
    
  else    
    
    ui = pm3 + F.*(pm1 - pm2);       
  end
  ui = popold.*mpo + ui.*mui;     
  
  % opposition
  if without_improve > 90  % without_improve >>10 para F14
    prob_opposition= (D/60)*diversity_norm;
  else    
    prob_opposition= (D/33.3333334)*diversity_norm;
  end  
        
  for i=1:PS
    for j=1:D        
       if rand < prob_opposition(j)  
         miu = min(ui(:,j));
         mau = max(ui(:,j));
         ui(i,j) =  miu + mau - ui(i,j);
         Op_popVar(i,j) = miu +  mau - ui(i,j);
         M(i,j)         = (miu + mau) / 2;
         if ui(i,j) < M(i,j)
           ui(i,j) = M(i,j) + (Op_popVar(i,j) - M(i,j))*rand;
         else
           ui(i,j) = Op_popVar(i,j) + (M(i,j) - Op_popVar(i,j))*rand;
         end                               
      end 
    end
  end
         
  % bounds   
  ui=(ui>XVmax).*XVmax+(ui<=XVmax).*ui;     
  ui=(ui<XVmin).*XVmin+(ui>=XVmin).*ui;
  
  tempval=feval(fhd,ui',varargin{:});
  nfeval=nfeval+PS;
  
  [best_tempval,indice]=min(tempval);
  if best_tempval < bestval,     
    bestval = best_tempval;      
    bestmem = ui(indice,:);
  end
  
  for i=1:PS,    
      if tempval(i) <= val(i),         
          pop(i,:) = ui(i,:);         
          val(i)   = tempval(i);      
      end
  end
  
  for i=1:PS,
      x(i)=x__axis+1.2*rand()-0.5;
      y(i)=y__axis+1.2*rand()-0.5;
      
      distance(i)=(x(i)^2+y(i)^2)^0.5;
      s(i)=1/distance(i);
      
      smell(i)=3*s(i)^2+5;        
  end
      
     [bestsmell bestindex]=max(smell);
     
     x__axis=x(bestindex);
     y__axis=y(bestindex);
     smellbest=bestsmell;
     
     for g=1:itermax
     
         for i=1:PS
              x(i)=x__axis+1.2*rand()-0.5;
              y(i)=y__axis+1.2*rand()-0.5;
              
             distance(i)=(x(i)^2+y(i)^2)^0.5;
             s(i)=1/distance(i);
             
              smell(i)=3*s(i)^2+5; 
         end
         [bestsmell bestindex]=max(smell);
     
         if bestsmell>smellbest;
         x__axis=x(bestindex);
         y__axis=y(bestindex);
         smellbest=bestsmell;
         end
         
         %%added newly
           yy(g)=smellbest; 
           Xbest(g)=x__axis;
           Ybest(g)=y__axis;

           %%added
         
           %%recently added to graph 
           
 %{   
    figure(1)
    plot(yy)
    grid  on;
    title('Optimization Process','fontsize',14)
    xlabel('Iteration Number','fontsize',12);ylabel('Smell','fontsize',14);

    figure(2)
    plot(Xbest,Ybest,'b.');
    grid on;
    title('Fruit fly flying route ','fontsize',14)
    xlabel('X-axis','fontsize',12);ylabel('Y-axis','fontsize',12);
    %pause(0.5)
 %}          
           
           
     end
     
  if (bestval-fopt) < VTR
    bestval = fopt;
  end
  
  if (nfeval >= Record_MaxFES(cont)) && (previous_nfeval < Record_MaxFES(cont))
      
    Inf_matrix(cont,j2) = bestval-fopt;
    cont =  cont + 1;  
  end    
  
  if (nfeval >= Record_MaxFES(end)) && (previous_nfeval <= Record_MaxFES(end))
    Inf_matrix(cont,j2) = bestval-fopt;      
    cont =  cont + 1;  
  end    
  
  if bestval == bestval_ant
    without_improve =  without_improve + 1;
  else  
    without_improve =  0;  
  end
  bestval_ant = bestval;
  iter = iter + 1;
end 
