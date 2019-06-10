clear all
clc

Xmin=-100;
Xmax=100;
runs=5;
rand('state',sum(100*clock));
x__axis=1.2*rand();
y__axis=1.2*rand();


fhd=str2func('cec13_func');
foptimal=[-1400:100:-100 100:100:1400];                                                                           

tipoF      = 2;
opposition = 2;
tipoCR     = 2;

fname_FFOA = strcat('FFOA','_opp_',num2str(opposition),'_tf_',num2str(tipoF),'_CR_',num2str(tipoCR),'_');

 for D=10:20:50

  clear xbest fbest
  
  for i=1:28 
      
    switch D
      case 10, pop_size=50; iter_max=2000; 
      case 30, pop_size=100; iter_max=3000; 
      case 50, pop_size=150; iter_max=3333;
    end
    
    MaxFES=pop_size*iter_max;
    Record_MaxFES=[0.01 0.1:0.1:1]*MaxFES;    
    Inf_matrix=zeros(11,5);
    
    clear error_best xbest fbest error_best
    
    for j=1:runs
        [x_best,f_bestval,FES,Inf_matrix] = FFOA_sub(fhd,D,pop_size,iter_max,Xmin,Xmax,Record_MaxFES,foptimal(i),j,Inf_matrix,tipoF,tipoCR,opposition,x__axis,y__axis,i);
        xbest(j,:)=x_best;
        fbest(i,j)=f_bestval;
        error_best(i,j) =f_bestval -foptimal(i);
        fbest_mean = mean(fbest(i,:)-foptimal(i));
        
        fprintf('Run: %d/%d\t\tFunc: %d D: %d\t\tFbest: %f\tFopt: %f\t\tError: %e\t\tFES: %d/%d\tMeanF: %f\n', ...
                j,runs,i,D,fbest(i,j),foptimal(i),fbest(i,j)-foptimal(i),FES,MaxFES,fbest_mean);
    end
    statistics_f(i,:) = [min(error_best(i,:)) max(error_best(i,:)) median(error_best(i,:)) mean(error_best(i,:)) std(error_best(i,:))];
    
    fname = strcat([fname_FFOA,'Preliminar_FFOA_',num2str(i),'_',num2str(D),'.txt']);
    periodic_f = fopen(fname, 'at');       
    for cont=1:11
      for c_run = 1:runs
          
   
        fprintf(periodic_f,'%e ',Inf_matrix(cont,c_run));
      end
      fprintf(periodic_f,'\n');
    end
    fclose(periodic_f);

    fname = strcat([fname_FFOA,'Preliminar__stat_5_runs_f',num2str(i),'_',num2str(D),'.txt']);    
    sta_f = fopen(fname, 'at');  
    for c_run = 1:runs
      fprintf(sta_f,'%e\n',Inf_matrix(11,c_run));
    end
    fprintf(sta_f,'%e\n',min(Inf_matrix(11,:)));
    fprintf(sta_f,'%e\n',max(Inf_matrix(11,:)));
    fprintf(sta_f,'%e\n',median(Inf_matrix(11,:)));
    fprintf(sta_f,'%e\n',mean(Inf_matrix(11,:)));
    fprintf(sta_f,'%e\n',std(Inf_matrix(11,:)));
    fprintf(sta_f,'\n');
    fclose(sta_f);
    
    fname = strcat([fname_FFOA,'Preliminar__Table_stat_5_runs_f',num2str(i),'_',num2str(D),'.txt']);    
    sta_f2 = fopen(fname, 'at');  
    fprintf(sta_f2,'%d\t',i);    
    fprintf(sta_f2,'%e\t',min(Inf_matrix(11,:)));
    fprintf(sta_f2,'%e\t',max(Inf_matrix(11,:)));
    fprintf(sta_f2,'%e\t',median(Inf_matrix(11,:)));
    fprintf(sta_f2,'%e\t',mean(Inf_matrix(11,:)));
    fprintf(sta_f2,'%e\t',std(Inf_matrix(11,:)));
    fprintf(sta_f2,'\n');
    fclose(sta_f2);
    
    fname2 = strcat([fname_FFOA,'results_f',num2str(i),'_',num2str(D),'D.txt']);
    save(fname2);
  end
  fname3 = strcat([fname_FFOA,'results_ALL_',num2str(D),'D.txt']);
  save(fname3);
 end
 fname4 = strcat([fname_FFOA,'results_ALL_f.txt']);
 save(fname4);

