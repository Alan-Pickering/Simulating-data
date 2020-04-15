%reg_by_mvnrnd1.m
%
%v2; 14.4.20
%
%fully functional version
%@@@denotes areas for code improvement

clear variables;
clc;

myseed=2000;
rng(myseed);

disp('Matlab pgm: reg_by_mvnrnd1.m; version 2 (13.4.2020)');
disp('***************************************************');
disp('This programme let''s you simulate a large number of regression datasets with known properties.');
%@@@add default info and selection here
disp('Please select the properties of the datasets by typing the required info + <Enter> in each case.');
disp(' ');

maxvars=8; %maximum number of variables allowed; which could be changed at later date

%now select the options you want
nsims=0; %the numbers of simulations
while nsims<1 || nsims>1000
    nsims=input('How many simulations do you want to run? (1-1000) ');
    if isempty(nsims)
        nsims=1000;
    end
end

disp([ 'In each of these ' num2str(nsims) ' simulations:- ' ])
nvars=0; %number of variables per simulation
while nvars <2 || nvars>maxvars
    msg=['How many variables do you want to simulate? (2-' num2str(maxvars) ') '];
    nvars=input(msg);
    %specify the default if hit <return> only
    if isempty(nvars)
        nvars=3;
    end
end
dlist=1:nvars; %creates a numbered list of code for data variables (used for filtering later on)

n=0; %the number of cases per simulation
while n<10 || n>10000
    n=input('How many cases do you want to simulate? (10-10000) ');
    if isempty(n)
        n=200;
    end
end

myalpha=0;  %the type 1 error rate (significance criterion)
while myalpha<=0 || myalpha >=1
    myalpha=input('What type 1 error rate (aka alpha or significance level) do you want? (0<alpha<1) ');
    if isempty(myalpha)
        myalpha=0.05;
    end
end

%next is the variable information
vnames=cell(1,nvars);
my_var=zeros(1,nvars);
my_mu= -10001.*ones(1,nvars);
disp(' ');
disp('Now give details of your to-be-simulated variables:- ');
for v=1:nvars
    msg=['Type a short name for variable ' num2str(v) ' '];
    myname=input(msg,'s');
    if isempty(myname)
        myname=['v' num2str(v)];
    end
    vnames{v}=myname;
    while my_mu(v)<-10000 || my_mu(v)>10000
        msg=['Type value for mean of ' char(vnames(v)) ' (-10000< mean <=10000) '];
        muval=input(msg);
        if isempty(muval)
            muval=0;
        end
        my_mu(v)=muval;
    end
    while my_var(v)<=0 || my_var(v)>10000
        msg=['Type value for variance of ' char(vnames(v)) ' (0< var <=10000) '];
        varval=input(msg);
        if isempty(varval)
            varval=1;
        end
        my_var(v)=varval;
    end
end

disp(' ');
for j=1:nvars
    disp(['Variable #' num2str(j) ' is ' char(vnames(j))]);
end

%@@@select variable by name here maybe
whichisdv=0;
while whichisdv<1 || whichisdv>maxvars
    msg=['Type number of variable which is DV (1-' num2str(maxvars) ') '];
    whichisdv=input(msg);
    if isempty(whichisdv)
        whichisdv=1;
    end
end

%default values are not in operation from this point
%specify the desired correlations and create the var-covar (vc) matrix
init_c=-1.1*ones(nvars); %create initial correlation matrix with illegal values
set_corr=init_c;
set_vc=diag(my_var); %put the variances along the main diagonal of var-covar matrix
for i=1:nvars-1
    for j=i+1: nvars
        msg=['Type value for correlation, r, between ' char(vnames(i)) ' and ' char(vnames(j)) ' (-1 =< r <=1) '];
        while set_corr(i,j)<-1 || set_corr(i,j)>1
            set_corr(i,j)=input(msg);
        end
        %compute off diag elements of variance-covariance matrix from corr matrix and variances
        set_vc(i,j)=set_corr(i,j).*sqrt(my_var(i)).*sqrt(my_var(j));
    end
end

%complete matrices using symmetry properties
set_corr=set_corr-tril(init_c);
set_corr=set_corr+set_corr'+eye(nvars);
set_vc=set_vc+set_vc' -diag(diag(set_vc));

%@@@select variables to act as predictors
%all that are not DV by default
%use var names to select

disp(' ');
disp('Simulating...');
%set up arrays to capture simulation outputs
record_corr=zeros(nvars,nvars,nsims);
record_data=zeros(n, nvars,nsims);
Fpvals=zeros(nsims,1);
Fvals=zeros(nsims,1);
npreds=nvars-1; %-1 because 1 of variables is DV
tvals=zeros(nsims,npreds+1); %+1 because 1 intercept plus predictors
tpvals=zeros(nsims,npreds+1);
bvals=zeros(nsims,npreds+1);
R2vals=zeros(nsims,1);
%controls freq of disp messages during loop
printit=100; %every 100 simulations
for s=1:nsims
    
    if mod(s,printit)==0
        disp(['Simulation no. ' num2str(s) ' completed']);
    end
 
    simdata=mvnrnd(my_mu,set_vc,n); %all data generated using mvnrnd (multivariate normal RNG)
    
    %@@@here we could add a way of creating binary variables here
    %based on a standard normal varianle with correls [r+] and then median split
    %values in [r+] would need to be increasd above intended point biserial r 
    %to allow for loss via median split
    
    %make record of data values
    %first data
    record_data(:,:,s)=simdata;
    %then correlation between vars
    record_corr(:,:,s)=corr(simdata);

    %compute OLS regression info
    %@@@restrict predictors according to selection above
    mystats = regstats(record_data(:,dlist==whichisdv,s), record_data(:,dlist~=whichisdv,s), 'linear', {'tstat', 'fstat', 'rsquare'}); %OLS
    %record reg results for each sim
    bvals(s,:)=mystats.tstat.beta;
    tvals(s,:)=mystats.tstat.t;
    tpvals(s,:)=mystats.tstat.pval;
    Fvals(s)=mystats.fstat.f;
    Fpvals(s)=mystats.fstat.pval;
    R2vals(s)=mystats.rsquare;
    %@@@compute betas from b and store
    
end

disp(' ');
disp('Overall results, across all simulations')
disp('This is the mean simulated correlation matrix');
%first the correlation matrics
mcorr=mean(record_corr,3); %average over 3rd dimension = over simulations
%tablenames=cell(1,nvars+1);
%tablenames{1}='Variable';
%tablenames(2:nvars+1)=vnames(1:nvars);
%t=table(char(vnames'),mcorr(:,1),mcorr(:,2),mcorr(:,3),'VariableNames',tablenames);
%t.Properties.Description = 'Mean simulated correlation matrix';
%disp(t);
disp(mcorr);

disp('This is the intended correlation matrix');
disp(set_corr);

%now some regression summary data
%@@@compute the expected value of R-sq and disp it here
%also pr and sr values too perhaps
%give mean R-sq next
disp([ 'The mean R2 value= ' num2str(mean(R2vals)) ]);
%compute power of regression effects, across simulations
%count up the p values significant at level myalpha, for various reg tests stats
%first for F-test of R2
numsig=sum(Fpvals < myalpha); 
disp(['The Monte Carlo estimated power of the F-test overall reg model = ' num2str(numsig/nsims)]);

ctr=0;
termlab=cell(1,npreds); %create cell array to hold labels for the predictors
for j=1:nvars
    if dlist(j)~=whichisdv
        ctr=ctr+1;
        termlab{ctr}=vnames{j};
        numsig=sum(tpvals(:,ctr+1) < myalpha); %+1 to start after the intercept
        disp([ 'The mean b value for predictor ' termlab{ctr} ' = ' num2str(mean(bvals(:,ctr+1))) ]);
        disp(['The Monte Carlo estimated power of the t-test for predictor ' termlab{ctr} ' = ' num2str(numsig/nsims)]);
    end
end

%some plots of the whole simulation
figure;
histogram(R2vals);
ylabel('Frequency');
xlabel('R-squared values');
msg=[ 'Histogram for ' num2str(nsims) ' simulations' ];
title(msg);


%now add visualisation of an individual simulation
disp(' ');
disp('Now, for the results of particular simulation(s) ...');
checkout=0;
while checkout==0

    %select and give details of a specific simulation
    whichsim=-1;
    disp(' ');
    while (whichsim<0 || whichsim>nsims)
        whichsim=input([ 'Type number of simulation (1-' num2str(nsims) ') to see its results, followed by <Enter>; or 0 to move on ' ]);
        if isempty(whichsim)
            whichsim=0;
        end
    end
    if whichsim>0
        %give some data from the specific 
        %first correlation info
        disp(' ');
        disp([ 'The correlation matrix for simulation ' num2str(whichsim) ' was:-' ]);
        disp(record_corr(:,:,whichsim));
        %next regression info
        disp([ 'Simulation no. ' num2str(whichsim) ]);
        disp([ 'The R2 value for this simulation= ' num2str(R2vals(whichsim)) ]);
        disp([ 'The F value (and associated p val) for this simulation= ' num2str(Fvals(whichsim)) ' , ' num2str(Fpvals(whichsim))]);
        for j=1:npreds
            disp([ 'The mean b value (and associated t and p vals) for predictor ' termlab{j} ' = ' num2str(bvals(whichsim,j+1)) ' , ' num2str(tvals(whichsim,j+1)) ' , ' num2str(tpvals(whichsim,j+1)) ]); %+1 because of intercept
        end
        whichsim=-1;
    elseif whichsim==0
        %no info given
        checkout=1;
    end

end

%save specific simulation
%@@@could add a save all data option
wantsave=-1;
disp(' ');
disp('Data-saving options:-');
disp('0 do not save any data');
disp('1 save data from specific simulation');
while wantsave<0 ||wantsave >1
    wantsave=input('Type number for selection required, then <Enter> ');
    if isempty(wantsave)
        wantsave=0;
    end
end   
    
if wantsave==1
    %select a specific simulation for saving
    whichsim2save=-1;
    while (whichsim2save<0 || whichsim2save>nsims)
        whichsim2save=input([ 'Type number of simulation (1-' num2str(nsims) ') to save data from, followed by <Enter>; or 0 to skip saving ' ]);
        if isempty(whichsim2save)
            whichsim2save=0;
        end
        if whichsim2save==0
            wantsave=0;
        end
    end
end

if wantsave==1
    
    disp([ 'Saving data for simulation no. ' num2str(whichsim2save) ' ...' ]);
    %@@@could save a file (.MAT type) containing the key simulation parameters
    %that were used to generate the file
    
    %@@@could use table and writetable to shorten this section
    %example below, although have to write out all the vars one by one I think
    %t = table(group,meas(:,1),meas(:,2),meas(:,3),meas(:,4),'VariableNames',{'Gender','t1','t2','t3','t4'});
    %writetable(t,'testit.xlsx','Filetype','spreadsheet')
    
    %create sim number and case number variable
    sim_num=whichsim2save.*ones(n,1);
    case_num=1:n;
    case_num=case_num';
   
    %merge data to save into mydata array
    mydata=[sim_num, case_num, record_data(:,:,33)];
    
    %create variable labels for saved datafile
    nc=2+nvars; %number of columns in saved data file
    vlabs=cell(1,nc);
    vlabs{1}='sim_num';
    vlabs{2}='case_num';
    for j=1:nvars
        vlabs{j+2}=vnames{j};
    end
    %write the data into a cell array, alldata, with variable labels, vlabs, in the
    %first row and saved part of mydata in subsequent rows
    %we use a cell array so that we can save text labels along with numerical
    %data when we write to Excel file
    %set up empty cell string array, with the right dimensions
    if wantsave==1
        alldata=cell(n+1,nc); %+1 for labels row
        alldata(1,:)=vlabs; %add labels to first row
        data_as_cell=mat2cell(mydata,ones(1,n),ones(1,nc)); %#ok<MMTC> %convert mydata to a cell array of same dimensions
        alldata(2:n+1,:)=data_as_cell;
    end
    
    %save mydata to Excel
    %first select file
    DataOutDirectory=pwd;  %@@@might put ths under user control at top
    defaultfilename='my_sim_reg.xlsx';
    [myoutfilename, myoutpathname, filtindex ] = myfileselecta(defaultfilename, DataOutDirectory, 'save');
    
    %now save fitting data to selected file, if we chose a filename
    if filtindex==1
        xlswrite([myoutpathname myoutfilename],alldata); %myoutpathname ends in \ so no need to add \ between path and filenames
    end
else
    disp('No data were saved');
end

%end of programme; function follows
%**********************************

function [filename, pathname, filterindex ] = myfileselecta(default_file, startpath, in_or_out)
    %this just selects a file

    cd(startpath);

    filterindex=0;
    while filterindex==0
    
        if strcmp(in_or_out,'save')
            [filename, pathname, filterindex] = uiputfile( '*.xlsx',  'Select File for Saving Simulated Data', default_file);
        elseif strcmp(in_or_out,'read')
            [filename, pathname, filterindex] = uigetfile( '*.xlsx',  'Select File for Saving Simulated Data', default_file);
        end
        
        if filterindex==0
            menu('You did not select a filename','ok');
            filterindex=2;
        end
    
    end

end

