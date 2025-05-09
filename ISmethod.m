% this script is implementing the iterated simulations (IS) method for
% recovering the diffusive coefficient from a Brownian confined trajectory (original dataset)
% removing errors associated to non-harmonicity of the confining potential.
% For this, simulations are made with different diffusion coefficients
% that the user has to define (Dstore parameter) within a reasonable range
% of values.
% For each element of Dstore, a set of brownian dynamics simulations is made. The  output
% is a plot of the fitted D from both the original dataset and the set of
% simulations.
% the output diffusion coefficient has to be identified as the input
% element of Dstore so that the dataset-simulation discrepancy in the 
% fitted diffusion coefficients are minimize. An interpolation can be made 
% to identify the ottimal D. Looking at the diffusion plot, the proper diffusion coefficient
% is te value of D_{in} where "D_{fit} original dataset" and "D_{fit} IS
% method" intersect

% ORIGINAL DATASET
% in the current version, the original dataset is simulated using an 
% asymmetric potential as defined in the reference paper. 
% For using an experimental dataset, replace the
% 'generate original dataset' section with the upload of the experimental
% dataset. 


%% parameters

scriptfolder = '/Users/svilla/MEGAsync/Montpellier_backup/paper_drag_simulations';
mainfolder = '/Users/svilla/Desktop/data/Montpellier/recover_D_simul/test_github';
savename = [mainfolder,'simulated_trajectory.mat'];



% parameters for simulating the original dataset 
asym_data = 0.8;
Ud = 1;%kB*T;
hw0 = 1;
Noriginal = 2e7*hw0;
D2=1;%1.76;

% parameters for the IS method's simulations
Nsim = 500; % number of simulations: should be large enough so that the std.dev. of the fitted D converges
N = 2e6*hw0;
Dstore = (0.8:0.025:1.2)'; % diffusion coefficients range to explore during simulations


cd(scriptfolder)
%% generate original dataset
plotpotential_and_msd = 1;
clr = [0.5 0.5 0];


b = Ud./(asym_data*hw0);


dtD2_over_hw=11.1e-6; %cioe dt=5e-5 per D2=1.76 e hw=0.5
dt = dtD2_over_hw*hw0/D2;



[xteo,U1]=Uasym_vs_x(b,Ud,hw0);
  

w = normrnd(0,1,Noriginal,1);
x = zeros(Noriginal,1);

para = [b,Ud,hw0];


for t=2:Noriginal
    x(t) = x(t-1) - (Fasym(x(t-1),para)*D2)*dt + sqrt(2*D2*dt)*w(t);
end

data.original.x = x;
data.original.Uteo = U1;
data.original.xteo = xteo;
data.original.Dreal = D2;

% if real dataset has to be used, replace 'data' variable with a struct
% element of same structure containing the real data.
% NB: pay attention to properly make non-dimensional the quantities



%% get MSD original dataset 


D_guess = D2;

x = data.original.x;

decim = 2000*hw0;
dtmsd = decim*dt;
Nmsdo = Noriginal/decim;
Nmsd = N/decim;

MSD_original = myMsd(x(1:decim:end),floor(Nmsdo));
MSD_err = MSD_original./sqrt(length(MSD_original)./(1:length(MSD_original))');
tau = (0:length(MSD_original)-1)';

k=0;
drag_guess  = 1/D_guess;
fig_msd=figure;
[msd_fit,cv] = myFit_MSDconfined(tau*dtmsd,MSD_original,hw0,k,drag_guess,plotpotential_and_msd,fig_msd,clr);

k_fit = 2/cv(1);
drag_fit = k_fit/cv(2);
drag_fit_dataset = drag_fit;
data.original.Dfit = 1/drag_fit_dataset;


%% get potential and force from the original dataset
% using Boltzmann relation

[ap1,ap2]=histcounts(x);
bin=(ap2(2:end)+ap2(1:end-1))/2;
UkBT = -log(ap1/max(ap1)); 

U_smooth = movmean(UkBT,15);
xteo = linspace(10*min(bin),10*max(bin),1000);

bin2 = [-1e3,bin,1e3];
U_smooth2 = [1e15,U_smooth,1e15];
Uteo = interp1(bin,U_smooth,xteo,'spline','extrap');
if Uteo(1)<0 || Uteo(end)<0
    Uteo = interp1(bin2,U_smooth2,xteo,'spline','extrap');
end

figU = figure; hold on
plot(bin,UkBT,'o')
plot(bin,U_smooth,'x')
plot(xteo,Uteo)

Fteo = zeros(size(Uteo))+NaN;
Fteo(2:end-1) = (Uteo(3:end)-Uteo(1:end-2))/(xteo(3)-xteo(1));

figure, plot(xteo,Fteo)
title('F from ')
%% simulation different for D

F_sim_store = cell(length(Dstore),1);
drag_fit_sim = cell(length(Dstore),1);

x = data.original.x; 
data.simulation.Dfit = cell(length(Dstore),1);

for kd=1:length(Dstore)
    D2=Dstore(kd);

    x_sim = zeros(N,1);
    drag_fit_sim{kd} = zeros(Nsim,1);
    for ki=1:Nsim
        myPrintProgressLine(['D: ',num2str(kd),'/',num2str(length(Dstore)),', status: '],ki,Nsim)
        w = normrnd(0,1,N,1);
        Fasym_store = zeros(N,1);
        noise_store = zeros(N,1);
        
        for t=2:N
            
            [~,inde]=min(abs(x_sim(t-1)-xteo));
            Fmeas = Fteo(inde);

            x_sim(t) = x_sim(t-1) - (Fmeas*D2)*dt + sqrt(2*D2*dt)*w(t);
            
        end
        
        % MSD: evaluate and fit
        x = x_sim;

        MSD_sim = myMsd(x(1:decim:end),floor(Nmsd));
        MSD_err = MSD_sim./sqrt(length(MSD_sim)./(1:length(MSD_sim))');
        tau = (0:length(MSD_sim)-1)';
        
        k=0;
        drag_guess  = 1/D2;
        fig_msd=0;
        [msd_fit,cv] = myFit_MSDconfined(tau*dtmsd,MSD_sim,hw0,k,drag_guess,0,fig_msd,clr);
        
        k_fit = 2/cv(1);
        drag_fit = k_fit/cv(2);
        
        drag_fit_sim{kd}(ki) = drag_fit;
        

    end
    
    figure, plot(1./drag_fit_sim{kd},'o')
    hold on
    plot(zeros(Nsim,1)+(1/drag_fit_dataset),'-')
    ylabel('D from fit')
    xlabel('sim. iteration')
    legend({'D_{fit} simulations';'D_{fit} original dataset'})
    title(['input D=',num2str(Dstore(kd))])

    Dfitsim = 1./drag_fit_sim{kd};
    data.simulation.Din(kd) = D2;
    data.simulation.Dfit{kd} = Dfitsim;
 
    save(savename,'data');
    %% plot std.dev
    
    sigma = zeros(Nsim,1);
    for ki=1:Nsim
        app = Dfitsim(1:ki);
        sigma(ki) = std(app);
    end
    figure, plot(sigma)
    title('Dev.Std. D_{fit}')
end

%%
Dinput = data.simulation.Din;
Dfit_OS = zeros(size(Dinput))+data.original.Dfit;
Dfit_IS = zeros(size(Dinput));

for kd=1:length(Dinput)
    Dfit_IS(kd) = mean(data.simulation.Dfit{kd});
end

p = polyfit(Dinput,Dfit_IS,1);
Dfit_IS_lin = Dinput*p(1)+p(2);
D_out = (data.original.Dfit - p(2))/p(1);
figure, hold on
plot(Dinput,Dfit_OS)
plot(Dinput,Dfit_IS,'o')
plot(Dinput,Dfit_IS_lin,'--')
xlabel('D_{in}')
ylabel('D_{fit}')
legend({'D_{fit} original dataset';'D_{fit} IS method';'lin.fit of D_{fit} IS method'})
title(['output D = ',num2str(D_out)])


%% functions

function msd=myMsd(x,maxtau)

    msd = zeros(maxtau,1);
    for tau = 1:maxtau
        msd(tau) = mean((x(tau:end)-x(1:end-tau+1)).^2);
    end
end

function [msd_fit,cv] = myFit_MSDconfined(tau,msd,hw_bin,k_bin,drag,to_plot,figref,colore)

    msd_err = msd./sqrt(tau(end)./tau);    
    if to_plot==1
        figure(figref), hold on
        errorbar(tau,msd,msd_err,'o','color',colore,'markerfacecolor',colore)
    end
    msd = msd/hw_bin^2;
    msd_err = msd./sqrt(tau(end)./tau);
    fitfun = fittype(@(a,b,x) a*(1-exp(-b*x)));
    cv0 = [1,k_bin/drag];
    
    [fitted_curve,~] = fit(tau(2:end),msd(2:end),fitfun,'Weight',msd_err(2:end).^(-2),'StartPoint',cv0,'lower',[0 0]);
    cv = coeffvalues(fitted_curve);
    cv(1) = cv(1)*hw_bin^2;
    msd_fit = cv(1)*(1-exp(-cv(2)*tau));
    if to_plot==1
        plot(tau,msd_fit,'color',colore)
        set(gca,'xscale','log','yscale','log')
    end
end


function [x,U1]=Uasym_vs_x(b,Ud,hw0)


    a = (b^3*hw0^4/Ud^2) - (b*hw0^2/2) + (Ud^2/(16*b));
    x = linspace(-1-sqrt(a/b),100*hw0+sqrt(a/b),50000)';
    U = a./(x+sqrt(a/b)) + b*(x+sqrt(a/b)) - 2*sqrt(a*b);
    F = -a./(x+sqrt(a/b)).^2 + b;
    % prevent divergence in x=-x0
    Ulim=5*Ud;
    xlim = (2*sqrt(a*b)+Ulim-sqrt(Ulim^2+(4*Ulim*sqrt(a*b))))/(2*b);
    indi1 = find(x+ sqrt(a/b)<xlim);
    m = F(indi1(end)+1);
    U1=U; 
    U1(indi1)=m*(x(indi1)-x(indi1(end)+1)) + U(indi1(end)+1);

  
end

    


function F= Fasym(x,para)

    b = para(1);
    Ud = para(2);
    hw0=    para(3);
    a = (b^3*hw0^4/Ud^2) - (b*hw0^2/2) + (Ud^2/(16*b));
    x=x+sqrt(a/b);
    Ulim=5*Ud;
    xlim = (2*sqrt(a*b)+Ulim-sqrt(Ulim^2+(4*Ulim*sqrt(a*b))))/(2*b);
    
    if x>=xlim
    U = a./(x) + b*(x) - 2*sqrt(a*b);
    F = -a./(x).^2 + b;
    else
        % prevent divergence in x=-x0
        Fdx = -a./(xlim).^2 + b;
        Udx = a./(xlim) + b*(xlim) - 2*sqrt(a*b);
        m=Fdx;
        U = m*(x) + Udx;
        F = m;
    end

end


function myPrintProgressLine(namestep,k,imax)

    if k == 1
    
        prog = ( 100*(k/imax) );
        fprintf('\n'); % To go to a new line
        fprintf(1,['Computation Progress ',namestep,':    %3d%%\n',prog]);
    else
    
	    prog = ( 100*(k/imax) );
	    fprintf(1,'\b\b\b\b%3.0f%%',prog); pause(0.01);
    end

end