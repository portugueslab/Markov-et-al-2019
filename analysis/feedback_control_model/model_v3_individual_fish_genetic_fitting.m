clc; clear all;
addpath(genpath('C:\Markov_et_al_2021_Nat_Commun_data&code\analysis_code\feedback_control_model'));
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_model=[pathname_MS_data 'feedback_control_model\'];
cd(pathname_model);

%% ground truth (load training data)
dt=0.01;
gt_train=load('gt.mat','gt_train','lat_train');
lat_train=gt_train.lat_train;
gt_train=gt_train.gt_train;
n_fish=length(lat_train);
weights=ones(1,36);
weights(24)=5; % open loop interbouts
weights(23)=5; %lag 300 interbouts
weights(15)=5; % gain 2 bouts

%% parameters and ranges
% par(1) - wf      - weight between forward motion sensor and sensory integrator
% par(2) - wr      - weight between reverse motion sensor and sensory integrator
% par(3) - taus    - time constant of sensory integrator [s]
% par(4) - wi      - weight between motor integrator and motor output generator
% par(5) - ws      - weight of feed-forward self-excitation of motor output command cell
% par(6) - t       - threshold of motor output command
% par(7) - wm      - weight between motor output command cell and motor integrator
% par(8) - taum    - time constant of motor integrator [s]
list_of_pars={'wf','wr','taus','wi','ws','t','wm','taum'};
n_pars=length(list_of_pars);

%% optimization parameters
min_n_gen=5;     % number of generations
pop_size=100000;  % population size
surv_rate=0.01;  % fraction of survivers
mut_rate=0.1;   % mutation rate

% number of survivers
n_surv=pop_size*surv_rate;
n_offsrping=pop_size/n_surv;

%% output
best_par=nan(n_pars,n_fish);

%% name
my_name=[datestr(datetime('now'),'yymmddHHMM') '_fitting_v3_n_fish_' num2str(n_fish) '_min_n_gen_' num2str(min_n_gen) '_pop_size_' num2str(pop_size) '_surv_rate_' num2str(surv_rate) '_mut_rate_' num2str(mut_rate)];

%% display stuff
disp(['********* Starting evolution of ' num2str(n_fish) ' fish *********']);
disp(['* Population size:                        ' num2str(pop_size)]);
disp(['* Min number of generations:              ' num2str(min_n_gen)]);
disp(['* Survival rate:                          ' num2str(surv_rate*100) '%']);
disp(['* Mutation rate:                          ' num2str(mut_rate*100) '%']);
disp(['* Number of survivers in each generation: ' num2str(n_surv)]);
disp(['* Number of offspring from each surviver: ' num2str(n_offsrping)]);
disp('**************************************************');

%% evolve populations one by one in parallel
lat_range=1;
tlat=0.01;
tmot=0.01;
for f=1:n_fish
    disp(['Fish ' num2str(f) ' / ' num2str(n_fish)]);
    reaf=load('gt.mat','reaf');
    reaf=reaf.reaf';
    for i=[2 4 5]
        reaf(:,i)=round(reaf(:,i)/dt);
    end
    gt=gt_train(f,:);
    weights0=weights./gt;
    
    lat_min=max(lat_train(f)-lat_range,0.25);
    lat_max=min(lat_train(f)+lat_range,10);
    
    best_par_this_fish=nan(8,1);
    best_perf=[];
    % original parameter ranges
    min_taus=0.005;
    max_taus=10;
    min_taum=0.005;
    max_taum=10;
    min_wi=0;
    max_wi=10;
    min_t=0;
    max_t=1;
    min_wr=0;
    min_ws=0;
    min_wm=0;
    % generate original population (independent parameters)
    taus=rand(1,pop_size).*(max_taus-min_taus)+min_taus;
    wi=rand(1,pop_size).*(max_wi-min_wi)+min_wi;
    t=rand(1,pop_size).*(max_t-min_t)+min_t;
    taum=rand(1,pop_size).*(max_taum-min_taum)+min_taum;
    % generate dependent parameters
    min_wf=compute_wf_v3(t, taus, lat_max);
    max_wf=compute_wf_v3(t, taus, lat_min);
    wf=rand(1,pop_size).*(max_wf-min_wf)+min_wf;
    max_wr=compute_max_wr_v3(taus, tlat);
    wr=rand(1,pop_size).*(max_wr-min_wf)+min_wr;
    ws=rand(1,pop_size).*(t-min_ws)+min_ws;
    max_wm=compute_wm_v3(taum, tmot);
    wm=rand(1,pop_size).*(max_wm-min_wm)+min_wm;
    
    % create a population
    this_pop=[wf;wr;taus;wi;ws;t;wm;taum];
    % convert taus to unitless
    this_pop(3,:)=dt./this_pop(3,:);
    this_pop(8,:)=dt./this_pop(8,:);
            
    process1=true;
    process2=true;
    gen=0;
    super_gen=0;
    super_gen_switch_id=0;
    
    while process1
        super_gen=super_gen+1;
        super_gen_switch_id(super_gen)=gen;
        while process2
            gen=gen+1;
            perf=ones(1,pop_size)*inf;
            for i=1:pop_size
                data=zeros(1,36);
                for j=1:18
                    [data(j), data(j+18)] = model_compute_parameters_v3 (model_short_trial_v3(this_pop(:,i),dt,reaf(j,:)),dt);
                end
                perf(i) = sum(abs(gt-data).*weights0);
            end
            % convert taus back to seconds
            this_pop(3,:)=dt./this_pop(3,:);
            this_pop(8,:)=dt./this_pop(8,:);
            % find best performers
            [~,surv]=sort(perf, 'ascend');
            surv=surv(1:n_surv);
            best_par_this_fish=this_pop(:,surv(1));
            best_perf(gen)=min(perf);
            % exit the inner loop
            if gen-super_gen_switch_id(super_gen)>=min_n_gen
                if all(best_perf((gen-min_n_gen+1):gen)==best_perf(gen))
                    process2=false;
                    super_gen_switch_id(super_gen)=gen;
                    do_the_breeding=false;
                else
                    do_the_breeding=true;
                end
            else
                do_the_breeding=true;
            end
            if do_the_breeding
                % let them breed and mutate
                    mut_pars=rand(n_pars,pop_size)>=mut_rate;
                    mut_pars(:,1:n_offsrping:pop_size)=true; % to make sure that at least one kid doesn't mutate and population doesn't exstinct
                    % generate new random population (independent parameters)
                    taus=rand(1,pop_size).*(max_taus-min_taus)+min_taus;
                    wi=rand(1,pop_size).*(max_wi-min_wi)+min_wi;
                    t=rand(1,pop_size).*(max_t-min_t)+min_t;
                    taum=rand(1,pop_size).*(max_taum-min_taum)+min_taum;
                    % generate dependent parameters
                    min_wf=compute_wf_v3(t, taus, lat_max);
                    max_wf=compute_wf_v3(t, taus, lat_min);
                    wf=rand(1,pop_size).*(max_wf-min_wf)+min_wf;
                    max_wr=compute_max_wr_v3(taus, tlat);
                    wr=rand(1,pop_size).*(max_wr-min_wf)+min_wr;
                    ws=rand(1,pop_size).*(t-min_ws)+min_ws;
                    max_wm=compute_wm_v3(taum, tmot);
                    wm=rand(1,pop_size).*(max_wm-min_wm)+min_wm;
                    new_pop=[wf;wr;taus;wi;ws;t;wm;taum];
                    % put non-mutated parametrs from the old population
                    count=0;
                    for i=1:n_surv
                        for ii=1:n_offsrping
                            count=count+1;
                            new_pop(mut_pars(:,count),count)=this_pop(mut_pars(:,count),surv(i));
                        end
                    end
                    this_pop=new_pop;
                    % convert taus to unitless
                    this_pop(3,:)=dt./this_pop(3,:);
                    this_pop(8,:)=dt./this_pop(8,:);
            end
        end
        
        % if last min_n_gen gen were the same we appear here!
        if super_gen>=min_n_gen
            if all(best_perf(super_gen_switch_id((super_gen-min_n_gen+1):super_gen))==best_perf(gen))
                process1=false;
                do_the_breeding=false;
            else
                do_the_breeding=true;
            end
        else
            do_the_breeding=true;
        end
        if do_the_breeding
            % let them breed and mutate
                mut_pars=rand(n_pars,pop_size)>=mut_rate;
                mut_pars(:,1:n_offsrping:pop_size)=true;
                % narrow the ranges by 25 %
                %         'wf','wr','taus','wi','ws','t','wm','taum';
                min_taus=min_taus+0.25*(best_par_this_fish(3)-min_taus);
                max_taus=max_taus-0.25*(max_taus-best_par_this_fish(3));
                min_taum=min_taum+0.25*(best_par_this_fish(8)-min_taum);
                max_taum=max_taum-0.25*(max_taum-best_par_this_fish(8));
                min_wi=min_wi+0.25*(best_par_this_fish(4)-min_wi);
                max_wi=max_wi-0.25*(max_wi-best_par_this_fish(4));
                min_t=min_t+0.25*(best_par_this_fish(6)-min_t);
                max_t=max_t-0.25*(max_t-best_par_this_fish(6));
                min_wr=min_wr+0.25*(best_par_this_fish(2)-min_wr);
                min_ws=min_ws+0.25*(best_par_this_fish(5)-min_ws);
                min_wm=min_wm+0.25*(best_par_this_fish(7)-min_wm);
                % generate new random population (independent parameters)
                if max_taus-min_taus>0 && ...
                        max_wi-min_wi>0 && ...
                        max_t-min_t>0 && ...
                        max_taum-min_taum>0
                    taus=rand(1,pop_size).*(max_taus-min_taus)+min_taus;
                    wi=rand(1,pop_size).*(max_wi-min_wi)+min_wi;
                    t=rand(1,pop_size).*(max_t-min_t)+min_t;
                    taum=rand(1,pop_size).*(max_taum-min_taum)+min_taum;
                    % generate dependent parameters
                    min_wf=compute_wf_v3(t, taus, lat_max);
                    max_wf=compute_wf_v3(t, taus, lat_min);
                    wf=rand(1,pop_size).*(max_wf-min_wf)+min_wf;
                    max_wr=compute_max_wr_v3(taus, tlat);
                    wr=rand(1,pop_size).*(max_wr-min_wf)+min_wr;
                    ws=rand(1,pop_size).*(t-min_ws)+min_ws;
                    max_wm=compute_wm_v3(taum, tmot);
                    wm=rand(1,pop_size).*(max_wm-min_wm)+min_wm;
                    new_pop=[wf;wr;taus;wi;ws;t;wm;taum];
                    % put non-mutated parametrs from the old population
                    count=0;
                    for i=1:n_surv
                        for ii=1:n_offsrping
                            count=count+1;
                            new_pop(mut_pars(:,count),count)=this_pop(mut_pars(:,count),surv(i));
                        end
                    end
                    this_pop=new_pop;
                    % convert taus to unitless
                    this_pop(3,:)=dt./this_pop(3,:);
                    this_pop(8,:)=dt./this_pop(8,:);
                    process2=true;
                else
                    process2=false;
                    process1=false;
                end
        end
    end
    best_par(:,f)=best_par_this_fish;
end
toc
save([pathname_model '\fitting\' my_name '.mat'],'best_par');