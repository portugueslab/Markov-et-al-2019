%% initial stuff
clc; close all; clear all;
addpath(genpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code'));
% folders
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data_code\data\';
pathname_behavior=[pathname_MS_data 'behavior\'];
pathname_model=[pathname_MS_data 'feedback_control_model\'];
pathname_PC_imaging=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_reference_brains = [pathname_MS_data 'reference_brain_stacks\'];
pathname_whole_brain_imaging_integrators = [pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_whole_brain_imaging_lta = [pathname_MS_data 'whole_brain_imaging_long_term_adaptation\'];

% colors used in all figures
col_motor=[0 0.5828 0];
col_sensory=[0.7 0 0.7];
col_activity=[255 127 0]/255;
col_sensors = [0 0.3730 1];
col_integrators = [0.8975 0 0.2390];
temp_fig=figure;
aa=colormap(parula(37));
col_lt_trials=[[0 0 0]; aa(1:21,:); [1 0.7 0.2]];
close(temp_fig); clear temp_fig aa;
% colormap(col_lt_trials);
% colorbar;
col_pre=col_lt_trials(1,:);
col_post=col_lt_trials(end,:);
col_post_start=[255 220 100]/255;
col_post_end=[255 130 0]/255;
col_adapt_start=col_lt_trials(2,:);
col_adapt_end=col_lt_trials(end-1,:);
col_adapt=col_lt_trials(12,:);
col_norm_reaf=[0 0 0];
col_lag_trained=[0.8 0 0];
col_lag_trained_bad=[0.8 0.6 0.6];
col_red_line=[0.8 0 0];
col_model_fitting=[0 0.6 0.6];
col_PC=[1 0.5 0];
r=linspace(col_sensory(1),1,6)';
g=linspace(col_sensory(2),0.7,6)';
b=linspace(col_sensory(3),0,6)';
cols_lags=[r g b];
r=linspace(0, col_sensory(1),4)';
g=linspace(0.7, col_sensory(2),4)';
b=linspace(1, col_sensory(3),4)';
r=[r(1:end-1); linspace(col_sensory(1),1,4)'];
g=[g(1:end-1); linspace(col_sensory(2),0.7,4)'];
b=[b(1:end-1); linspace(col_sensory(3),0,4)'];
cols_gains=flip([r g b],1);
r=linspace(col_sensory(1),0,4)';
g=linspace(col_sensory(2),0.7,4)';
b=linspace(col_sensory(3),1,4)';
r=[linspace(col_sensory(1),1,5)'; r(2:end)];
g=[linspace(col_sensory(2),0.7,5)'; g(2:end)];
b=[linspace(col_sensory(3),0,5)'; b(2:end)];
cols_gain_drops=[r g b];
clear r g b;


%% Fig. 1: Closed-loop experimental assay to study optomotor behavior in larval zebrafish
% % Fig. 1a: Fish swimming with respect to observational and fish reference frames
% % created in allustrator
% 
% % Fig. 1b: Tail movement, fish position and velocity with respect to observational and fish reference frames
% figure('name','Fig. 1b');
% cartoon_pad_length=1;
% [cartoon_time, cartoon_bout, cartoon_vigor] = make_cartoon_bout(cartoon_pad_length);
% est_pos = cumsum(cartoon_vigor);
% axes('xlim',[0 2.4],'ycolor','none','xcolor','none'); hold on;
% fill_bout (gca,cartoon_pad_length+0.001,cartoon_pad_length+0.4,-1,0,0.3);
% plot(cartoon_time,cartoon_bout,'color',col_motor);
% plot(cartoon_time, cartoon_vigor,'color',col_sensory);
% plot(cartoon_time, est_pos/200,'color',col_sensory);
% dm_fix_fig_fonts;
% clear cartoon_bout cartoon_vigor est_pos cartoon_pad_length cartoon_time;
% 
% % Fig. 1c: Experimental rig
% % created in illustrator
% 
% % Fig. 1d: Example trial with all traces
% create_single_trial_axes('Fig. 1d');
% acute_TL=load([pathname_behavior 'acute_reaction_experiment/WT_TL_group/pooled_data.mat'],'behavior');
% acute_TL=acute_TL.behavior;
% fish_num=8;
% fish_id=acute_TL(fish_num).fish_id;
% [time_be, tail, grspeed]=get_example_fish(fish_id, pathname_behavior, 'acute_reaction_experiment/WT_TL_group/');
% ex_trial_num=20;
% [this_time, this_tail, this_gr, bout_starts, bout_ends] = get_example_trial_data(ex_trial_num,time_be,tail,grspeed,acute_TL(fish_num));
% for b=1:length(bout_starts)
%     fill_bout (gca,bout_starts(b),bout_ends(b),-55,0,0.3);
% end
% base_gr=zeros(1,length(this_time));
% base_gr(this_time>0 & this_time<15)=10;
% plot(this_time,base_gr,'color',col_sensory);
% plot(this_time,this_tail*3-4,'color',col_motor);
% fish_vel=this_gr;
% fish_vel(this_time<0.1 | this_time>14.9)=10;
% plot(this_time,-fish_vel-30,'color',col_motor);
% plot(this_time,this_gr-55,'color',col_sensory);
% dm_fix_fig_fonts;
% clear b acute_TL fish_num fish_id time_be tail grspeed ex_trial_num this_time this_tail this_gr bout_starts bout_ends base_gr fish_vel;


%% Fig. 2: Acute reaction to unexpected perturbations in visual feedback can be implemented by a feedback controller
% % Fig. 2a: All reafference conditions
% cartoon_pad_length = 0.5;
% [cartoon_time, cartoon_bout, cartoon_vigor, cartoon_swim] = make_cartoon_bout(cartoon_pad_length);
% cartoon_dt=cartoon_time(2)-cartoon_time(1);
% fig_names = {'Fig. 2ai','Fig. 2aii top','Fig. 2aii bottom','Fig. 2aiii top','Fig. 2aiii bottom'};
% reaf_conds = {[[0;0;0;0;0],[0.33;0;0;0;0],[0.66;0;0;0;0],[1;0;0;0;0],[1.33;0;0;0;0],[1.66;0;0;0;0],[2;0;0;0;0]],...
%     [[1;0;0;0;0],[1;0.075/cartoon_dt;0;0;0],[1;0.15/cartoon_dt;0;0;0],[1;0.225/cartoon_dt;0;0;0],[1;0.3/cartoon_dt;0;0;0],[0;0;0;0;0]],...
%     [[1;0;0;0;0],[1;0.075/cartoon_dt;1;0;0],[1;0.15/cartoon_dt;1;0;0],[1;0.225/cartoon_dt;1;0;0],[1;0.3/cartoon_dt;1;0;0],[0;0;0;0;0]],...
%     [[1;0;0;0;0],[1;0;0;0;0.075/cartoon_dt],[1;0;0;0;0.150/cartoon_dt],[1;0;0;0;0.225/cartoon_dt],[1;0;0;0;0.3/cartoon_dt]],...
%     [[1;0;0;0;0],[1;0;0;0.225/cartoon_dt;0.3/cartoon_dt],[1;0;0;0.15/cartoon_dt;0.3/cartoon_dt],[1;0;0;0.075/cartoon_dt;0.3/cartoon_dt],[1;0;0;0;0.3/cartoon_dt]]};
% for i=1:length(fig_names)
%     figure('name',fig_names{i});
%     axes('ycolor','none','xcolor','none'); hold on;
%     plot_reaf_cond(reaf_conds{i}, cartoon_time, cartoon_bout, cartoon_swim, cartoon_vigor, cartoon_pad_length, col_motor, col_sensory)
%     dm_fix_fig_fonts;
% end
% clear cartoon_time cartoon_bout cartoon_swim cartoon_vigor cartoon_pad_length cartoon_dt reaf_conds
% 
% % Fig. 2b: The model (short little circuits for insets)
% dt=0.005;
% par_schema=[0.27; 0.4; 1; 5; 0.6; 0.6; 0.5; 0.6];
% par_schema(3)=dt/par_schema(3);
% par_schema(8)=dt/par_schema(8);
% [swim,grspeed,brain_state] = model_one_bout_trial_v3 (par_schema, dt);
% figure('name','Fig. 2b');
% subplot(4,4,1,'xcolor','none','ycolor','none','ylim',[-10 10]); hold on
% plot(grspeed,'color',col_sensory,'linewidth',1.5)
% pbaspect([2 1 1])
% subplot(4,4,2,'xcolor','none','ycolor','none','ylim',[-1 1]); hold on
% plot(swim,'color',col_motor,'linewidth',1.5)
% pbaspect([2 1 1])
% y_lims=[-10 10; -10 10; -1 1; -0.6 0.6; -0.2 0.2];
% for i=1:5
%     subplot(4,4,2+i,'xcolor','none','ycolor','none','ylim',y_lims(i,:)); hold on
%     plot(brain_state(i,:),'color',col_activity,'linewidth',1.5);
%     pbaspect([2 1 1])
% end
% subplot(4,4,2+i+1,'xcolor','none','ycolor','none','ylim',[-1 1]); hold on
% plot(swim,'color',col_activity,'linewidth',1.5);
% pbaspect([2 1 1])
% subplot(4,4,2+i+2,'xcolor','none','ycolor','none','ylim',[-10 10]); hold on;
% grspeed(grspeed<0)=10;
% plot(grspeed,'color',col_sensory,'linewidth',1.5);
% pbaspect([2 1 1]);
% dm_fix_fig_fonts;
% clear par_schema swim grspeed brain_state y_lims i;
% 
% % Fig. 2c-d: bout and interbout duration in acute reaction experiment
% acute_TL=load([pathname_behavior 'acute_reaction_experiment/WT_TL_group/pooled_data.mat'],'behavior');
% acute_TL=acute_TL.behavior;
% reaf=load([pathname_model 'gt.mat'],'reaf');
% reaf=reaf.reaf';
% for i=[2 4 5]
%     reaf(:,i)=round(reaf(:,i)/dt);
% end
% reaf_cond_ids={[6 11 12 1 13 14 15];... % gain
%     [1 2 3 4 5 6];... % lag
%     [1 7 8 9 10 6];... % shunted lag
%     [1 7 8 9 10 16 17 18]}; % gain drop
% best_par=load([pathname_model 'fitted_parameters.mat'],'best_par');
% best_par=best_par.best_par;
% best_par(3,:)=dt./best_par(3,:);
% best_par(8,:)=dt./best_par(8,:);
% panel_names={'Fig. 2ci', 'Fig. 2cii', 'Fig. 2ciii', 'Fig. 2civ'; 'Fig. 2di', 'Fig. 2dii', 'Fig. 2diii', 'Fig. 2div'};
% cond_names={'gain', 'lag', 'shunted_lag', 'gain_drop'};
% par_names={'bout_duration';'next_interbout_duration'};
% x_labels={'gain','lag [ms]','shunted lag [ms]','gain profile'};
% y_labels={'mean bout duration [s]'; 'mean interbout duration [s]'};
% y_lims=[0.25 0.55; 0.8 1.7];
% for i=1:2
%     for c=1:4
%         prepare_acute_axes(panel_names{i,c},acute_TL,cond_names{c},x_labels{c});
%         data_exp=prepare_data_acute(acute_TL,cond_names{c},par_names{i});
%         data_mod=prepare_data_acute_model(reaf(reaf_cond_ids{c},:),best_par,dt,par_names{i});
%         plot_acute({data_exp,data_mod},[0 0 0; col_model_fitting],y_labels{i},y_lims(i,:));
%         dm_fix_fig_fonts;
%     end
% end
% 
% % Extended Data Fig. 1b: fitted model parameters
% best_par=load([pathname_model 'fitted_parameters.mat'],'best_par');
% best_par=best_par.best_par;
% figure('name','Extended Data Fig. 1b');
% titles = {'wf','wr','taus','wi','ws','t','wm','taum'};
% my_edges = {0:1/10:1, 0:1/10:1, 0:1:10, 0:1:10, 0:1/10:1, 0:1/10:1, 0:4/10:4, 0:1:10};
% for i=1:8
%     subplot(2,4,i,'xtick',my_edges{i},'ylim',[0 40]); hold on;
%     title(titles{i});
%     a = best_par(i,:);
%     [b, tf] = rmoutliers(a);
%     disp([titles{i} ': outliers removed: ' num2str(a(tf))]);
%     histogram(b,my_edges{i},'facecolor',[0.5 0.5 0.5]);
%     drawnow;
%     x_tick_label = get(gca,'xticklabel');
%     for j=2:5:length(my_edges{i})
%         x_tick_label{j}='';
%         x_tick_label{j+1}='';
%         x_tick_label{j+2}='';
%         x_tick_label{j+3}='';
%     end
%     set(gca,'xticklabel',x_tick_label);
% end
% dm_fix_fig_fonts;
% 
% % Fig. 2e: bout power under different reafference conditions
% time_power=-0.1:dt:1-dt;
% h_acute = find_ballistic_end_acyte(acute_TL,[4 18],time_power);
% panel_names={'Fig. 2ei','Fig. 2eii','Fig. 2eiii','Fig. 2eiv'};
% cond_names={'gain', 'lag', 'shunted_lag', 'gain_drop'};
% cols={cols_gains, cols_lags, cols_lags, cols_gain_drops};
% for i=1:4
%     plot_bout_power_acute(panel_names{i},acute_TL,cond_names{i},cols{i},time_power,h_acute,[0 0 0],[0 2.5]);
%     dm_fix_fig_fonts;
% end
% clear time_power dt fig_names data_mod data_exp panel_names cond_names cols i acute_TL h_acute ballistic_end_acute_gaindrop0000 ballistic_end_acute_gaindrop0001 ballistic_end_acute_gaindrop0011 ballistic_end_acute_gaindrop0111 ballistic_end_acute_gaindrop1000 ballistic_end_acute_gaindrop1100 ballistic_end_acute_gaindrop1110 acute_TL best_par i c reaf reaf_cond_ids panel_names cond_names par_names x_labels y_labels y_lims;


%% Fig. 3: Larval zebrafish are able to integrate the optic flow
% % load the data
% [sz, rez] = read_nrrd_metadata([pathname_reference_brains 'PortuguesLab_wholebrain_ref_gamma_masked.nrrd']);
% all_fish = dm_dir([pathname_whole_brain_imaging_integrators 'behavior\*_f*_behavior.mat']);
% all_fish = strrep(all_fish,'_behavior.mat','');
% n_fish = length(all_fish);
% ROI_coord_all = cell(n_fish,1);
% n_ROIs_all = zeros(n_fish,1);
% sensmot_clust_all = cell(n_fish,1);
% traces_gr_trig_all = cell(n_fish,1);
% traces_bout_trig_all = cell(n_fish,1);
% traces_gr_trig_tau_all = cell(n_fish,1);
% time_constants_all = cell(n_fish,1);
% load([pathname_whole_brain_imaging_integrators 'triggered_traces\' all_fish{1} '_trig_traces.mat'],'time_trig','time_trig_tau')
% progressbar('Loading data for Fig. 3...');
% for f = 1:n_fish
%     fish_id = all_fish{f};
%     load([pathname_whole_brain_imaging_integrators 'ROIs\' fish_id '_ROIs.mat'],'ROI_coord');
%     ROI_coord_all{f} = ROI_coord;
%     n_ROIs_all(f) = length(ROI_coord);
%     load([pathname_whole_brain_imaging_integrators 'triggered_traces\' fish_id '_trig_traces.mat'],'traces_gr_trig_mean','traces_bout_trig_mean','traces_gr_trig_mean_tau');
%     traces_gr_trig_all{f} = traces_gr_trig_mean;
%     traces_bout_trig_all{f} = traces_bout_trig_mean;
%     traces_gr_trig_tau_all{f} = traces_gr_trig_mean_tau;
%     load([pathname_whole_brain_imaging_integrators 'clustering\' fish_id '_clustering.mat'],'sensmot_clust');
%     sensmot_clust_all{f} = sensmot_clust;
%     load([pathname_whole_brain_imaging_integrators 'time_constants\' fish_id '_time_constants.mat'],'time_constants');
%     time_constants_all{f} = time_constants;
%     progressbar(f/n_fish);
% end
% my_colormap=hot(64);
% my_colormap=my_colormap(round(linspace(1,64,n_fish+1)),:);
% ex_fish_id = '190523_f0';
% ex_fish_num = find(contains(all_fish,ex_fish_id));
% ex_trial = 45;
% ex_sens_ROI_id = 21644;
% ex_mot_ROI_id = 26246;
% ex_trial_int = 43;
% ex_sensor_id = 21645;
% ex_integrator_id = 26774;
% clear all_fish fish_id ROI_coord traces_gr_trig_mean traces_gr_trig_mean_tau sensmot_clust time_constants
% % Fig. 3a: Schematics of the light-sheet setup
% % created in illustrator
% 
% % Fig. 3b: colormap of all imaged ROIs with location of 4 example ROIs
% % compute stack with fraction of ROIs
% A = zeros(sz,'uint8');
% for f = 1:n_fish
%     for i = 1:n_ROIs_all(f)
%         coord = ROI_coord_all{f}{i};
%         A(coord) = A(coord) + 1;
%     end
% end
% % compute top and side views
% [top_view, side_view]=three_orthogonal_views(A,rez(1));
% top_view=1-top_view/n_fish;
% side_view=1-round(side_view)/n_fish;
% % plot top and side views of the fractions of all ROIs
% figure('name','Fig. 3b top');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/sz(2) 1]);
% imagesc(top_view);
% colormap(my_colormap);
% figure('name','Fig. 3b side');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
% imagesc(side_view);
% colormap(my_colormap);
% % make a colorbar for fish-fraction plots
% figure('name','Fig. 3b colormap');
% axes('ycolor','none','xcolor','none');
% colormap(flip(my_colormap,1));
% clrbr=colorbar;
% clrbr.Ticks=1/(n_fish+1)/2:1/(n_fish+1):1-1/(n_fish+1)/2;
% clrbr.TickLabels=round((0:1/n_fish:1)'*100);
% ylabel(clrbr,{'percentage of fish';'with activity'});
% dm_fix_fig_fonts;
% % plot top and side views of the reference brain mask (for outlines)
% ref_brain_mask = nrrdread([pathname_reference_brains 'PortuguesLab_wholebrain_ref_mask.nrrd']);
% [top_view,side_view] = three_orthogonal_views(ref_brain_mask,rez(2));
% figure('name','Fig. 3b reference brain mask top');
% imshow(1-top_view);
% figure('name','Fig. 3b reference brain mask side');
% imshow(1-side_view);
% % plot 100 um scale bar
% a=ones(sz(1:2));
% a(50:51,50:50+round(100/rez(1)))=0;
% figure('name','Fig. 3b 100 micron scale bar');
% imshow(a);
% % plot top and side views of locations of example ROIs shown in other panels
% A_ex_locations = zeros([sz 3]);
% for i=1:3
%     temp = zeros(sz);
%     temp(ROI_coord_all{ex_fish_num}{ex_sens_ROI_id}) = col_sensory(i);
%     temp(ROI_coord_all{ex_fish_num}{ex_mot_ROI_id}) = col_motor(i);
%     temp(ROI_coord_all{ex_fish_num}{ex_sensor_id}) = col_sensors(i);
%     temp(ROI_coord_all{ex_fish_num}{ex_integrator_id}) = col_integrators(i);
%     A_ex_locations(:,:,:,i) = temp;
% end
% top_view_ex_locations = squeeze(max(A_ex_locations,[],3));
% figure('name','Fig. 3b example roi locations top');
% imshow(top_view_ex_locations);
% side_view_ex_locations=squeeze(double(max(A_ex_locations,[],2)));
% side_view_ex_locations=imresize(side_view_ex_locations,[sz(1),round(sz(3)/rez(1))]);
% side_view_ex_locations=flip(side_view_ex_locations,2);
% figure('name','Fig. 3b example roi locations side');
% imshow(side_view_ex_locations);
% % plot raw anatomy to show that these ROIs are cell bodies
% ex_anatomy = nrrdread([pathname_whole_brain_imaging_integrators 'anatomies\' ex_fish_id '_anatomy.nrrd']);
% figure('name','Fig. 3b example roi locations zoomed in');
% subplot(2,2,1);
% imshow(small_anatomy_image(ex_anatomy,sz,ROI_coord_all,ex_fish_num,ex_sens_ROI_id,col_sensory));
% subplot(2,2,2);
% imshow(small_anatomy_image(ex_anatomy,sz,ROI_coord_all,ex_fish_num,ex_mot_ROI_id,col_motor));
% subplot(2,2,3);
% imshow(small_anatomy_image(ex_anatomy,sz,ROI_coord_all,ex_fish_num,ex_sensor_id,col_sensors));
% subplot(2,2,4);
% imshow(small_anatomy_image(ex_anatomy,sz,ROI_coord_all,ex_fish_num,ex_integrator_id,col_integrators));
% figure('name','Fig. 3b example roi locations zoomed in scale bar 50 um');
% a = small_anatomy_image(ex_anatomy,sz,ROI_coord_all,ex_fish_num,ex_integrator_id,col_integrators);
% a(:,:,:)=0;
% a(4:5,5:5+round(10/rez(1)),:)=255;
% subplot(2,2,1);
% imshow(a);
% clear ex_anatomy A top_view side_view clrbr ref_brain_mask a A_ex_locations temp top_view_ex_locations side_view_ex_locations
% 
% % Fig. 3c: activity of example sensory and motor ROIs in one example trial
% % load example fish data
% load([pathname_whole_brain_imaging_integrators 'behavior\' ex_fish_id '_behavior.mat'],'tail','bouts','time_be','meta');
% load([pathname_whole_brain_imaging_integrators 'traces\' ex_fish_id '_traces.mat'],'traces','time_offsets');
% % compute traces for the example trial
% dt = round((time_be(2)-time_be(1))*1000)/1000;
% trial_starts=120+7.5:30:time_be(end);
% this_trial_starts = trial_starts(ex_trial);
% this_tail = tail(time_be>=this_trial_starts-5 & time_be<this_trial_starts+20);
% this_time_be = -5+dt:dt:20;
% this_grspeed = zeros(1,length(this_time_be));
% this_grspeed(this_time_be>=0 & this_time_be<15)=10;
% bout_starts = bouts.start(bouts.trial==ex_trial)-this_trial_starts;
% bout_ends = bouts.end(bouts.trial==ex_trial)-this_trial_starts;
% dt_im = round(1/meta.fs*1000)/1000;
% time_im = dt_im:dt_im:dt_im*size(traces,2);
% this_trace_mot = traces(ex_mot_ROI_id,time_im+time_offsets(ex_mot_ROI_id)>=this_trial_starts-5 & time_im+time_offsets(ex_mot_ROI_id)<this_trial_starts+20);
% this_trace_sens = traces(ex_sens_ROI_id,time_im+time_offsets(ex_sens_ROI_id)>=this_trial_starts-5 & time_im+time_offsets(ex_sens_ROI_id)<this_trial_starts+20);
% this_time_im = -5+dt_im:dt_im:20;
% % plot
% figure('name','Fig. 3c');
% axes('ycolor','none','xlim',[-5 20],'xtick',-5:5:20,'ylim',[-8.5 2]); hold on
% plot(this_time_be,this_grspeed/8,'color',col_sensory);
% plot(this_time_be,this_tail/2-1,'color',col_motor);
% plot(this_time_im,this_trace_sens-6,'color',col_sensory);
% plot(this_time_im,this_trace_mot-6,'color',col_motor);
% line([this_time_be(2) this_time_be(2)],[-6-1 -6+1],'color','k'); % line indicating 2sd
% line([0 0],ylim,'color','k','linestyle',':');
% for i=1:length(bout_starts)
%     fill([bout_starts(i) bout_ends(i) bout_ends(i) bout_starts(i)],[0 0 1 1],'k');
% end
% xlabel('time relative to trial onset [s]');
% dm_fix_fig_fonts;
% clear this_trace_sens this_trace_mot meta dt bout_starts bout_ends
% 
% % Fig. 3d: triggered activity of example traces
% load([pathname_whole_brain_imaging_integrators 'triggered_traces\' ex_fish_id '_trig_traces.mat'],'traces_gr_trig_ste','traces_bout_trig_ste');
% figure('name','Fig. 3d');
% subplot(1,2,1,'ylim',[-1 2],'xlim',[-1 4],'xtick',-1:4,'ycolor','none'); hold on
% xlabel('time relative to trigger [s]');
% dm_fill_plot(time_trig,traces_gr_trig_all{ex_fish_num}(ex_sens_ROI_id,:),traces_gr_trig_ste(ex_sens_ROI_id,:),col_sensory,0.2);
% dm_fill_plot(time_trig,traces_gr_trig_all{ex_fish_num}(ex_mot_ROI_id,:),traces_gr_trig_ste(ex_mot_ROI_id,:),col_motor,0.2);
% line([0 0],[-1 2],'color','k','linestyle',':');
% line([-0.5 -0.5],[-0.5 0.5],'color','k'); % line indicating 1 sd
% subplot(1,2,2,'ylim',[-1 2],'xlim',[-1 4],'xtick',-1:4,'ycolor','none'); hold on
% xlabel('time relative to trigger [s]');
% dm_fill_plot(time_trig,traces_bout_trig_all{ex_fish_num}(ex_sens_ROI_id,:),traces_bout_trig_ste(ex_sens_ROI_id,:),col_sensory,0.2);
% dm_fill_plot(time_trig,traces_bout_trig_all{ex_fish_num}(ex_mot_ROI_id,:),traces_bout_trig_ste(ex_mot_ROI_id,:),col_motor,0.2);
% line([0 0],[-1 2],'color','k','linestyle',':');
% dm_fix_fig_fonts;
% clear traces_bout_trig_ste traces_gr_trig_ste
% 
% % Fig. 3e: triggered activity of all sensory and motor ROIs
% % prepare the data
% traces_gr_trig_sens = [];
% traces_gr_trig_mot = [];
% traces_bout_trig_sens = [];
% traces_bout_trig_mot = [];
% for f=1:n_fish
%     traces_gr_trig_sens = [traces_gr_trig_sens; traces_gr_trig_all{f}(sensmot_clust_all{f}==1,:)];
%     traces_gr_trig_mot = [traces_gr_trig_mot; traces_gr_trig_all{f}(sensmot_clust_all{f}==2,:)];
%     traces_bout_trig_sens = [traces_bout_trig_sens; traces_bout_trig_all{f}(sensmot_clust_all{f}==1,:)];
%     traces_bout_trig_mot = [traces_bout_trig_mot; traces_bout_trig_all{f}(sensmot_clust_all{f}==2,:)];
% end
% traces_gr_trig_sens = traces_gr_trig_sens(randperm(size(traces_gr_trig_sens,1)),:);
% traces_gr_trig_mot = traces_gr_trig_mot(randperm(size(traces_gr_trig_mot,1)),:);
% traces_bout_trig_sens = traces_bout_trig_sens(randperm(size(traces_bout_trig_sens,1)),:);
% traces_bout_trig_mot = traces_bout_trig_mot(randperm(size(traces_bout_trig_mot,1)),:);
% data_gr_trig = [traces_gr_trig_sens;traces_gr_trig_mot];
% data_bout_trig = [traces_bout_trig_sens;traces_bout_trig_mot];
% min_val=min([data_gr_trig(:); data_bout_trig(:)]);
% max_val=max([data_gr_trig(:); data_bout_trig(:)]);
% my_colormap_gray=flip(gray);
% to_add=ones(15,3).*linspace(1,my_colormap_gray(2,1),15)';
% my_colormap_gray(1:2,:)=[];
% my_colormap_gray=[to_add; my_colormap_gray];
% % plot
% figure('name','Fig. 3e');
% subplot(1,2,1,'ydir','reverse','xlim',[-1 4],'ylim',0.5+[0,size(data_gr_trig,1)],'ycolor','none','layer','top','clim',[min_val max_val],'xtick',-1:4); hold on;
% xlabel('time relative to trigger [s]');
% imagesc([time_trig(1) time_trig(end)],[1,size(data_gr_trig,1)],data_gr_trig);
% line([time_trig(10) time_trig(10)],0.5+[0 size(traces_gr_trig_sens,1)],'linewidth',4,'color',col_sensory);
% line([time_trig(10) time_trig(10)],0.5+[size(traces_gr_trig_sens,1) size(data_gr_trig,1)],'linewidth',4,'color',col_motor);
% line([0 0],ylim,'color','k','linestyle',':');
% colormap(my_colormap_gray);
% colorbar;
% subplot(1,2,2,'ydir','reverse','xlim',[-1 4],'ylim',[1,size(data_bout_trig,1)],'ycolor','none','layer','top','clim',[min_val max_val],'xtick',-1:4); hold on;
% xlabel('time relative to trigger [s]');
% imagesc([time_trig(1) time_trig(end)],[1,size(data_bout_trig,1)],data_bout_trig);
% line([0 0],ylim,'color','k','linestyle',':');
% colormap(my_colormap_gray);
% colorbar;
% dm_fix_fig_fonts;
% clear my_colormap_gray max_val min_val data_bout_trig data_gr_trig traces_bout_trig_mot traces_bout_trig_sens traces_gr_trig_mot traces_gr_trig_sens
% 
% % Fig. 3f: location of sensory and motor ROIs
% % compute stacks with fraction of ROIs
% A_sens = zeros(sz,'uint8');
% A_mot = zeros(sz,'uint8');
% for f = 1:n_fish
%     for i = 1:n_ROIs_all(f)
%         coord = ROI_coord_all{f}{i};
%         if sensmot_clust_all{f}(i)==1
%             A_sens(coord) = A_sens(coord) + 1;
%         elseif sensmot_clust_all{f}(i)==2
%             A_mot(coord) = A_mot(coord) + 1;
%         end
%     end
% end
% % compute top and side views
% [top_view_sens, side_view_sens]=three_orthogonal_views(A_sens,rez(1));
% top_view_sens=1-top_view_sens/n_fish;
% side_view_sens=1-round(side_view_sens)/n_fish;
% [top_view_mot, side_view_mot]=three_orthogonal_views(A_mot,rez(1));
% top_view_mot=1-top_view_mot/n_fish;
% side_view_mot=1-round(side_view_mot)/n_fish;
% % plot top and side views of the fractions of sensory ROIs
% figure('name','Fig. 3f sensory top');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/sz(2) 1]);
% imagesc(top_view_sens);
% colormap(my_colormap);
% figure('name','Fig. 3f sensory side');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
% imagesc(side_view_sens);
% colormap(my_colormap);
% figure('name','Fig. 3f motor top');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/sz(2) 1]);
% imagesc(top_view_mot);
% colormap(my_colormap);
% figure('name','Fig. 3f motor side');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
% imagesc(side_view_mot);
% colormap(my_colormap);
% clear top_view_mot side_view_mot top_view_sens side_view_sens A_mot A_sens
% 
% % Fig. 3g: activity of example sensor and integrator ROIs in one example trial
% % compute traces for the example trial
% this_trial_starts = trial_starts(ex_trial_int);
% this_tail = tail(time_be>=this_trial_starts-5 & time_be<this_trial_starts+20);
% this_trace_sensor = traces(ex_sensor_id,time_im+time_offsets(ex_sensor_id)>=this_trial_starts-5 & time_im+time_offsets(ex_sensor_id)<this_trial_starts+20);
% this_trace_integrator = traces(ex_integrator_id,time_im+time_offsets(ex_integrator_id)>=this_trial_starts-5 & time_im+time_offsets(ex_integrator_id)<this_trial_starts+20);
% this_time_im = -5+dt_im:dt_im:20;
% % plot
% figure('name','Fig. 3g');
% axes('ycolor','none','xlim',[-5 20],'xtick',-5:5:20,'ylim',[-8.5 2]); hold on
% plot(this_time_be,this_grspeed/8,'color',col_sensory);
% plot(this_time_be,this_tail/2-1,'color',col_motor);
% plot(this_time_im,this_trace_sensor-6,'color',col_sensors);
% plot(this_time_im,this_trace_integrator-6,'color',col_integrators);
% line([0 0],ylim,'color','k','linestyle',':');
% line([this_time_be(2) this_time_be(2)],[-6-1 -6+1],'color','k'); % line indicating 2sd
% xlabel('time relative to trial onset [s]');
% dm_fix_fig_fonts;
% % inset showing distribution of time constants
% edges = 0.5:0.5:10.5;
% temp_fig = figure('name','temp_fig');
% ax = axes;
% data = zeros(n_fish,length(edges)-1);
% for f=1:n_fish
%     temp = time_constants_all{f}(~isnan(time_constants_all{f}));
%     h = histogram(ax,temp,edges,'Normalization','probability');
%     data(f,:)=h.Values;
% end
% close(temp_fig);
% figure('name','Fig. 3g inset');
% axes('xlim',[0 6],'xtick',0:2:6); hold on;
% xlabel('time constant, s');
% ylabel('fraction of sensory ROIs');
% dm_fill_plot(0.25:0.5:9.75,mean(data,1),std(data,0,1)/sqrt(6),[0 0 0],0.2);
% dm_fix_fig_fonts;
% clear edges temp_fig ax data f temp h this_trace_integrator this_trace_sensor this_time_im time_im dt_im this_grspeed this_time_be tail time_be traces time_offsets trial_starts this_trial_starts this_tail
% 
% % Fig. 3h: triggered activity of example traces
% load([pathname_whole_brain_imaging_integrators 'triggered_traces\' ex_fish_id '_trig_traces.mat'],'traces_gr_trig_ste_tau');
% figure('name','Fig. 3h');
% axes('ylim',[-1 3],'xlim',[-1 20],'xtick',[-1 0:5:20],'ycolor','none'); hold on
% xlabel('time relative to trial onset [s]');
% dm_fill_plot(time_trig_tau,traces_gr_trig_tau_all{ex_fish_num}(ex_sensor_id,:),traces_gr_trig_ste_tau(ex_sensor_id,:),col_sensors,0.2);
% dm_fill_plot(time_trig_tau,traces_gr_trig_tau_all{ex_fish_num}(ex_integrator_id,:),traces_gr_trig_ste_tau(ex_integrator_id,:),col_integrators,0.2);
% line([0 0],[-1 3],'color','k','linestyle',':');
% line([15 15],[-1 3],'color','k','linestyle',':');
% line([-0.5 -0.5],[-0.5 0.5],'color','k'); % line indicating 1 sd
% dm_fix_fig_fonts;
% clear traces_gr_trig_ste_tau
% 
% % Fig. 3i: triggered activity of all sensory ROIs sorted based on tau
% % prepare the data
% data_time_constants = [];
% data_gr_trig_tau = [];
% for f=1:n_fish
%     data_time_constants = [data_time_constants; time_constants_all{f}(sensmot_clust_all{f}==1)];
%     data_gr_trig_tau = [data_gr_trig_tau; traces_gr_trig_tau_all{f}(sensmot_clust_all{f}==1,:)];
% end
% [data_time_constants,ids_tau_sorted] = sort(data_time_constants);
% data_gr_trig_tau = data_gr_trig_tau(ids_tau_sorted,:);
% temp = data_gr_trig_tau(data_time_constants<=1.5,:);
% temp = temp(randperm(size(temp,1)),:);
% data_gr_trig_tau(data_time_constants<=1.5,:) = temp;
% temp = data_gr_trig_tau(data_time_constants>1.5,:);
% temp = temp(randperm(size(temp,1)),:);
% data_gr_trig_tau(data_time_constants>1.5,:) = temp;
% % plot
% figure('name','Fig. 3i');
% axes('ydir','reverse','xlim',[-5 20],'ylim',0.5+[0,size(data_gr_trig_tau,1)],'ycolor','none','layer','top','xtick',-5:5:20); hold on;
% xlabel('time relative to trial onset [s]');
% imagesc([time_trig_tau(1) time_trig_tau(end)],[1,size(data_gr_trig_tau,1)],data_gr_trig_tau);
% line([0 0],ylim,'color','k','linestyle',':');
% colormap(flip(gray));
% colorbar;
% line([-2 -2],[0 find(data_time_constants>=0 & data_time_constants<=1.5,1,'last')]+0.5,'color',col_sensors,'linewidth',4);
% line([-2 -2],[find(data_time_constants>=0 & data_time_constants<=1.5,1,'last') size(data_gr_trig_tau,1)]+0.5,'color',col_integrators,'linewidth',4);
% dm_fix_fig_fonts;
% clear temp ids_tau_sorted data_gr_trig_tau data_time_constants
% 
% % Fig. 3j: location of sensors and integrators
% % compute stacks with fraction of ROIs
% A_sensors = zeros(sz,'uint8');
% A_integrators = zeros(sz,'uint8');
% for f = 1:n_fish
%     for i = 1:n_ROIs_all(f)
%         coord = ROI_coord_all{f}{i};
%         if time_constants_all{f}(i)>1.5
%             A_integrators(coord) = A_integrators(coord) + 1;
%         elseif time_constants_all{f}(i)>0
%             A_sensors(coord) = A_sensors(coord) + 1;
%         end
%     end
% end
% % compute top and side views
% [top_view_sensors, side_view_sensors]=three_orthogonal_views(A_sensors,rez(1));
% top_view_sensors=1-top_view_sensors/n_fish;
% side_view_sensors=1-round(side_view_sensors)/n_fish;
% [top_view_integrators, side_view_integrators]=three_orthogonal_views(A_integrators,rez(1));
% top_view_integrators=1-top_view_integrators/n_fish;
% side_view_integrators=1-round(side_view_integrators)/n_fish;
% % plot top and side views of the fractions of sensory ROIs
% figure('name','Fig. 3j sensors top');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/sz(2) 1]);
% imagesc(top_view_sensors);
% colormap(my_colormap);
% figure('name','Fig. 3j sensors side');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
% imagesc(side_view_sensors);
% colormap(my_colormap);
% figure('name','Fig. 3j integrators top');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/sz(2) 1]);
% imagesc(top_view_integrators);
% colormap(my_colormap);
% figure('name','Fig. 3j integrators side');
% axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
% pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
% imagesc(side_view_integrators);
% colormap(my_colormap);
% clear n_fish ROI_coord_all n_ROIs_all sensmot_clust_all traces_gr_trig_all traces_bout_trig_all traces_gr_trig_tau_all ...
%     time_constants_all time_trig time_trig_tau f my_colormap ...
%     ex_fish_id ex_fish_num ex_trial ex_sens_ROI_id ex_mot_ROI_id ex_trial_int ex_sensor_id ex_integrator_id...
%     top_view_integrators side_view_integrators top_view_sensors side_view_sensors coord A_sensors A_integrators

%% Fig. 4: Larval zebrafish adapt their behavior in response to a long-lasting perturbation in visual reafference
% % Fig. 4a: Experimental protocol (insets showing normal and lagged reafference)
% figure('name','Fig. 4a');
% axes('ycolor','none','xcolor','none'); hold on;
% cartoon_pad_length = 0.5;
% [cartoon_time, cartoon_bout, cartoon_vigor, cartoon_swim] = make_cartoon_bout(cartoon_pad_length);
% reaf_cond=[[1;0;0;0;0], [1;0.225/0.001;0;0;0]];
% plot_reaf_cond(reaf_cond, cartoon_time, cartoon_bout, cartoon_swim, cartoon_vigor, cartoon_pad_length, col_motor, col_sensory)
% line([0 0.5],[-20 -20],'color','k');
% dm_fix_fig_fonts;
% clear i cartoon_pad_length cartoon_time cartoon_bout cartoon_vigor cartoon_swim cartoon_gr reaf_cond;
% 
% % Fig. 4b: Example fish: individual trials
% lt_lag_TL=load([pathname_behavior 'long_term_adaptation_experiment/lag_trained/WT_TL_group/pooled_data.mat'],'behavior');
% lt_lag_TL=lt_lag_TL.behavior;
% fish_num=40;
% fish_id=lt_lag_TL(fish_num).fish_id;
% [time_be, tail, grspeed]=get_example_fish(fish_id, pathname_behavior, 'long_term_adaptation_experiment/lag_trained/WT_TL_group/');
% dt=time_be(2)-time_be(1);
% create_single_trial_axes('Fig. 4b');
% ex_trial_num = [19:22 225:226 232:233];
% set(gca,'ycolor','k','ytick',1:length(ex_trial_num),'yticklabel',flip(ex_trial_num),'ylim',[0.5 length(ex_trial_num)+0.5]); 
% ylabel('trial #');
% set(gcf,'position',[411 100 560 500]);
% c=0;
% for i=ex_trial_num
%     c=c+1;
%     [this_time, this_tail, ~, bout_starts, bout_ends] = get_example_trial_data(i,time_be,tail,grspeed,lt_lag_TL(fish_num));
%     this_tail=this_tail/5+length(ex_trial_num)+1-c;
%     fill_bout (gca,bout_starts(1),bout_ends(1),-0.5+length(ex_trial_num)+1-c,0.5+length(ex_trial_num)+1-c,0.3);
%     plot(this_time,this_tail,'color',col_motor,'linewidth',1);
% end
% show_lt_trials_lines(-7.4,-7.4,6.5,8.5,col_pre);
% show_lt_trials_lines(-7.4,-7.4,4.5,6.5,col_adapt_start);
% show_lt_trials_lines(-7.4,-7.4,2.5,4.5,col_adapt_end);
% show_lt_trials_lines(-7.4,-7.4,0.5,2.5,col_post);
% show_trial_start_and_end;
% dm_fix_fig_fonts;
% clear i fish_id c this_time this_tail bout_starts bout_ends time_be tail grspeed ex_trial_num;
% 
% % Fig. 4c: All fish: first bout duration during all trials
% lt_norm_TL=load([pathname_behavior 'long_term_adaptation_experiment/normal_reafference_control/WT_TL_group/pooled_data.mat'],'behavior');
% lt_norm_TL=lt_norm_TL.behavior;
% plot_long_term_adaptation('Fig. 4c',{lt_norm_TL, lt_lag_TL},'bout_duration','first',[col_norm_reaf; col_lag_trained],{'-', '-'},...
%     col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post,[0.25 0.65]);
% dm_fix_fig_fonts;
% 
% % Fig. 4d-f: Quantificication of acute reaction, back-to-baseline and after-effect
% p_value_lt_TL = nan(1,3);
% fig_names = {'Fig. 4d', 'Fig. 4e', 'Fig. 4f'};
% temp_labels = {'(first 10 adaptation trials - pre) [s]','(last 10 adaptation trials - first 10 adaptation trials) [s]','(post - pre) [s]'};
% temp_trial_num = {21:30, 221:230, 231:240};
% temp_trial_num0 = {11:20, 21:30, 11:20};
% for i=1:3
%     p_temp=make_lt_fbd_panel(fig_names{i},{'normalized first bout duration'; temp_labels{i}},{lt_norm_TL, lt_lag_TL},'bout_duration','first',temp_trial_num{i},[0 0 0],fish_num,temp_trial_num0{i});
%     dm_fix_fig_fonts;
%     p_value_lt_TL(i) = p_temp(1,2);
% end
% clear temp_trial_num0 p_temp fish_num fig_names temp_labels temp_trial_num;
% 
% % Fig. 4g: Bout power
% h=create_bout_power_axes('Fig. 4g',2);
% set(gcf,'position',[206 10 1034 420]);
% lt_norm_TL_bp=extract_bout_power(lt_norm_TL);
% lt_lag_TL_bp=extract_bout_power(lt_lag_TL);
% time_power=-0.1:dt:1-dt;
% h_bp_norm=plot_bout_power(h(1),time_power,lt_norm_TL_bp,[2 3 24],[col_pre; col_adapt_start; col_post]);
% h_bp_lag=plot_bout_power(h(2),time_power,lt_lag_TL_bp,[2 3 24],[col_pre; col_adapt_start; col_post]);
% y_lim=[0 4];
% set(h,'ylim',y_lim,'ytick',0:4);
% y=[y_lim(2)-diff(y_lim)/10 y_lim(2)-diff(y_lim)/12];
% h_bp=cat(3,h_bp_norm,h_bp_lag);
% cols=[col_adapt_start; col_post];
% for k=1:2 % norm or lag
%     for j=1:2 % acute or long-term
%         for i=1:220 % time
%             if h_bp(j,i,k)
%                 line(h(k),[time_power(i)-dt/2 time_power(i)+dt/2],[y(j) y(j)],'color',cols(j,:),'linewidth',3);
%             end
%         end
%     end
% end
% ballistic_end=time_power(find(h_bp_lag(1,:),1));
% ballistic_change_start=time_power(find(h_bp_lag(2,:),1));
% ballistic_change_end=time_power(find(h_bp_lag(2,time_power<0.15),1,'last'));
% dm_fix_fig_fonts;
% clear lt_norm_TL lt_lag_TL y_lim y h_bp h_bp_norm h_bp_lag cols k j i h;
% 
% % Fig. 4h: Quantification of acute reaction in ballistic power
% prepare_quantif_axes2('Fig. 4h');
% data = prepare_data_for_bout_power_quantif({lt_norm_TL_bp, lt_lag_TL_bp},3,0.1/dt:(ballistic_end+0.1)/dt);
% p_value_lt_TL_acute_ballistic_power = plot_quantif(data,[0 0 0],{'normalized first bout ballistic power'; '(first 10 adaptation trials - pre) [au]'});
% set(gca,'ylim',[-1 1]);
% dm_fix_fig_fonts;
% clear data;
% 
% % Fig. 4i: Quantification of acute reaction in reactive power
% prepare_quantif_axes2('Fig. 4i');
% data = prepare_data_for_bout_power_quantif({lt_norm_TL_bp, lt_lag_TL_bp},3,(ballistic_end+0.1)/dt+1:220);
% p_value_lt_TL_acute_reactive_power = plot_quantif(data,[0 0 0],{'normalized first bout reactive power'; '(first 10 adaptation trials - pre) [au]'});
% dm_fix_fig_fonts;
% clear data;
% 
% % Fig. 4j: Quantification of long-term adaptation of ballistic power
% prepare_quantif_axes2('Fig. 4j');
% data = prepare_data_for_bout_power_quantif({lt_norm_TL_bp, lt_lag_TL_bp},24,(0.1+ballistic_change_start)/dt+1:(0.1+ballistic_change_end)/dt+1);
% p_value_lt_TL_aftereffect_ballistic_power = plot_quantif(data,[0 0 0],{'normalized first bout ballistic power'; '(post - pre) [au]'});
% set(gca,'ylim',[-2 6]);
% dm_fix_fig_fonts;
% clear data lt_norm_TL_bp lt_lag_TL_bp;


%% Fig. 5: Long-term adaptation, but not acute reaction, is impaired after PC ablation
% % depends on Figure 4
% % Fig. 5a: Experimental flow
% % created in illustrator
% 
% % Fig. 5bi: All treatment control fish: first bout duration during all trials
% lt_norm_PC_neg=load([pathname_behavior 'long_term_adaptation_experiment/normal_reafference_control/treatment_control_group/pooled_data.mat'],'behavior');
% lt_norm_PC_neg=lt_norm_PC_neg.behavior;
% lt_lag_PC_neg=load([pathname_behavior 'long_term_adaptation_experiment/lag_trained/treatment_control_group/pooled_data.mat'],'behavior');
% lt_lag_PC_neg=lt_lag_PC_neg.behavior;
% plot_long_term_adaptation('Fig. 5bi',{lt_norm_PC_neg, lt_lag_PC_neg},'bout_duration','first',[col_norm_reaf; col_lag_trained],{'-', '-'},...
%     col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post,[0.25 0.6]);
% dm_fix_fig_fonts;
% 
% % Fig. 5bii: All PC-ablated fish: first bout duration during all trials
% lt_norm_PC_pos=load([pathname_behavior 'long_term_adaptation_experiment/normal_reafference_control/PC_ablated_group/pooled_data.mat'],'behavior');
% lt_norm_PC_pos=lt_norm_PC_pos.behavior;
% lt_lag_PC_pos=load([pathname_behavior 'long_term_adaptation_experiment/lag_trained/PC_ablated_group/pooled_data.mat'],'behavior');
% lt_lag_PC_pos=lt_lag_PC_pos.behavior;
% plot_long_term_adaptation('Fig. 5bii',{lt_norm_PC_pos, lt_lag_PC_pos},'bout_duration','first',[col_norm_reaf; col_lag_trained],{'-', '-'},...
%     col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post,[0.25 0.6]);
% dm_fix_fig_fonts;
% 
% % Fig. 5c: Quantification of acute reaction
% p_value_lt_PC_acute=make_lt_fbd_panel_PC('Fig. 5c',{'normalized first bout duration'; '(first 10 adaptation trials - pre) [s]'},{lt_norm_PC_neg, lt_lag_PC_neg, lt_norm_PC_pos, lt_lag_PC_pos},'bout_duration','first',21:30,[col_norm_reaf; col_lag_trained],'both');
% dm_fix_fig_fonts;
% 
% % Fig. 5d: Quantification of back to baseline effect
% p_value_lt_PC_back_to_baseline=make_lt_fbd_panel_PC('Fig. 5d',{'normalized first bout duration'; '(last 10 adaptation trials - first 10 adaptation trials) [s]'},{lt_norm_PC_neg, lt_lag_PC_neg, lt_norm_PC_pos, lt_lag_PC_pos},'bout_duration','first',221:230,[col_norm_reaf; col_lag_trained],'both',21:30);
% dm_fix_fig_fonts;
% 
% % Fig. 5e: Quantification of after-effect
% p_value_lt_PC_after_effect=make_lt_fbd_panel_PC('Fig. 5e',{'normalized first bout duration'; '(post - pre) [s]'},{lt_norm_PC_neg, lt_lag_PC_neg, lt_norm_PC_pos, lt_lag_PC_pos},'bout_duration','first',231:240,[col_norm_reaf; col_lag_trained],'both');
% dm_fix_fig_fonts;
% 
% % Fig. 5f: Bout power
% h=create_bout_power_axes('Fig. 5f',2);
% set(gcf,'position',[206 10 1034 420]);
% y_lim=[0 4];
% set(h,'ylim',y_lim);
% lt_norm_PC_neg_bp=extract_bout_power(lt_norm_PC_neg);
% lt_lag_PC_neg_bp=extract_bout_power(lt_lag_PC_neg);
% lt_PC_neg_bp=cat(3,lt_norm_PC_neg_bp,lt_lag_PC_neg_bp);
% lt_norm_PC_pos_bp=extract_bout_power(lt_norm_PC_pos);
% lt_lag_PC_pos_bp=extract_bout_power(lt_lag_PC_pos);
% lt_PC_pos_bp=cat(3,lt_norm_PC_pos_bp,lt_lag_PC_pos_bp);
% clear lt_norm_PC_neg_bp lt_lag_PC_neg_bp lt_norm_PC_pos_bp lt_lag_PC_pos_bp lt_lag_PC_neg lt_lag_PC_pos lt_norm_PC_neg lt_norm_PC_pos;
% h_lt_bp_PC_neg=plot_bout_power(h(1),time_power,lt_PC_neg_bp,[2 24],[col_pre; col_post]);
% h_lt_bp_PC_pos=plot_bout_power(h(2),time_power,lt_PC_pos_bp,[2 24],[col_pre; col_post]);
% y1=y_lim(2)-diff(y_lim)/10;
% y2=y_lim(2)-diff(y_lim)/12;
% for i=1:220
%     if h_lt_bp_PC_neg(i)
%         line(h(1),[time_power(i)-dt/2 time_power(i)+dt/2],[y1 y1],'color',col_post,'linewidth',3);
%     end
%     if h_lt_bp_PC_pos(i)
%         line(h(2),[time_power(i)-dt/2 time_power(i)+dt/2],[y2 y2],'color',col_post,'linewidth',3);
%     end
% end
% dm_fix_fig_fonts;
% clear h_lt_bp_PC_neg h_lt_bp_PC_pos y_lim y1 y2 i h
% 
% % Fig. 5g: Quantification of long-term adaptation of ballistic power
% prepare_quantif_axes2('Fig. 5g');
% set(gca,'xticklabel',{'treatment control','PC-ablated'});
% data = prepare_data_for_bout_power_quantif({lt_PC_neg_bp, lt_PC_pos_bp},24,(0.1+ballistic_change_start)/dt+1:(0.1+ballistic_change_end)/dt+1);
% p_value_lt_PC_aftereffect_ballistic_power = plot_quantif(data,[0 0 0],{'normalized first bout ballistic power'; '(post - pre) [au]'},[],'both');
% set(gca,'ylim',[-2 6]);
% dm_fix_fig_fonts;
% clear data lt_PC_neg_bp lt_PC_pos_bp;


%% Fig. 6: Activity of a subpopulation of PCs can represent the output of an internal model
% % extract the data
% [all_fish, n_fish, n_control, n_good_lag, n_bad_lag,...
%     tf_control_fish, tf_good_lag_fish, tf_bad_lag_fish,...
%     tf_control_rois, tf_good_lag_rois, tf_bad_lag_rois, fish_id_rois,...
%     Scores_be, Scores_im, ~,...
%     Crit_im, Crit_im_signif,...
%     Traces_bout_trig, ~, time_trig,...
%     ROI_coord,sz_PC_ref,rez_PC_ref] = extract_PC_imaging_data(pathname_PC_imaging,pathname_reference_brains);
% fish_id_ex = '200303_f0';
% fish_num_ex = find(strcmp(all_fish,fish_id_ex));
% tf_mot_rois = Crit_im_signif(:,1)==1;
% tf_exp_rois = Crit_im_signif(:,1)==0 & Crit_im_signif(:,2)==-1 & Crit_im_signif(:,3)==0 & Crit_im_signif(:,4)==1;
% Data_temp_exp_rois = nan(n_fish,120);
% trig_data = compute_mean_trig_data(Traces_bout_trig);
% Trig_data_exp_rois = nan(n_fish,size(trig_data,2),size(trig_data,3));
% for i=1:n_fish
%     Data_temp_exp_rois(i,:)=nanmean(Scores_im(fish_id_rois==i & tf_exp_rois,:),1);
%     Trig_data_exp_rois(i,:,:) = nanmean(trig_data(fish_id_rois==i & tf_exp_rois,:,:),1);
% end
% exp_roi_ids=find(tf_exp_rois & fish_id_rois==fish_num_ex);
% mot_roi_ids=find(tf_mot_rois & fish_id_rois==fish_num_ex);
% ex_roi_ids = [exp_roi_ids(16) mot_roi_ids(13)];
% clear exp_roi_ids mot_roi_ids fish_num_ex all_fish
% 
% % Fig. 6a: Schematics of the light-sheet setup
% % created in illustrator
% % inset with location of example ROIs
% tf_roi1=false(sz_PC_ref);
% tf_roi1(ROI_coord{ex_roi_ids(1)})=true;
% tf_roi1 = max(tf_roi1,[],3);
% tf_roi2=false(sz_PC_ref);
% tf_roi2(ROI_coord{ex_roi_ids(2)})=true;
% tf_roi2 = max(tf_roi2,[],3);
% r=255*ones(sz_PC_ref(1:2),'uint8');
% g=r; b=r;
% g(tf_roi1)=0;
% b(tf_roi1)=0;
% r(tf_roi2)=0;
% g(tf_roi2)=0;
% A = cat(3,r,g,b);
% figure('name','Fig. 6a location of example ROIs');
% imshow(A);
% dm_fix_fig_fonts;
% clear A tf_roi1 tf_roi2 r g b z_mm z1 z2 anat;
% 
% % Fig. 6b: Experimental protocol
% % created in illustrator
% 
% % Fig. 6c: Behavioral plots, divided into three groups
% p=plot_long_term_adaptation_PC_im('Fig. 6c',{Scores_be(tf_control_fish,:), Scores_be(tf_bad_lag_fish,:), Scores_be(tf_good_lag_fish,:)},[col_norm_reaf; col_lag_trained_bad; col_lag_trained],{'-','-','-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0.2 0.8]);
% legend(p(1:3),{['control, N = ' num2str(sum(tf_control_fish))],['lag-trained non-adapting, N = ' num2str(sum(tf_bad_lag_fish))],['lag-trained adapting, N = ' num2str(sum(tf_good_lag_fish))]},'box','off','color','none');
% dm_fix_fig_fonts;
% clear p Scores_be;
% 
% % Fig. 6d: Computing criteria on two example ROIs
% % extract the data
% trials_to_show=[19 22 70 75 113];
% load([pathname_PC_imaging 'behavior\' fish_id_ex '_behavior.mat'],'time_be','grspeed','tail','bouts','meta');
% temp_struct.bouts=bouts;
% load([pathname_PC_imaging 'traces\' fish_id_ex '_traces.mat'],'traces','time_offsets');
% time_im=(1/meta.F:1/meta.F:size(traces,2)/meta.F)+1/meta.F;
% time_im_roi1=time_im + time_offsets(ex_roi_ids(1));
% time_im_roi2=time_im + time_offsets(ex_roi_ids(2));
% clear time_im;
% % Fig. 6di: Two raw ROIs with tail and grating speed
% create_single_trial_axes('Fig. 6di'); 
% set(gca,'ycolor','k','ytick',1:5,'yticklabel',flip(trials_to_show),'ylim',[0 length(trials_to_show)+1]); 
% ylabel('trial #');
% set(gcf,'position',[411 10 560 853]);
% c=0;
% for i=trials_to_show
%     c=c+1;
%     [this_time, this_tail, this_grspeed, bout_starts, bout_ends] = get_example_trial_data(i,time_be,tail,grspeed,temp_struct);
%     this_tail=this_tail/20+length(trials_to_show)+1-c;
%     this_grspeed=this_grspeed/120+length(trials_to_show)+1.1-c;
%     fill_bout (gca,bout_starts(1),bout_ends(1),-0.5+length(trials_to_show)+1-c,0.5+length(trials_to_show)+1-c,0.3);
%     plot(this_time,this_tail,'color',col_motor,'linewidth',1);
%     plot(this_time,this_grspeed,'color',col_sensory,'linewidth',0.5);
%     
%     tf=time_im_roi1>(i-1)*30 & time_im_roi1<=i*30;
%     this_trace = traces(ex_roi_ids(1),tf);
%     this_trace = this_trace/10+length(trials_to_show)+0.9-c;
%     this_time = time_im_roi1(tf);
%     this_time=this_time-this_time(1)-7.5-1/meta.F/2;
%     plot(this_time,this_trace,'color','r','linewidth',1);
%     tf=time_im_roi2>(i-1)*30 & time_im_roi2<=i*30;
%     this_trace = traces(ex_roi_ids(2),tf);
%     this_trace = this_trace/10+length(trials_to_show)+0.7-c;
%     this_time = time_im_roi2(tf);
%     this_time=this_time-this_time(1)-7.5-1/meta.F/2;
%     plot(this_time,this_trace,'color','b','linewidth',1);
% end
% show_lt_trials_lines(-7.4,-7.4,4.5,5.5,col_pre);
% show_lt_trials_lines(-7.4,-7.4,3.5,4.5,col_adapt_start);
% show_lt_trials_lines(-7.4,-7.4,2.5,3.5,col_adapt_end);
% show_lt_trials_lines(-7.4,-7.4,1.5,2.5,col_post_start);
% show_lt_trials_lines(-7.4,-7.4,0.5,1.5,col_post_end);
% line([-6 -6],[0 10]/120+length(trials_to_show)+1.1-c,'color','k');
% line([-5 -5],[0 1]/10+length(trials_to_show)+0.7-c,'color','k');
% show_trial_start_and_end;
% dm_fix_fig_fonts;
% clear this_trace tf this_time this_tail this_gr bout_starts bout_ends time_im_roi1 time_im_roi2 temp_struct;
% % Fig. 6dii: Triggered activity of these two ROIs
% figure('name','Fig. 6dii'); axes; hold on;
% set(gcf,'position',[411 10 300 853]);
% ylabel('trial #');
% xlabel('time relative to first bout onset [s]');
% set(gca,'xlim',[-0.8 1.2],'xtick',-0.8:0.4:1.2,'ycolor','k','ytick',1:5,'yticklabel',flip(trials_to_show),'ylim',[0 length(trials_to_show)+1]);
% c=0;
% for i=trials_to_show
%     c=c+1;
%     this_trace = Traces_bout_trig(ex_roi_ids(1),:,i);
%     this_trace = this_trace/4+length(trials_to_show)+0.6-c;
%     plot(time_trig,this_trace,'color','r','linewidth',1);
%     this_trace = Traces_bout_trig(ex_roi_ids(2),:,i);
%     this_trace = this_trace/4+length(trials_to_show)+0.6-c;
%     plot(time_trig,this_trace,'color','b','linewidth',1);
%     line([-0.8 1.2],[0 0]+length(trials_to_show)+0.6-c,'color','k','linestyle',':');
% end
% show_lt_trials_lines(-0.95,-0.95,4.5,5.5,col_pre);
% show_lt_trials_lines(-0.95,-0.95,3.5,4.5,col_adapt_start);
% show_lt_trials_lines(-0.95,-0.95,2.5,3.5,col_adapt_end);
% show_lt_trials_lines(-0.95,-0.95,1.5,2.5,col_post_start);
% show_lt_trials_lines(-0.95,-0.95,0.5,1.5,col_post_end);
% line([-0.5 -0.5],[0 1]/4+length(trials_to_show)+0.6-c,'color','k');
% show_trial_start_and_end;
% dm_fix_fig_fonts;
% clear this_trace trials_to_show Traces_bout_trig;
% % Fig. 6diii: Scores of these two ROIs
% cols=['r';'b'];
% plot_long_term_adaptation_PC_im('Fig. 6diii',{Scores_im(ex_roi_ids(1),:), Scores_im(ex_roi_ids(2),:)},cols,{'-','-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0 4]);
% ylabel('bout-triggered response [sd]');
% set(gca,'ytick',0:4);
% temp_filt=[Scores_im(ex_roi_ids(1),:); Scores_im(ex_roi_ids(2),:)];
% for j=1:2
%     for i=5:length(temp_filt)-4
%         temp_filt(j,i)=nanmean(Scores_im(ex_roi_ids(j),i-4:i+4));
%     end
%     plot(1:120,temp_filt(j,:),'color',cols(j),'linewidth',2);
% end
% dm_fix_fig_fonts;
% clear j temp_filt cols;
% % Fig. 6div: Criteria of these two ROIs
% figure('name','Fig. 6div');
% axes('xlim',[0 4]+0.5,'xtick',1:4,'ylim',[0 2]+0.5,'ytick',1:2,'yticklabel',{'ROI 1','ROI 2'},'ydir','reverse'); hold on;
% xlabel('criterion #');
% dm_imagesc(1:4,[Crit_im(ex_roi_ids(1),:); Crit_im(ex_roi_ids(2),:)]);
% line(xlim,[1.5 1.5],'color','k','linewidth',2);
% ylabel(colorbar,{'difference in responses [sd]'});
% dm_fix_fig_fonts;
% clear Crit_im;
% % Fig. 6dv: Significant criteria of these two ROIs
% figure('name','Fig. 6dv');
% axes('xlim',[0 4]+0.5,'xtick',1:4,'ylim',[0 2]+0.5,'ytick',1:2,'yticklabel',{'ROI 1','ROI 2'},'ydir','reverse'); hold on;
% xlabel('criterion #');
% dm_imagesc(1:4,[Crit_im_signif(ex_roi_ids(1),:); Crit_im_signif(ex_roi_ids(2),:)]);
% line(xlim,[1.5 1.5],'color','k','linewidth',2);
% cl1=colormap;
% cl1=[cl1(1,:); cl1(32,:); cl1(end,:)];
% colormap(gca,cl1);
% colorbar('ytick',[-0.6 0 0.6],'yticklabel',{'decrease','no change','increase'});
% dm_fix_fig_fonts;
% clear cl1;
% % Fig. 6dvi: triggered responses of ROI1
% plot_trigaver_PC('Fig. 6dvi',time_trig,trig_data(ex_roi_ids(1),:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2]);
% dm_fix_fig_fonts;
% % Fig. 6dvii: triggered responses of ROI2
% plot_trigaver_PC('Fig. 6dvii',time_trig,trig_data(ex_roi_ids(2),:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 4]);
% dm_fix_fig_fonts;
% 
% % Fig. 6e: Clustering of ROIs
% % extract data
% clust_code='+-0';
% signif_code=[1 -1 0];
% clust_data=cell(81,7);
% % 1 - cluster code
% % 2 - ids of control ROIs
% % 3 - ids of bad lag ROIs
% % 4 - ids of good lag ROIs
% % 5 - fractions in control fish
% % 6 - fractions in bad lag fish
% % 7 - fractions in good lag fish
% data_fraction=nan(n_fish,81);
% p_values_PC_fraction = nan(81,1);
% c=0;
% for c1=1:3
%     for c2=1:3
%         for c3=1:3
%             for c4=1:3
%                 c=c+1;
%                 clust_data{c,1}=[clust_code(c1) clust_code(c2) clust_code(c3) clust_code(c4)];
%                 temp = ...
%                     Crit_im_signif(:,1)==signif_code(c1) & ...
%                     Crit_im_signif(:,2)==signif_code(c2) & ...
%                     Crit_im_signif(:,3)==signif_code(c3) & ...
%                     Crit_im_signif(:,4)==signif_code(c4);
%                 clust_data{c,2} = find(temp & tf_control_rois);
%                 clust_data{c,3} = find(temp & tf_bad_lag_rois);
%                 clust_data{c,4} = find(temp & tf_good_lag_rois);
%                 for f=1:n_fish
%                     temp2 = temp & fish_id_rois==f;
%                     if tf_control_fish(f)
%                         id = 5;
%                     elseif tf_bad_lag_fish(f)
%                         id = 6;
%                     elseif tf_good_lag_fish(f)
%                         id = 7;
%                     end
%                     clust_data{c,id} = [clust_data{c,id} sum(temp2)/sum(fish_id_rois==f)*100]; 
%                 end
%                 data_fraction(:,c)=[clust_data{c,5}'; clust_data{c,6}'; clust_data{c,7}'];
%                 p_values_PC_fraction(c,1) = kruskalwallis(data_fraction(:,c),[ones(n_control,1);2*ones(n_bad_lag,1);3*ones(n_good_lag,1)],'off');
%             end
%         end
%     end
% end
% Crit_im_signif_sorted=[];
% Scores_im_sorted=[];
% interesting_cluster_id = find(strcmp(clust_data(:,1),'0-0+'));
% line_start = [];
% line_end = [];
% cc=0;
% n=zeros(3,1);
% for g=2:4
%     cc=cc+1;
%     for c=1:81
%         n(cc)=n(cc)+length(clust_data{c,g});
%         if c==interesting_cluster_id
%             line_start = [line_start size(Crit_im_signif_sorted,1)]; %#ok<AGROW>
%         end
%         Crit_im_signif_sorted=[Crit_im_signif_sorted; Crit_im_signif(clust_data{c,g},:)]; %#ok<AGROW>
%         Scores_im_sorted=[Scores_im_sorted; Scores_im(clust_data{c,g},:)]; %#ok<AGROW>
%         if c==interesting_cluster_id
%             line_end = [line_end size(Crit_im_signif_sorted,1)]; %#ok<AGROW>
%         end
%     end
% end
% n=[1; n];
% clear g cc interesting_cluster_id id temp2 temp c1 c2 c3 c4 signif_code clust_code Crit_im_signif;
% 
% % Fig. 6ei: Clustering of all ROIs within 3 groups of fish (plot crit signif)
% figure('name','Fig. 6ei');
% set(gcf,'position',[680 20 159 420]);
% for i=1:3
%     subplot(3,1,i,'xlim',[0 4]+0.5,'xtick',1:4,'xticklabel',[],'ytick',[],'ydir','reverse'); hold on;
%     ylabel('ROIs');
%     dm_imagesc(1:4,Crit_im_signif_sorted(sum(n(1:i)):sum(n(2:i+1)),:));
% end
% xlabel('criterion #');
% line([0 4]+0.5,[line_start(3) line_start(3)]-(n(2)+n(3))+0.5,'color','k','linewidth',1);
% line([0 4]+0.5,[line_end(3) line_end(3)]-(n(2)+n(3))+0.5,'color','k','linewidth',1);
% dm_fix_fig_fonts;
% % Fig. 6eii: Fraction of ROIs in individual fish
% c_lim=[0 10];
% useless_clust = max([mean(data_fraction(1:n_control,:),1); mean(data_fraction(n_control+1:n_control+n_bad_lag,:),1); mean(data_fraction(n_control+n_bad_lag+1:n_fish,:),1)],[],1)<2;
% data_fraction(:,useless_clust)=[];
% clust_data(useless_clust,:)=[];
% p_values_PC_fraction(useless_clust)=[];
% good_cluster_id=find(p_values_PC_fraction<0.05);
% figure('name','Fig. 6eii');
% set(gcf,'position',[990 20 290 420]);
% n_temp = [1 n_control n_bad_lag n_good_lag];
% for i=1:3
%     subplot(3,1,i,'xlim',[0 size(data_fraction,2)]+0.5,'ylim',[0 n_temp(i)]+0.5,'ydir','reverse','xtick',1:size(data_fraction,2),'xticklabel',[],'ytick',[]); hold on;
%     ylabel('fish');
%     dm_imagesc(1:size(data_fraction,2),data_fraction(sum(n_temp(1:i)):sum(n_temp(2:i+1)),:));
%     line([good_cluster_id good_cluster_id]-0.5,ylim,'color','k');
%     line([good_cluster_id good_cluster_id]+0.5,ylim,'color','k');
%     colormap(flip(bone));
%     set(gca,'clim',c_lim);
%     colorbar;
% end
% set(gca,'xticklabel',clust_data(:,1));
% xtickangle(45);
% dm_fix_fig_fonts;
% clear good_cluster_id useless_clust line_start line_end Crit_im_signif_sorted Scores_im_sorted data_fraction clust_data;
% 
% % Fig. 6f: Responses of 0-0+ ROIs
% % Fig. 6fi top: Scores of individual ROIs
% c_lim=[-1 3];
% data_temp = Scores_im(tf_exp_rois,:);
% data_temp(isnan(data_temp))=nanmean(data_temp(:));
% data_temp=data_temp(randperm(size(data_temp,1)),:);
% figure('name','Figure 6fi top');
% imagesc_scores(1:2,data_temp,c_lim,col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end);
% xlabel('trial #');
% set(gca,'xticklabelmode','auto');
% dm_fix_fig_fonts;
% % Fig. 6fi bottom: Scores averaged across fish
% plot_long_term_adaptation_PC_im('Fig. 6fi bottom',{Data_temp_exp_rois(tf_good_lag_fish,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0 1.5]);
% ylabel('bout-triggered response [sd]');
% dm_fix_fig_fonts;
% clear data_temp c_lim Data_temp_exp_rois Scores_im tf_control tf_bad_lag_fish;
% % Fig. 6fii: Triggered avearges, mean across good lag fish
% plot_trigaver_PC('Fig. 6fii',time_trig,Trig_data_exp_rois(tf_good_lag_fish,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2]);
% dm_fix_fig_fonts;
% clear ex_roi_ids Trig_data_exp_rois trig_data tf_good_lag_fish time_trig;
% 
% % Fig. 6g: Spatial organisation of 0-0+ ROIs
% % control fish
% A=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_control,tf_control_rois,tf_exp_rois,fish_id_rois,ROI_coord,1);
% show_PC_map('Fig. 6g control',A,rez_PC_ref,50);
% dm_fix_fig_fonts;
% % bad lag fish (nonadapting)
% A=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_bad_lag,tf_bad_lag_rois,tf_exp_rois,fish_id_rois,ROI_coord,1);
% show_PC_map('Fig. 6g bad lag',A,rez_PC_ref,50);
% dm_fix_fig_fonts;
% % good lag fish (adapting)
% A=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_good_lag,tf_good_lag_rois,tf_exp_rois,fish_id_rois,ROI_coord,1);
% show_PC_map('Fig. 6g good lag',A,rez_PC_ref,50);
% dm_fix_fig_fonts;


%% Fig. 7: A cerebellar internal model calibrates a feedback controller involved in sensorimotor control
% load the data
[sz, rez] = read_nrrd_metadata([pathname_reference_brains 'PortuguesLab_wholebrain_ref_gamma_masked.nrrd']);
all_fish = dm_dir([pathname_whole_brain_imaging_lta 'processed_data\*_f*_processed_data.mat']);
all_fish = strrep(all_fish,'_processed_data.mat','');
n_fish = length(all_fish);
ROI_coord_all = cell(n_fish,1);
n_ROIs_all = zeros(n_fish,1);
sensmot_clust_all = cell(n_fish,1);
traces_gr_trig_all = cell(n_fish,1);
traces_bout_trig_all = cell(n_fish,1);
time_constants_all = cell(n_fish,1);
time_constants_trials_all = cell(n_fish,1);
group = nan(n_fish,1);
load([pathname_whole_brain_imaging_lta 'processed_data\' all_fish{1} '_processed_data.mat'],'time_trig')
progressbar('Loading data for Fig. 7...');
for f = 1:n_fish
    fish_id = all_fish{f};
    load([pathname_whole_brain_imaging_lta 'ROIs\' fish_id '_ROIs.mat'],'ROI_coord');
    ROI_coord_all{f} = ROI_coord;
    n_ROIs_all(f) = length(ROI_coord);
    load([pathname_whole_brain_imaging_lta 'processed_data\' fish_id '_processed_data.mat'],'traces_gr_trig_mean','traces_bout_trig_mean','sensmot_clust','time_constants');
    traces_gr_trig_all{f} = traces_gr_trig_mean;
    traces_bout_trig_all{f} = traces_bout_trig_mean;
    sensmot_clust_all{f} = sensmot_clust;
    time_constants_all{f} = nanmean(time_constants,2);
    time_constants_trials_all{f} = time_constants;
    load([pathname_whole_brain_imaging_lta 'behavior\' fish_id '_behavior.mat'],'meta')
    switch meta.group
        case 'control'
            group(f) = 1;
        case 'lag-trained non-adapting'
            group(f) = 2;
        case 'lag-trained adapting'
            group(f) = 3;
    end      
    progressbar(f/n_fish);
end
my_colormap=hot(64);
my_colormap=my_colormap(round(linspace(1,64,sum(group==3)+1)),:);

% Fig. 7a: change in taus
tau_change1 = [];
tau_change2 = [];
tau_change3 = [];
x = -10:0.01:10;
for f=1:n_fish 
    tau_change = time_constants_trials_all{f}(:,7) - time_constants_trials_all{f}(:,3);
    pd = fitdist(tau_change,'kernel','width',0.1);
    tau_change = pdf(pd,x);
    tau_change = tau_change/sum(tau_change);
     switch group(f)
        case 1
            tau_change1 = [tau_change1;  tau_change];
        case 2
            tau_change2 = [tau_change2;  tau_change];
        case 3
            tau_change3 = [tau_change3;  tau_change];
    end
end
figure('name','Fig. 7a');
axes; hold on;
fill([x flip(x)], [nanmean(tau_change1,1) - nanstd(tau_change1,[],1)/sqrt(size(tau_change1,1)) flip(nanmean(tau_change1,1) + nanstd(tau_change1,[],1)/sqrt(size(tau_change1,1)))],col_norm_reaf,'edgecolor','none','facealpha',0.3);
plot(x,nanmean(tau_change1,1),'color',col_norm_reaf,'linewidth',2,'linestyle',':');
fill([x flip(x)], [nanmean(tau_change2,1) - nanstd(tau_change2,[],1)/sqrt(size(tau_change2,1)) flip(nanmean(tau_change2,1) + nanstd(tau_change2,[],1)/sqrt(size(tau_change2,1)))],col_lag_trained_bad,'edgecolor','none','facealpha',0.3);
plot(x,nanmean(tau_change2,1),'color',col_lag_trained_bad,'linewidth',2,'linestyle',':');
fill([x flip(x)], [nanmean(tau_change3,1) - nanstd(tau_change3,[],1)/sqrt(size(tau_change3,1)) flip(nanmean(tau_change3,1) + nanstd(tau_change3,[],1)/sqrt(size(tau_change3,1)))],col_lag_trained,'edgecolor','none','facealpha',0.3);
plot(x,nanmean(tau_change3,1),'color',col_lag_trained,'linewidth',2,'linestyle',':');
xlabel('change in tau (criterion 2) [s]');
ylabel('probability');
set(gca,'xlim',[-4 4]);
dm_fix_fig_fonts;

% maps of cells with crit 2<-0.4
for gg=1:3
    % Fig. 7b: map of all sensors
    % compute stacks with fraction of ROIs
    A_sens = zeros(sz,'uint8');
    for f = 1:n_fish
        if group(f)==gg
            tau_change = time_constants_trials_all{f};
            crit2 = tau_change(:,7) - tau_change(:,3);
            for i = 1:n_ROIs_all(f)
                coord = ROI_coord_all{f}{i};
                if crit2(i)<-0.4
                    A_sens(coord) = A_sens(coord) + 1;
                end
            end
        end
    end
    % compute top and side views
    [top_view_sens, side_view_sens]=three_orthogonal_views(A_sens,rez(1));
    top_view_sens=1-top_view_sens/sum(group==gg);
    side_view_sens=1-round(side_view_sens)/sum(group==gg);
    % plot top and side views of the fractions of sensory ROIs
    figure('name',['Fig. 7b top g' num2str(gg)]);
    axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 sz(2)],'ylim',[1 sz(1)]); hold on;
    pbaspect([1 sz(1)/sz(2) 1]);
    imagesc(top_view_sens);
    colormap(my_colormap);
    figure('name',['Fig. 7b side g' num2str(gg)]);
    axes('xcolor','none','ycolor','none','ydir','reverse','clim',[0 1],'xlim',[1 round(sz(3)/rez(2))],'ylim',[1 sz(1)]); hold on;
    pbaspect([1 sz(1)/round(sz(3)/rez(2)) 1]);
    imagesc(side_view_sens);
    colormap(my_colormap);
end
figure('name','Fig. 7b colormap');
axes('ycolor','none','xcolor','none');
colormap(flip(hot(64),1));
clrbr=colorbar;
clrbr.Ticks=0:0.2:1;
ylabel(clrbr,{'percentage of fish';'with activity'});
dm_fix_fig_fonts;


%% Extended Data Fig. 1: Behavior of the feedback control model of acute reaction
% figure('name','Extended Data Fig. 1');
% set(gcf,'position',[188 20 1547 365]);
% par_trial=[0.161, 0.15, 2.8, 2.5, 0.8, 0.9, 0.5, 0.6];
% dt=0.005;
% par_trial(3)=dt/par_trial(3);
% par_trial(8)=dt/par_trial(8);
% reaf_trial=ones(7,5).*[1 0 0 0 0]; % all bouts with normal reafference
% [swim,grspeed,brain_state]=model_v3_real_trial (par_trial, dt, reaf_trial, 3);
% int_gr=zeros(1,length(swim));
% int_gr(3/dt:18/dt)=10;
% axes('ycolor','none','xlim',[1 length(swim)],'xtick',[1 3/dt:5/dt:(3+15)/dt length(swim)],'xticklabel',[-3 0 5 10 15 15+3]); hold on
% 
% plot(int_gr/10*2,'color',col_sensory);
% plot(swim*2-4,'color',col_motor);
% plot(brain_state(1,:)/10*2 - 8,'color',col_activity)
% plot(brain_state(2,:)/20*2 - 12,'color',col_activity)
% for i=3:5
%     plot(brain_state(i,:)*5 - 4*(i+1),'color',col_activity)
% end
% plot(grspeed/10*2-4*(i+2),'color',col_sensory);
% drawnow;
% bs=find(diff([0 swim])==1);
% be=find(diff([swim 0])==-1);
% for i=1:length(bs)
%     fill([bs(i) be(i) be(i) bs(i)],[-1 -1 1 1],'k');
% end
% dm_fix_fig_fonts;


%% Extended Data Fig. 2: Anatomical location of sensory- and motor-related ROIs is consistent across fish
% % load the data
% load([pathname_whole_brain_imaging_integrators 'significant_ROIs.mat'],'signif_ROIs');
% [sz, rez] = read_nrrd_metadata([pathname_reference_brains 'PortuguesLab_wholebrain_ref_gamma_masked.nrrd']);
% all_fish = dm_dir([pathname_whole_brain_imaging_integrators 'behavior\*_f*_behavior.mat']);
% all_fish = strrep(all_fish,'_behavior.mat','');
% n_fish = length(all_fish);
% ROI_coord_all = [];
% sensmot_clust_all = [];
% time_constants_all = [];
% progressbar('Loading data for Extended Data Fig. 2...');
% for f = 1:n_fish
%     fish_id = all_fish{f};
%     load([pathname_whole_brain_imaging_integrators 'ROIs\' fish_id '_ROIs.mat'],'ROI_coord');
%     ROI_coord_all = [ROI_coord_all; ROI_coord];
%     load([pathname_whole_brain_imaging_integrators 'clustering\' fish_id '_clustering.mat'],'sensmot_clust');
%     sensmot_clust_all = [sensmot_clust_all; sensmot_clust];
%     load([pathname_whole_brain_imaging_integrators 'time_constants\' fish_id '_time_constants.mat'],'time_constants');
%     time_constants_all = [time_constants_all; time_constants];
%     progressbar(f/n_fish);
% end
% tf_clust = false(length(ROI_coord_all),1);
% tf_clust(sensmot_clust_all==1 & time_constants_all>1.5,1)=true;
% tf_clust(sensmot_clust_all==2,2)=true;
% tf_clust(sensmot_clust_all==1 & time_constants_all<=1.5,3)=true;
% clear time_constants_all sensmot_clust_all time_constants sensmot_clust ROI_coord fish_id n_fish all_fish
% 
% % Extended Data Fig. 2a: Anatomical reference
% % morphing anatomical regions from the Z-brain atlas
% all_regions=dm_dir([pathname_reference_brains 'morphed_regions_from_ZBrain_atlas\*.mat']);
% ref_brain_mask = nrrdread([pathname_reference_brains 'PortuguesLab_wholebrain_ref_mask.nrrd']);
% temp_tf=ref_brain_mask==0;
% for i=1:length(all_regions)
%     load([pathname_reference_brains 'morphed_regions_from_ZBrain_atlas\' all_regions{i}],'A');
%     A(temp_tf)=false;
%     [top_view,side_view]=three_orthogonal_views(A,rez(2));
%     top_view=255-uint8(top_view*255);
%     side_view=255-uint8(side_view*255);
%     panel_name=strrep(all_regions{i},'.mat','');
%     figure('name',['Extended Data Fig. 2a ' strrep(all_regions{i},'.mat','') ' top']);
%     imshow(top_view);
%     figure('name',['Extended Data Fig. 2a ' strrep(all_regions{i},'.mat','') ' side']);
%     imshow(side_view);
% end
% clear A temp_tf ref_brain_mask all_regions A coord top_view side_view; 
% 
% % Extended Data Fig. 2b: Consistent anatomical regions
% r=ones(sz);
% g=ones(sz);
% b=ones(sz);
% cols=[col_integrators;col_motor;col_sensors];
% for i=1:3
%     coord=ROI_coord_all(tf_clust(:,i) & signif_ROIs);
%     for ii=1:length(coord)
%         this_coord=coord{ii};
%         r(this_coord)=cols(i,1);
%         g(this_coord)=cols(i,2);
%         b(this_coord)=cols(i,3);
%     end
% end
% top_view = cat(3,sum(r,3),sum(g,3),sum(b,3));
% top_view = top_view/max(top_view(:));
% HSV = rgb2hsv(top_view);
% HSV(:,:,2) = min(HSV(:,:,2)*1.7,1);
% top_view = hsv2rgb(HSV);
% side_view = cat(3,squeeze(sum(r,2)),squeeze(sum(g,2)),squeeze(sum(b,2)));
% side_view=imresize3(side_view,[sz(1),round(sz(3)/rez(1)),3]);
% side_view=flip(side_view,2);
% side_view = side_view/max(side_view(:));
% HSV = rgb2hsv(side_view);
% HSV(:,:,2) = min(HSV(:,:,2)*1.5,1);
% side_view = hsv2rgb(HSV);
% figure('name','Extended Data Fig. 2b top');
% imshow(top_view);
% figure('name','Extended Data Fig. 2b side');
% imshow(side_view);
% clear HSV this_coord coord cols r g b side_view top_view tf_clust ROI_coord_all sz rez;


%% Extended Data Fig. 3: Treatment of Tg(PC:epNtr-tagRFP) larvae with metronidazole ablates the PCs
% % Extended Data Fig. 3a-c: confocal images
% % created in illustrator
% % Extended Data Fig. 3d: entropy
% A = readmatrix([pathname_MS_data 'PC_ablation_quantification\PC_ablation_quantification_entropy_x.csv'],'range','E2:F13');
% g = readmatrix([pathname_MS_data 'PC_ablation_quantification\PC_ablation_quantification_entropy_x.csv'],'range','B2:B13','output','string');
% A(8,:) = []; % discarded due to unstable acquisition that affects entropy computations
% g(8) = [];
% g = g=='-';
% figure('name','Extended Data Fig. 3d');
% axes('xlim',[0.6 2.4],'xtick',[1 2],'xticklabel',{'before treatment','after treatment'}); hold on;
% ylabel('entropy [bits]');
% ids = find(g)';
% scatter(ones(length(ids),1)-0.05,A(ids,1),10,'markerfacecolor',[0 0 0],'markeredgecolor','none');
% scatter(2*ones(length(ids),1)-0.05,A(ids,2),10,'markerfacecolor',[0 0 0],'markeredgecolor','none');
% for i=ids
%     line([1 2]-0.05,[A(i,1) A(i,2)],'color','k','linestyle',':');
% end
% line([-0.05 0.05]+1-0.05,[1 1]*median(A(ids,1)),'color','k');
% line([0 0]+1-0.05,[prctile(A(ids,1),25) prctile(A(ids,1),75)],'color','k');
% line([-0.05 0.05]+2-0.05,[1 1]*median(A(ids,2)),'color','k');
% line([0 0]+2-0.05,[prctile(A(ids,2),25) prctile(A(ids,2),75)],'color','k');
% 
% ids = find(~g)';
% scatter(ones(length(ids),1)+0.05,A(ids,1),10,'markerfacecolor',col_PC,'markeredgecolor','none');
% scatter(2*ones(length(ids),1)+0.05,A(ids,2),10,'markerfacecolor',col_PC,'markeredgecolor','none');
% for i=ids
%     line([1 2]+0.05,[A(i,1) A(i,2)],'color',col_PC,'linestyle',':');
% end
% line([-0.05 0.05]+1+0.05,[1 1]*median(A(ids,1)),'color',col_PC);
% line([0 0]+1+0.05,[prctile(A(ids,1),25) prctile(A(ids,1),75)],'color',col_PC);
% line([-0.05 0.05]+2+0.05,[1 1]*median(A(ids,2)),'color',col_PC);
% line([0 0]+2+0.05,[prctile(A(ids,2),25) prctile(A(ids,2),75)],'color',col_PC);
% dm_fix_fig_fonts;

%% Extended Data Fig. 4: Acute reaction is not impaired after PC ablation
acute_PC_neg=load([pathname_behavior 'acute_reaction_experiment/treatment_control_group/pooled_data.mat'],'behavior');
acute_PC_neg=acute_PC_neg.behavior;
acute_PC_pos=load([pathname_behavior 'acute_reaction_experiment/PC_ablated_group/pooled_data.mat'],'behavior');
acute_PC_pos=acute_PC_pos.behavior;
data={acute_PC_neg acute_PC_pos};
panel_names={'Extended Data Fig. 4i', 'Extended Data Fig. 4ii', 'Extended Data Fig. 4iii', 'Extended Data Fig. 4iv'; 'Extended Data Fig. 4v', 'Extended Data Fig. 4vi', 'Extended Data Fig. 4vii', 'Extended Data Fig. 4viii'};
cond_names={'gain', 'lag', 'shunted_lag', 'gain_drop'};
par_names={'bout_duration';'next_interbout_duration'};
x_labels={'gain','lag [ms]','shunted lag [ms]','gain profile'};
y_labels={'mean bout duration [s]'; 'mean interbout duration [s]'};
y_lims=[0.3 0.6; 0.8 2];
for i=1:2
    for c=1:4
        plot_acute_complete_PC(panel_names{i,c},data,cond_names{c},par_names{i},[0 0 0; col_PC],x_labels{c},y_labels{i},y_lims(i,:));
        dm_fix_fig_fonts;
    end
end
clear acute_PC_neg acute_PC_pos i c panel_names cond_names par_names x_labels y_labels;


%% Extended Data Fig. 5: Long-term adaptation effects are detectable in the light-sheet functional imaging experiment
% % extract the data
% [~, n_fish, n_control, n_good_lag, n_bad_lag,...
%     tf_control_fish, tf_good_lag_fish, tf_bad_lag_fish,...
%     ~, ~, ~, ~,...
%     Scores_be] = extract_PC_imaging_data(pathname_PC_imaging,pathname_reference_brains);
% Scores_be_blocks = nan(n_fish,12);
% c=0;
% for i=1:10:111
%     c=c+1;
%     Scores_be_blocks(:,c)=nanmean(Scores_be(:,i:i+9),2);
% end
% 
% % Extended Data Fig. 5a: experimental protocol
% % taken from Fig. 6b
% 
% % Extended Data Fig. 5b: First bout duration in all trials
% fig_names = {'Extended Data Fig. 5bi','Extended Data Fig. 5bii','Extended Data Fig. 5biii'};
% temp = {Scores_be(tf_control_fish,:),Scores_be(tf_bad_lag_fish,:),Scores_be(tf_good_lag_fish,:)};
% temp_cols = [col_norm_reaf;col_lag_trained_bad;col_lag_trained];
% for i=1:3
%     plot_long_term_adaptation_PC_im(fig_names{i},temp(i),temp_cols(i,:),{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0.2 0.9]);
%     set(gca,'ytick',0.2:0.1:0.9);
%     dm_fix_fig_fonts;
% end
% 
% % Extended Data Fig. 5c: First bout duration in all trials
% % Extended Data Fig. 5ci: Acute reaction
% figure('name','Extended Data Fig. 5ci');
% axes('ylim',[-0.5 0.5],'ytick',-0.5:0.25:0.5,'xlim',[0 4],'xtick',1:3,'xticklabel',{'control','lag-trained non adapting','lag-trained adapting'}); hold on;
% xtickangle(45);
% p_value_PC_LS_acute = plot_quantif({Scores_be_blocks(tf_control_fish,3)-Scores_be_blocks(tf_control_fish,2),Scores_be_blocks(tf_bad_lag_fish,3)-Scores_be_blocks(tf_bad_lag_fish,2),Scores_be_blocks(tf_good_lag_fish,3)-Scores_be_blocks(tf_good_lag_fish,2)},'k','normalized mean bout duration [s]',[],'right');
% dm_fix_fig_fonts;
% % Extended Data Fig. 5cii: Back-to-baseline affect
% figure('name','Extended Data Fig. 5cii');
% axes('ylim',[-0.5 0.5],'ytick',-0.5:0.25:0.5,'xlim',[0 4],'xtick',1:3,'xticklabel',{'control','lag-trained non adapting','lag-trained adapting'}); hold on;
% xtickangle(45);
% p_value_PC_LS_back_to_base = plot_quantif({Scores_be_blocks(tf_control_fish,7)-Scores_be_blocks(tf_control_fish,3),Scores_be_blocks(tf_bad_lag_fish,7)-Scores_be_blocks(tf_bad_lag_fish,3),Scores_be_blocks(tf_good_lag_fish,7)-Scores_be_blocks(tf_good_lag_fish,3)},'k','normalized mean bout duration [s]',[],'right');
% dm_fix_fig_fonts;
% % Extended Data Fig. 5ciii: After-effect
% figure('name','Extended Data Fig. 5ciii');
% axes('ylim',[-0.5 0.5],'ytick',-0.5:0.25:0.5,'xlim',[0 4],'xtick',1:3,'xticklabel',{'control','lag-trained non adapting','lag-trained adapting'}); hold on;
% xtickangle(45);
% p_value_PC_LS_aftereffect = plot_quantif({Scores_be_blocks(tf_control_fish,8)-Scores_be_blocks(tf_control_fish,2),Scores_be_blocks(tf_bad_lag_fish,8)-Scores_be_blocks(tf_bad_lag_fish,2),Scores_be_blocks(tf_good_lag_fish,8)-Scores_be_blocks(tf_good_lag_fish,2)},'k','normalized mean bout duration [s]',[],'right');
% dm_fix_fig_fonts;
% clear fig_names temp temp_cols n_fish n_control n_good_lag n_bad_lag tf_control_fish tf_good_lag_fish tf_bad_lag_fish Scores_be Scores_be_blocks c


%% Extended Data Fig. 6: Activity of 0-0+ ROIs cannot be explained by behavior
% % extract the data
% [all_fish, n_fish, ~, ~, ~,...
%     ~, tf_good_lag_fish, ~,...
%     ~, ~, ~, fish_id_rois,...
%     Scores_be, Scores_im, Scores_mr,...
%     ~, Crit_taus_signif,...
%     Traces_bout_trig, Traces_mr_bout_trig, time_trig] = extract_PC_imaging_data(pathname_PC_imaging,pathname_reference_brains);
% fish_id_ex = '200303_f0';
% fish_num_ex = find(strcmp(all_fish,fish_id_ex));
% load([pathname_PC_imaging 'behavior\' fish_id_ex '_behavior.mat'],'time_be','grspeed','tail','bouts','meta');
% tf_exp_rois = Crit_taus_signif(:,1)==0 & Crit_taus_signif(:,2)==-1 & Crit_taus_signif(:,3)==0 & Crit_taus_signif(:,4)==1;
% tf_mot_rois = Crit_taus_signif(:,1)==1;
% Data_temp_exp_rois = nan(n_fish,120);
% Data_temp_mot_rois = nan(n_fish,120);
% trig_data = compute_mean_trig_data(Traces_bout_trig);
% Trig_data_mot_rois = nan(n_fish,size(trig_data,2),size(trig_data,3));
% Trig_data_exp_rois = nan(n_fish,size(trig_data,2),size(trig_data,3));
% for i=1:n_fish
%     Data_temp_mot_rois(i,:)=nanmean(Scores_im(fish_id_rois==i & tf_mot_rois,:),1);
%     Data_temp_exp_rois(i,:)=nanmean(Scores_im(fish_id_rois==i & tf_exp_rois,:),1);
%     Trig_data_mot_rois(i,:,:) = nanmean(trig_data(fish_id_rois==i & tf_mot_rois,:,:),1);
%     Trig_data_exp_rois(i,:,:) = nanmean(trig_data(fish_id_rois==i & tf_exp_rois,:,:),1);
% end
% Trig_data_mr = compute_mean_trig_data(Traces_mr_bout_trig);
% clear trig_data n_fish all_fish fish_id_rois Scores_im Crit_im_signif Traces_bout_trig Traces_mr_bout_trig tf_exp_rois tf_mot_rois
% 
% % Extended Data Fig. 6a: Motor regressor in some example trials
% trials_to_show = 19:22;
% tf=time_be>(trials_to_show(1)-1)*30 & time_be<=trials_to_show(end)*30;
% dt=time_be(2)-time_be(1);
% this_tail=tail(tf);
% this_grspeed=grspeed(tf);
% this_grspeed(1:10)=0;
% this_grspeed(end-9:end)=0;
% this_grspeed(end/2-10:end/2+10)=0;
% this_time=time_be(tf);
% bout_starts = nan(length(trials_to_show),1);
% bout_ends=bout_starts;
% c=0;
% for i=trials_to_show
%     c=c+1;
%     id=find(bouts.trial==i,1);
%     bout_starts(c)=bouts.start(id);
%     bout_ends(c)=bouts.end(id);
% end
% load([pathname_PC_imaging 'traces\' fish_id_ex '_traces.mat'],'trace_motor_regr');
% time_im=1/meta.F:1/meta.F:length(trace_motor_regr)/meta.F;
% tf=time_im>(trials_to_show(1)-1)*30 & time_im<=trials_to_show(end)*30;
% trace_motor_regr=trace_motor_regr(tf);
% time_im=time_im(tf);
% figure('name','Extended Data Fig. 6a');
% set(gcf,'position',[109 20 1647 420]);
% axes('ycolor','none','xlim',[this_time(1)-dt this_time(end)]/60); hold on;
% xlabel('time in experiment [min]');
% plot(this_time/60,this_tail,'color',col_motor,'linewidth',1);
% plot(this_time/60,this_grspeed/7-5,'color',col_sensory,'linewidth',0.5);
% plot(time_im/60,trace_motor_regr*2-10,'color','k','linewidth',0.5);
% line([9.1 9.1],[-7 -5]*2,'color','k');
% drawnow; y_lim=get(gca,'ylim');
% show_lt_trials_lines(9,10,y_lim(2),y_lim(2),col_pre);
% show_lt_trials_lines(10,11,y_lim(2),y_lim(2),col_adapt_start);
% for i=1:length(trials_to_show)
%     fill_bout (gca,bout_starts(i)/60,bout_ends(i)/60,y_lim(1),y_lim(2),0.3);  
% end
% dm_fix_fig_fonts;
% clear y_lim tail id grspeed fish_id_ex tf meta time_im trace_motor_regr bouts bout_ends bout_starts time_be this_time this_gr this_tail trials_to_show tf dt
% 
% % Extended Data Fig. 6b: Behavior
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6 top',{Scores_be(fish_num_ex,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0.2 1.2]);
% set(gca,'ytick',0.2:0.2:1.2);
% dm_fix_fig_fonts;
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6 bottom',{Scores_be(tf_good_lag_fish,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0 1.2],true);
% set(gca,'ytick',0:0.2:1.2);
% dm_fix_fig_fonts;
% clear Scores_be;
% 
% % Extended Data Fig. 6c: Responses of motor regressors
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6c top left',{Scores_mr(fish_num_ex,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0.4 1.6]);
% set(gca,'ytick',0.4:0.2:1.6);
% ylabel('bout-triggered response [sd]');
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6c top right',time_trig,Trig_data_mr(fish_num_ex,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2])
% dm_fix_fig_fonts;
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6c bottom left',{Scores_mr(tf_good_lag_fish,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0 2],true);
% ylabel('bout-triggered response [sd]');
% set(gca,'ytick',0:0.5:2);
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6c bottom right',time_trig,Trig_data_mr(tf_good_lag_fish,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2])
% dm_fix_fig_fonts;
% 
% % Extended Data Fig. 6d: Responses of 0-0+ ROIs
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6d top left',{Data_temp_exp_rois(fish_num_ex,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0 2]);
% ylabel('bout-triggered response [sd]');
% set(gca,'ytick',0:0.4:2);
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6d top right',time_trig,Trig_data_exp_rois(fish_num_ex,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2])
% dm_fix_fig_fonts;
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6d bottom left',{Data_temp_exp_rois(tf_good_lag_fish,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[-0.5 2],true);
% ylabel('bout-triggered response [sd]');
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6d bottom right',time_trig,Trig_data_exp_rois(tf_good_lag_fish,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2])
% dm_fix_fig_fonts;
% 
% % Extended Data Fig. 6e: Responses of motor ROIs
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6e top left',{Data_temp_mot_rois(fish_num_ex,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[0.2 2.2]);
% ylabel('bout-triggered response [sd]');
% set(gca,'ytick',0.2:0.4:2.2);
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6e top right',time_trig,Trig_data_mot_rois(fish_num_ex,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2.5])
% dm_fix_fig_fonts;
% plot_long_term_adaptation_PC_im('Extended Data Fig. 6e bottom left',{Data_temp_mot_rois(tf_good_lag_fish,:)},'k',{'-'},col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,[-0.5 2.5],true);
% ylabel('bout-triggered response [sd]');
% dm_fix_fig_fonts;
% plot_trigaver_PC('Extended Data Fig. 6e bottom right',time_trig,Trig_data_mot_rois(tf_good_lag_fish,:,:),col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,[-0.5 2])
% dm_fix_fig_fonts;
% clear c Trig_data_mr Scores_mr Trig_data_mot_rois Data_temp_mot_rois ans time_trig tf_good_lag_fish Trig_data_exp_rois fish_num_ex Data_temp_exp_rois;


%% Extended Data Fig. 7: 0-0+ ROIs represent a spatially distributed subpopulation of PCs
% % extract the data
% [~, n_fish, n_control, n_good_lag, n_bad_lag,...
%     ~, ~, ~,...
%     tf_control_rois, tf_good_lag_rois, tf_bad_lag_rois, fish_id_rois,...
%     ~, ~, ~,...
%     ~, Crit_taus_signif,...
%     ~, ~, ~,...
%     ROI_coord,sz_PC_ref,rez_PC_ref] = extract_PC_imaging_data(pathname_PC_imaging,pathname_reference_brains);
% tf_exp_rois = Crit_taus_signif(:,1)==0 & Crit_taus_signif(:,2)==-1 & Crit_taus_signif(:,3)==0 & Crit_taus_signif(:,4)==1;
% tf_mot_rois = Crit_taus_signif(:,1)==1;
% 
% % shuffle ROI labels 100 times and build shuffled maps
% % n_boots=100;
% % tf_exp_rois_shuf = false(length(tf_exp_rois),n_boots);
% % tf_mot_rois_shuf = tf_exp_rois_shuf;
% % for f=1:n_fish
% %     ids_this_fish_rois = find(fish_id_rois==f);
% %     for i=1:n_boots
% %         tf_exp_rois_shuf(ids_this_fish_rois,i)=tf_exp_rois(ids_this_fish_rois(randperm(length(ids_this_fish_rois))));
% %         tf_mot_rois_shuf(ids_this_fish_rois,i)=tf_mot_rois(ids_this_fish_rois(randperm(length(ids_this_fish_rois))));
% %     end
% % end
% % A_exp_control = zeros(sz_PC_ref);
% % A_exp_bad_lag = zeros(sz_PC_ref);
% % A_exp_good_lag = zeros(sz_PC_ref);
% % A_mot_good_lag = zeros(sz_PC_ref);
% % progressbar('Building shuffled maps...');
% % for i=1:n_boots
% %     a=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_control,tf_control_rois,tf_exp_rois_shuf(:,i),fish_id_rois,ROI_coord,1,false);
% %     A_exp_control=A_exp_control+a;
% %     a=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_bad_lag,tf_bad_lag_rois,tf_exp_rois_shuf(:,i),fish_id_rois,ROI_coord,1,false);
% %     A_exp_bad_lag=A_exp_bad_lag+a;
% %     a=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_good_lag,tf_good_lag_rois,tf_exp_rois_shuf(:,i),fish_id_rois,ROI_coord,1,false);
% %     A_exp_good_lag=A_exp_good_lag+a;
% %     a=build_PC_map(sz_PC_ref,rez_PC_ref,n_fish,n_good_lag,tf_good_lag_rois,tf_mot_rois_shuf(:,i),fish_id_rois,ROI_coord,1,false);
% %     A_mot_good_lag=A_mot_good_lag+a;
% %     progressbar(i/n_boots);
% % end
% % A_exp_control = A_exp_control/n_boots;
% % A_exp_bad_lag = A_exp_bad_lag/n_boots;
% % A_exp_good_lag = A_exp_good_lag/n_boots;
% % A_mot_good_lag = A_mot_good_lag/n_boots;
% % save([path_to_data_PC_imaging 'shuffled_maps.mat'],'A_exp_control','A_exp_bad_lag','A_exp_good_lag','A_mot_good_lag');
% load([pathname_PC_imaging 'shuffled_maps.mat'],'A_exp_control','A_exp_bad_lag','A_exp_good_lag','A_mot_good_lag');
% 
% % plot the maps
% show_PC_map('Extended Data Fig. 7 control',A_exp_control,rez_PC_ref,50);
% dm_fix_fig_fonts;
% show_PC_map('Extended Data Fig. 7 bad lag',A_exp_bad_lag,rez_PC_ref,50);
% dm_fix_fig_fonts;
% show_PC_map('Extended Data Fig. 7 good lag',A_exp_good_lag,rez_PC_ref,50);
% dm_fix_fig_fonts;
% clear tf_mot_rois tf_exp_rois A_exp_control A_exp_bad_lag A_exp_good_lag A_mot_good_lag n_fish n_control n_good_lag n_bad_lag tf_control_rois tf_good_lag_rois tf_bad_lag_rois fish_id_rois Crit_im_signif ROI_coord sz_PC_ref rez_PC_ref
% 

%% Functions related to bout power
function [h, ballistic_end] = find_ballistic_end_acyte(data_struct,ids,time_power)
data_gain = extract_bout_power_acute(data_struct,'gain');
data_lag = extract_bout_power_acute(data_struct,'lag');
data_gain_drop = extract_bout_power_acute(data_struct,'gain_drop');
data = [data_gain;data_lag;data_gain_drop];
data=data(ids,:,:);
data = permute(data,[3,1,2]);
[~,N,T]=size(data);
p=nan(1,T);
if N>2
    for i=1:T
        p(i)=kruskalwallis(data(:,:,i),[],'off');
    end
else
    for i=1:T
        p(i)=signrank(data(:,1,i),data(:,2,i));
    end
end
h=p<0.05/T;
ballistic_end=time_power(find(h,1)-1);
end
function [h] = create_bout_power_axes(str,num_subplots)
figure('name',str);
for i=1:num_subplots
    h(i)=subplot(1,num_subplots,i,'xlim',[-0.1 1],'xtick',-0.1:0.1:1); hold on;
    xlabel('time relative to bout onset [s]');
    ylabel('mean bout power [au]')
end
end
function [data] = prepare_data_for_bout_power_quantif(cell_struct,trial_num,ids)
num_struct=length(cell_struct);
data=cell(num_struct,1);
for s=1:num_struct
    this_struct=cell_struct{s};
    data_temp=squeeze(this_struct(trial_num,:,:)-this_struct(2,:,:));
    data{s}=nanmean(data_temp(ids,:),1);
end
end
function [h] = plot_bout_power(h,time,data,ids,cols)
[m, p25, p75] = compute_median_and_quartilles (data);
num_plots=length(ids);
N=size(data,2);
x=nan(N,size(data,3),num_plots);
for i=1:num_plots
    fill(h,[time flip(time)],[p25(ids(i),:) flip(p75(ids(i),:))],cols(i,:),'edgecolor','none','facealpha',0.3);
    plot(h,time,m(ids(i),:),'color',cols(i,:));
    drawnow;
    line(h,[0 0],get(h,'ylim'),'color','k','linestyle',':');
    x(:,:,i)=squeeze(data(ids(i),:,:));
end
p=nan(num_plots-1,N);
for i=2:num_plots
    for j=1:N
        p(i-1,j)=signrank(x(j,:,1),x(j,:,i));
    end
end
h=p<0.05/N;%(N*(num_plots-1));
end
function [data]=extract_bout_power(data_struct)
N=length(data_struct);
[M,T]=size(data_struct(1).trials.power.first);
data=nan(M,T,N);
for i=1:N
    data(:,:,i)=data_struct(i).trials.power.first;
end
for i=1:M/10
    data(i,:,:)=nanmean(data((i-1)*10+1:i*10,:,:),1);
end
data=data(1:M/10,:,:);
end
function plot_bout_power_acute(panel_name,data_struct,par_name,cols,time_power,h,col_line,y_lim)
create_bout_power_axes(panel_name,1);
set(gca,'ylim',y_lim);
y=y_lim(2)-diff(y_lim)/10;
data=extract_bout_power_acute(data_struct,par_name);
data=nanmedian(data,3);
[M,T]=size(data);
for i=1:M
    plot(time_power,data(i,:),'color',cols(i,:),'linewidth',2);
end
dt=time_power(2)-time_power(1);
for i=1:T
    if h(i)
        line([time_power(i)-dt/2 time_power(i)+dt/2],[y y],'color',col_line,'linewidth',3);
    end
end
end
function [data] = extract_bout_power_acute(data_struct,par_name)
N=length(data_struct);
[M,T]=size(data_struct(1).bouts.condition_mean.(par_name).power.mean);
data=nan(M,T,N);
for i=1:N
    data(:,:,i)=data_struct(i).bouts.condition_mean.(par_name).power.mean;
end
end

%% other functions
function [anat_im] = small_anatomy_image(anatomy_stack,sz,ROI_coord_all,ex_fish_num,ROI_num,col)
temp = false(sz);
temp(ROI_coord_all{ex_fish_num}{ROI_num}) = true;
temp = max(temp,[],3);
temp = imfill(temp,'holes');
[x,y,z]=ind2sub(sz,ROI_coord_all{ex_fish_num}{ROI_num});
main_z = mode(z);
c = [round(mean(x(z==main_z))) round(mean(y(z==main_z)))];
anat_im = uint8(nanmean(anatomy_stack(:,:,unique(z)),3));
B = cell2mat(bwboundaries(temp));
anat_im = cat(3,anat_im,anat_im,anat_im);
for i=1:3
    for j=1:size(B,1)
        anat_im(B(j,1):B(j,1),B(j,2):B(j,2),i)=col(i)*255;
    end
end
anat_im = anat_im(c(1)-25:c(1)+25,c(2)-25:c(2)+25,:);
end
function [A] = build_PC_map(sz,rez,n_fish,n_group,tf_group_rois,tf_type_rois,fish_id_rois,ROI_coord,sigma,show_progressbar)
if nargin==9
    show_progressbar=true;
end
A=zeros(sz);
if show_progressbar
    progressbar('Building functional map...');
end
for f=1:n_fish
    ids=find(tf_group_rois & tf_type_rois & fish_id_rois==f)';
    if ~isempty(ids)
        A2=zeros(sz,'uint8');
        for i=ids
            A2(ROI_coord{i})=255;
        end
        A2=imgaussfilt3(A2,sigma./rez);
        A2=A2>0;
        A=A+A2;
    end
    if show_progressbar
        progressbar(f/n_fish);
    end
end
A=A/n_group*100;
end
function show_PC_map(str,A,rez_PC_ref,c_max)
un_vals=unique(A(:));
A=imresize3(A,[size(A,1),size(A,2),round(size(A,3)/rez_PC_ref(1))]);
if length(un_vals)<10
    A=roundtowardvec(A,un_vals);
end
figure('name',str,'position',[148 88 1462 866]);
subplot(2,2,1,'xcolor','none','ycolor','none','clim',[0 c_max],'ydir','reverse'); hold on;
im_top = max(A,[],3);
imagesc(im_top);
pbaspect([flip(size(im_top)) 1]);
line([10 10+100/rez_PC_ref(1)],[10 10],'color','k')
subplot(2,2,2,'xcolor','none','ycolor','none','clim',[0 c_max],'xdir','reverse','ydir','reverse'); hold on;
im_side = squeeze(max(A,[],2));
imagesc(im_side);
pbaspect([flip(size(im_side)) 1]);
subplot(2,2,3,'xcolor','none','ycolor','none','clim',[0 c_max]); hold on;
im_front = rot90(squeeze(max(A,[],1)));
imagesc(im_front);
pbaspect([flip(size(im_front)) 1]);
subplot(2,2,4,'xcolor','none','ycolor','none','clim',[0 c_max]);
clbr=colorbar;
ylabel(clbr,'percentage of fish with ROIs');
colormap(flip(hot));
end

function [all_fish, n_fish, n_control, n_good_lag, n_bad_lag,...
    tf_control_fish, tf_good_lag_fish, tf_bad_lag_fish,...
    tf_control_rois, tf_good_lag_rois, tf_bad_lag_rois, fish_id_rois,...
    Scores_be, Scores_im, Scores_mr,...
    Crit_im, Crit_im_signif,...
    Traces_bout_trig, Traces_mr_bout_trig, time_trig,...
    ROI_coord,sz_PC_ref,rez_PC_ref] = extract_PC_imaging_data(path_to_data_PC_imaging,path_to_reference_brain_stacks) %#ok<STOUT>
good_lag_fish_thresh = -0.04;
all_fish = dm_dir([path_to_data_PC_imaging 'behavior\' '*_f*_behavior.mat']);
all_fish=strrep(all_fish,'_behavior.mat','');
n_fish = length(all_fish);
tf_control_fish = false(n_fish,1);
tf_control_rois = false(0,1);
tf_good_lag_rois = false(0,1);
Scores_be=nan(n_fish,120);
Scores_im=[];
Scores_mr=nan(n_fish,120);
Crit_be=nan(n_fish,4);
Crit_im=[];
Crit_im_signif=[];
Traces_bout_trig=[];
Traces_mr_bout_trig=[];
fish_id_rois=[];
ROI_coord2={};
progressbar('Loading data for PC imaging experiment...');
for f=1:length(all_fish)
    fish_id = all_fish{f};
    load([path_to_data_PC_imaging 'behavior\' fish_id '_behavior.mat'],'meta','trials');
    tf_control_fish(f) = meta.lag_condition==0;
    load([path_to_data_PC_imaging 'criteria\' fish_id '_criteria.mat']); %#ok<LOAD>
    Crit_be(f,:)=crit_be;
    Crit_im=[Crit_im; crit_im]; %#ok<AGROW>
    Crit_im_signif = [Crit_im_signif; crit_im_signif]; %#ok<AGROW>
    n_rois_this_fish = size(crit_im_signif,1);
    tf_control_rois = [tf_control_rois; ones(n_rois_this_fish,1)*tf_control_fish(f)]; %#ok<AGROW>
    good_lag_fish=false;
    if ~tf_control_fish(f) && Crit_be(f,2)<good_lag_fish_thresh
        good_lag_fish=true;
    end
    tf_good_lag_rois = [tf_good_lag_rois; ones(n_rois_this_fish,1)*good_lag_fish]; %#ok<AGROW>
    load([path_to_data_PC_imaging 'scores\' fish_id '_scores.mat']); %#ok<LOAD>
    Scores_be(f,:)=scores_be;
    Scores_im = [Scores_im; scores_im]; %#ok<AGROW>
    Scores_mr(f,:) = scores_motor_regr;
    load([path_to_data_PC_imaging 'triggered_traces\' fish_id '_trig_traces.mat']); %#ok<LOAD>
    Traces_bout_trig = [Traces_bout_trig; traces_bout_trig]; %#ok<AGROW>
    Traces_mr_bout_trig = [Traces_mr_bout_trig; trace_motor_regr_bout_trig]; %#ok<AGROW>
    fish_id_rois=[fish_id_rois; f*ones(size(traces_bout_trig,1),1)]; %#ok<AGROW>
    load([path_to_data_PC_imaging 'ROIs\' fish_id '_ROIs.mat'],'ROI_coord');
    ROI_coord2=[ROI_coord2; ROI_coord]; %#ok<AGROW>
    progressbar(f/length(all_fish));
end
tf_good_lag_fish = ~tf_control_fish & Crit_be(:,2)<good_lag_fish_thresh;
tf_bad_lag_fish = ~tf_control_fish & ~tf_good_lag_fish;
tf_bad_lag_rois = ~tf_control_rois & ~tf_good_lag_rois;
n_control = sum(tf_control_fish);
n_good_lag = sum(tf_good_lag_fish);
n_bad_lag = sum(tf_bad_lag_fish);
ROI_coord=ROI_coord2;
[sz_PC_ref,rez_PC_ref]=read_nrrd_metadata([path_to_reference_brain_stacks 'PortuguesLab_PC_ref.nrrd']);
end
function [trig_data] = compute_mean_trig_data(traces_bout_trig)
trig_data_pre = nanmean(traces_bout_trig(:,:,11:20),3);
trig_data_adapt1 = nanmean(traces_bout_trig(:,:,21:30),3);
trig_data_adapt2 = nanmean(traces_bout_trig(:,:,61:70),3);
trig_data_post1 = nanmean(traces_bout_trig(:,:,71:80),3);
trig_data_post2 = nanmean(traces_bout_trig(:,:,111:120),3);
trig_data=cat(3,trig_data_pre,trig_data_adapt1,trig_data_adapt2,trig_data_post1,trig_data_post2);
end
function imagesc_scores(num_subplot,data,c_lim,col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end)
subplot(3,1,num_subplot,'xlim',[10 120]+0.5,'xtick',[11 20 30:10:120],'xticklabel',[],'ylim',[0 size(data,1)]+0.5,'ytick',[],'ydir','reverse'); hold on;
xtickangle(45);
ylabel('ROIs');
dm_imagesc(11:120,data(:,11:120));
colormap(gca,viridis);
colorbar;
set(gca,'clim',c_lim);
y_lim=get(gca,'ylim');
show_lt_trials_lines(10.5,20.5,y_lim(1),y_lim(1),col_pre,false);
show_lt_trials_lines(20.5,30.5,y_lim(1),y_lim(1),col_adapt_start,false);
show_lt_trials_lines(30.5,60.5,y_lim(1),y_lim(1),col_adapt,false);
show_lt_trials_lines(60.5,70.5,y_lim(1),y_lim(1),col_adapt_end,false);
show_lt_trials_lines(70.5,80.5,y_lim(1),y_lim(1),col_post_start,false);
show_lt_trials_lines(80.5,110.5,y_lim(1),y_lim(1),col_post,false);
show_lt_trials_lines(110.5,120.5,y_lim(1),y_lim(1),col_post_end,false);
end
function plot_acute_complete_PC(panel_name,data_struct,cond_name,par_name,cols,x_label,y_label,y_lim)
prepare_acute_axes(panel_name,data_struct{1},cond_name,x_label);
data_neg=prepare_data_acute(data_struct{1},cond_name,par_name);
data_pos=prepare_data_acute(data_struct{2},cond_name,par_name);
plot_acute({data_neg,data_pos},cols,y_label,y_lim);
end
function plot_acute(data,cols,y_label,y_lim)
n_cond=size(data{1},2);
n_datasets=length(data);
offsets=[-0.1 0.1];
x_ticks=get(gca,'xtick');
x_lim=get(gca,'xlim');
for i=1:n_datasets
    for c=1:n_cond
        x=x_ticks(c)+offsets(i);
        data_temp=data{i}(:,c);
        %         [m, p25, p75] = compute_median_and_quartilles (data_temp);
        %         h=plotSpread(data_temp,'spreadWidth',0.2,'xValues',x,'xMode','auto');
        %         h1=h{1};
        %         set(h1,'marker','o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','none','markersize',4);
        m=mean(data_temp,1);
        ste=nanstd(data_temp,[],1)/sqrt(sum(~isnan(data_temp),1));
        p25=m-ste; p75=m+ste;
        line([x-0.1 x+0.1],[m m],'color',cols(i,:),'linewidth',2);
        line([x x],[p25 p75],'color',cols(i,:),'linewidth',2);
    end
end
set(gca,'xlim',x_lim,'ylim',y_lim);
ylabel(y_label);
end
function [data] = prepare_data_acute_model(reaf_cond,best_par,dt,par_name)
N=size(best_par,2);
n_cond=size(reaf_cond,1);
data=nan(N,n_cond);
for i=1:N
    for c=1:n_cond
        switch par_name
            case 'bout_duration'
                data(i,c)=model_compute_parameters_v3 (model_short_trial_v3(best_par(:,i),dt,reaf_cond(c,:)),dt);
            case 'next_interbout_duration'
                [~, data(i,c)]=model_compute_parameters_v3 (model_short_trial_v3(best_par(:,i),dt,reaf_cond(c,:)),dt);
        end
    end
end
end
function [data] = prepare_data_acute(data_struct,cond_name,par_name)
N=length(data_struct);
n_cond=length(data_struct(1).bouts.condition_mean.(cond_name).conditions);
data=nan(N,n_cond);
for i=1:N
    data(i,:)=data_struct(i).bouts.condition_mean.(cond_name).(par_name).mean;
end

end
function [] = prepare_acute_axes(panel_name,data_struct,cond_name,x_label)
figure('name',panel_name);
x_ticks=data_struct(1).bouts.condition_mean.(cond_name).conditions;
n_cond=length(x_ticks);
temp=9/(n_cond+1);
axes('xtick',temp:temp:(9-temp),'xticklabel',x_ticks,'xlim',[0 9]); hold on;
xlabel(x_label);
end
function [] = plot_reaf_cond(reaf_cond, cartoon_time, cartoon_bout, cartoon_swim, cartoon_vigor, cartoon_pad_length, col_motor, col_sensory)
plot(cartoon_time,cartoon_bout*10,'color',col_motor);
for i=1:size(reaf_cond,2)
    cartoon_gr = compute_reaf (cartoon_swim, cartoon_vigor, reaf_cond(:,i))-10;
    plot(cartoon_time, cartoon_gr-i*5,'color',col_sensory);
end
fill_bout (gca,cartoon_pad_length+0.001,cartoon_pad_length+0.4,-20,0,0.3);
end
function [m, p25, p75] = compute_median_and_quartilles (data)
if length(size(data))==3
    m=nanmedian(data,3);
    p25=prctile(data,25,3);
    p75=prctile(data,75,3);
elseif length(size(data))==2
    m=nanmedian(data);
    p25=prctile(data,25);
    p75=prctile(data,75);
end
end

function [p] = make_lt_fbd_panel(str,y_label,cell_struct,str1,str2,trials_num,col,fish_num,trials_num0)
if nargin<=8
    trials_num0=11:20;
end
prepare_quantif_axes2(str);
data = prepare_data_for_quantif(cell_struct,str1,str2,trials_num,trials_num0);
if nargin>=8
    p = plot_quantif(data,col,y_label,fish_num);
else
    p = plot_quantif(data,col,y_label);
end
end

function [p] = make_lt_fbd_panel_PC(str,y_label,cell_struct,str1,str2,trials_num,col,tail,trials_num0)
if nargin <= 8
    trials_num0 = 11:20;
end
prepare_quantif_axes_PC(str);
data = prepare_data_for_quantif(cell_struct,str1,str2,trials_num,trials_num0);
p = plot_quantif_PC(data,col,y_label,tail);
end

function [] = prepare_quantif_axes2(str)
figure('name',str);
axes('xlim',[0 3],'xtick',1:2,'xticklabel',{'normal reafference control','lag-trained'}); hold on;
end

function [] = prepare_quantif_axes_PC(str)
figure('name',str);
axes('xlim',[0 6],'xtick',[1.5 4.5],'xticklabel',{'treatment control','PC-ablated'}); hold on;
end

function [data] = prepare_data_for_quantif(cell_struct,str1,str2,trials_num,trials_num0)
if nargin == 4
    trials_num0 = 11:20;
end
num_struct=length(cell_struct);
data=cell(num_struct,1);
for s=1:num_struct
    this_struct=cell_struct{s};
    N=length(this_struct);
    data_temp=nan(N,1);
    for i=1:N
        temp=this_struct(i).trials.(str1).(str2);
        data_temp(i)=nanmean(temp(trials_num)) - nanmean(temp(trials_num0));
    end
    data{s}=data_temp;
end
end

function [p, p0] = plot_quantif(data,col,y_label,fish_num,tail)
all_datasets=1:length(data);
h=plotSpread(data,'spreadWidth',0.7,'yLabel',y_label);
h1=h{1};
set(h1,'marker','o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','none','markersize',4);
if nargin==4
    if ~isempty(fish_num)
        a=get(h1(2));
        plot(a.XData(fish_num),a.YData(fish_num),'marker','o','markerfacecolor','none','markeredgecolor',[0.7 0 0],'markersize',6,'linewidth',2)
    end
end
p=nan(length(data));
p0=nan(length(data),1);
if nargin<5
    tail='both';
end
for s=all_datasets
    data_temp=data{s};
    [m, p25, p75] = compute_median_and_quartilles (data_temp);
    line([s-0.2 s+0.2],[m m],'color',col,'linewidth',2);
    line([s s],[p25 p75],'color',col,'linewidth',2);
    other_datasets=all_datasets;
    other_datasets(s)=[];
    for ss=other_datasets
        data_temp2=data{ss};
        p(s,ss)=ranksum(data_temp,data_temp2,'tail',tail);
    end
    p0(s)=signrank(data_temp);
end
end

function [p, p0] = plot_quantif_PC(data,col,y_label,tail)
% there must be 4 datasets: 
% 1. Treatment control: normal reafference control
% 2. Treatment control: lag-trained
% 3. PC-ablated: normal reafference control
% 4. PC-ablated control: lag-trained
all_datasets=1:length(data);
x_pos=[1 2 4 5];
cols_pos=[1 2 1 2];
h=plotSpread(data,'spreadWidth',0.7,'yLabel',y_label,'xValues',x_pos,'xMode','auto');
h1=h{1};
set(h1,'marker','o','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor','none','markersize',4);
p=nan(length(data));
p0=nan(length(data),1);
for s=all_datasets
    data_temp=data{s};
    [m, p25, p75] = compute_median_and_quartilles (data_temp);
    line([x_pos(s)-0.2 x_pos(s)+0.2],[m m],'color',col(cols_pos(s),:),'linewidth',2);
    line([x_pos(s) x_pos(s)],[p25 p75],'color',col(cols_pos(s),:),'linewidth',2);
    p0(s)=signrank(data_temp);
end
for i=1:4
    for j=1:4
        p(i,j) = ranksum(data{i},data{j},'tail',tail);
    end
end
end

function [cartoon_time, cartoon_bout, cartoon_vigor, cartoon_swim] = make_cartoon_bout(cartoon_pad_length)
cartoon_dt=0.001;
cartoon_bd=0.4;
F=20; % Hz
l=cartoon_bd/cartoon_dt;
cartoon_bout=sin(F*2*pi*cartoon_dt*(1:l));
bout_shape=zeros(1,l);
k=0.01;
t=0.98;
for i=2:l
    bout_shape(i)=t*k+t*bout_shape(i-1);
end
bout_shape(bout_shape>1)=1;
bout_shape(end/2+1:end)=flip(bout_shape(1:end/2));
bout_shape=bout_shape/max(bout_shape);
cartoon_bout=cartoon_bout.*bout_shape;
cartoon_bout=padarray(cartoon_bout,[0 cartoon_pad_length/cartoon_dt],0,'both');
cartoon_time=cartoon_dt:cartoon_dt:cartoon_bd+2*cartoon_pad_length;
cartoon_vigor=zeros(1,length(cartoon_time));
nn=0.04;
for i=round(nn/0.001):length(cartoon_time)
    cartoon_vigor(i)=std(cartoon_bout(i-(round(nn/0.001)-1):i));
end
[b,a] = butter(3,0.04,'low');
cartoon_vigor=filtfilt(b,a,cartoon_vigor);
cartoon_vigor(cartoon_vigor<0)=0;
cartoon_swim=false(1,length(cartoon_time));
cartoon_swim(cartoon_vigor>0.05)=true;
end

function [] = fill_bout (h,x1,x2,y1,y2,alpha)
fill(h,[x1 x2 x2 x1],[y1 y1 y2 y2],[0 0 0],'edgecolor','none','facealpha',alpha);
end

function [this_time, this_tail, this_gr, bout_starts, bout_ends] = get_example_trial_data(trial_num,time_be,tail,grspeed,this_struct)
tf=time_be>(trial_num-1)*30 & time_be<=trial_num*30;
dt=time_be(2)-time_be(1);
this_tail=tail(tf);
this_gr=grspeed(tf);
this_gr(1:10)=0;
this_gr(end-9:end)=0;
this_time=time_be(tf);
tf=this_struct.bouts.trial==trial_num;
bout_starts=this_struct.bouts.start(tf)-this_time(1)-7.5+dt;
bout_ends=this_struct.bouts.end(tf)-this_time(1)-7.5+dt;
this_time=this_time-this_time(1)-7.5+dt;
end

function [time_be, tail, grspeed] = get_example_fish(fish_id, path_to_data, path_str)
time_be=load([path_to_data 'time_array_behavior.mat'],'time_be');
time_be=time_be.time_be;
data=load([path_to_data path_str fish_id '_data.mat'],'tail','time_tail','swim','time_stim','grspeed');
tail=interp1(data.time_tail,data.tail,time_be);
swim=interp1(data.time_stim,double(data.swim),time_be)>0.5;
grspeed=interp1(data.time_stim,data.grspeed,time_be);
tail=tail-nanmedian(tail(~swim));
tail=tail/nanstd(tail(swim));
end

function [] = show_trial_start_and_end
line([0 0],ylim,'color','k','linestyle','--');
line([15 15],ylim,'color','k','linestyle','--');
end

function [] = create_single_trial_axes(str)
figure('name',str);
axes('xlim',[-7.5 22.5],'ycolor','none','xtick',[-7.5 0 5 10 15 22.5]); hold on;
xlabel('time relative to trial onset [s]');
end


function [] = plot_long_term_adaptation(str,cell_struct,str1,str2,col,str_line,col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post,y_lim)
figure('name',str);
axes('xlim',[10.5 240.5],'xtick',[11 20 30:10:220 230 240]); hold on;
xlabel('trial #');
ylabel('first bout duration [s]');
xtickangle(45);
num_struct=length(cell_struct);
M=length(cell_struct{1}(1).trials.(str1).(str2));
X=1:M;
N=zeros(1,num_struct);
for s=1:num_struct
    this_struct=cell_struct{s};
    N(s)=length(this_struct);
    data_temp=nan(N(s),M);
    for i=1:N(s)
        data_temp(i,:)=this_struct(i).trials.(str1).(str2);
    end
%         data_temp=data_temp-nanmean(data_temp(:,11:20),2);
    %     [m, p25, p75] = compute_median_and_quartilles (data_temp);
    m = nanmean(data_temp,1);
    er=nanstd(data_temp,[],1)./sqrt(sum(~isnan(data_temp),1));
    p25=m-er; p75=m+er;
    fill([X flip(X)],[p25 flip(p75)],col(s,:),'edgecolor','none','facealpha',0.2);
    plot(X,m,'color',col(s,:),'linewidth',1.5,'linestyle',str_line{s});
end
show_lt_trials_lines(10.5,20.5,y_lim(2),y_lim(2),col_pre);
show_lt_trials_lines(20.5,30.5,y_lim(2),y_lim(2),col_adapt_start);
show_lt_trials_lines(30.5,220.5,y_lim(2),y_lim(2),col_adapt);
show_lt_trials_lines(220.5,230.5,y_lim(2),y_lim(2),col_adapt_end);
show_lt_trials_lines(230.5,240.5,y_lim(2),y_lim(2),col_post);
set(gca,'ylim',y_lim);
end

function plot_trigaver_PC(str,time_trig,trig_data,col_pre, col_adapt_start, col_adapt_end, col_post_start, col_post_end,y_lim)
cols = [col_pre; col_adapt_start; col_adapt_end; col_post_start; col_post_end];
figure('name',str);
axes('xlim',[-0.8 1.2],'ylim',y_lim,'xtick',-0.8:0.4:1.2); hold on;
xlabel('time relative to first bout onset [s]');
ylabel({'mean z-scored fluorescence'; '(baseline subtracted) [sd]'});
for i=1:5
    m = nanmean(trig_data(:,:,i),1);
    if size(trig_data,1)>1
        er = nanstd(trig_data(:,:,i),[],1)./sqrt(sum(~isnan(trig_data(:,:,i)),1));
        p25=m-er; p75=m+er;
        fill([time_trig flip(time_trig)],[p25 flip(p75)],cols(i,:),'edgecolor','none','facealpha',0.2);
    end
    plot(time_trig,m,'color',cols(i,:),'linewidth',2);
end
line([0 0],ylim,'color','k','linestyle',':');
end

function [p] = plot_long_term_adaptation_PC_im(str,cell_data,col,str_line,col_pre,col_adapt_start,col_adapt,col_adapt_end,col_post_start,col_post,col_post_end,y_lim,plot_individual)
if nargin==12
    plot_individual=false;
end
figure('name',str);
axes('xlim',[10.5 120.5],'xtick',[11 20 30:10:120],'ylim',y_lim); hold on;
xtickangle(45);
xlabel('trial #');
ylabel('first bout duration [s]');
xtickangle(45);
num_datasets=length(cell_data);
M=size(cell_data{1},2);
X=1:M;
N=zeros(1,num_datasets);
for s=1:num_datasets
    data_temp=cell_data{s};
    N(s)=size(data_temp,1);
    m = nanmean(data_temp,1);
    if N(s)>1
        if plot_individual
            for i=1:N(s)
                plot(X,data_temp(i,:),'linewidth',1,'linestyle',str_line{s},'color',[0.75 0.75 0.75]);
            end
        end
        er=nanstd(data_temp,[],1)./sqrt(sum(~isnan(data_temp),1));
        p25=m-er; p75=m+er;
        fill([X flip(X)],[p25 flip(p75)],col(s,:),'edgecolor','none','facealpha',0.2);
    end
    p(s)=plot(X,m,'color',col(s,:),'linewidth',1.5,'linestyle',str_line{s}); %#ok<AGROW>
end
show_lt_trials_lines(10.5,20.5,y_lim(2),y_lim(2),col_pre);
show_lt_trials_lines(20.5,30.5,y_lim(2),y_lim(2),col_adapt_start);
show_lt_trials_lines(30.5,60.5,y_lim(2),y_lim(2),col_adapt);
show_lt_trials_lines(60.5,70.5,y_lim(2),y_lim(2),col_adapt_end);
show_lt_trials_lines(70.5,80.5,y_lim(2),y_lim(2),col_post_start);
show_lt_trials_lines(80.5,110.5,y_lim(2),y_lim(2),col_post);
show_lt_trials_lines(110.5,120.5,y_lim(2),y_lim(2),col_post_end);
end

function [] = show_lt_trials_lines(x1,x2,y1,y2,col,show_fill)
if nargin==5
    show_fill=true;
end
if x1==x2
    line([x1 x1],[y1 y2],'color',col,'linewidth',2);
    x_lim=get(gca,'xlim');
    if show_fill
        h=fill([x_lim flip(x_lim)],[y1 y1 y2 y2],col,'edgecolor','none','facealpha',0.2);
        uistack(h,'bottom');
    end
    line(x_lim,[y2 y2],'color','k','linestyle',':','linewidth',0.5);
elseif y1==y2
    line([x1 x2],[y1 y1],'color',col,'linewidth',2);
    y_lim=get(gca,'ylim');
    if show_fill
        h=fill([x1 x1 x2 x2],[y_lim flip(y_lim)],col,'edgecolor','none','facealpha',0.2);
        uistack(h,'bottom');
    end
    line([x2 x2],y_lim,'color','k','linestyle',':','linewidth',0.5);
end
end

function [top_view,side_view,front_view] = three_orthogonal_views(stack,xy_rez)
sz=size(stack);
top_view=double(max(stack,[],3));

side_view=squeeze(double(max(stack,[],2)));
side_view=imresize(side_view,[sz(1),round(sz(3)/xy_rez)]);
side_view=flip(side_view,2);

% front_view = squeeze(double(max(stack,[],1)));
%     
% front_view=imresize(front_view,[sz(2),round(sz(3)/xy_rez)]);
% front_view=flip(imrotate(front_view,-90),2);
end