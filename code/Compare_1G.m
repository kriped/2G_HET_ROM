clear variables; close all;

load("../data/tempFile.mat");
load("../data/time_signal.mat")
diffusion_term_2G = 1./LAMBDA.*(1/K_VALUE(1)-1./K_VALUE');
feedback_term_2G = -1./LAMBDA*FB.*PHID_PHI_eq_mat_PHI./PHID_F_PHI;
xenon_absorbtion_term_2G_1 = -1./LAMBDA.*sigmaX.*PHID_PHILOWER_PHI./(PHID_F_PHI.^2).*PHID_PHI;
iodine_creation_2G = PHID_GAMMAI_PHI./PHID_PHI;
xenon_creation_2G = PHID_GAMMAX_PHI./PHID_PHI;
xenon_absorption_2G_2 = -sigmaX * PHID_X0_PHI./PHID_PHI;
xenon_absorption_2G_3 = -sigmaX * PHID_PHIUPPER_PHI./PHID_F_PHI;

% hold on
% yyaxis left
% [ d_2G, ix_2G ] = min( abs( time_2G-3 ) ); %Find closest value to three seconds to avoid prompt drop
% state_values_2G(ix_2G,(1-1)*3+1);
% plot(time_2G(ix_2G:end)/3600,state_values_2G(ix_2G:end,([4 1 9]-1)*3+1),"LineWidth",2)
% xlim([0 70])
% ylim([-4E7,2E7])
% ax1 = gca;
% ax1.YColor = 'b';
% ax1.FontSize = 15;
% ylabel("Amplitude 2G-HET [cm^{-2}s^{-1}]",'Fontsize', 18)
hold on
yyaxis left
plot(time_2G(100:end)/3600,state_values_2G(100:end,([4 1 9]-1)*3+1),"LineWidth",2)
xlim([0 70])
%ylim([-4E7,2E7])
ylim([-3e6,2e6])
ax.YColor = 'b';
ylabel("Amplitude 2-G HET [cm^{-2}s^{-1}]",'Fontsize', 22)

load("C:/Users/kriped/Chalmers/Christophe Demaziere - XEROM/Matlab code/1G_HOM_REAL_MAIN_SC/input/ROM_input.mat")
load("C:/Users/kriped/Chalmers/Christophe Demaziere - XEROM/Matlab code/1G_HOM_REAL_MAIN_SC/Results/time_dev.mat")

% Diffusion_term_1G = v*(-XS_HOMO.D*bsq+sigmaE);
Diffusion_term_1G = v.*XS_HOMO.NF*(1/kn(1)-1./kn); 
%Diffusion_term_test_1 =  v.*XS_HOMO.NF*(1/K_VALUE(1)-1./K_VALUE');
%Diffusion_term_test_2 = 1./LAMBDA.*(1/kn(1)-1./kn);
feedback_term_1G = -v*FB_1G.*phi0mat;
xenon_absorbtion_term_1G_1 = -v*sigmaX.*phi0mat;
iodine_creation_1G = gammaI*XS_HOMO.NF/(nu*keff);
xenon_creation_1G = gammaX*XS_HOMO.NF/(nu*keff);
xenon_absorption_1G_2 = -sigmaX.*X0mat;
xenon_absorption_1G_3 = -sigmaX*phi0mat;

C11_2G = diag(diffusion_term_2G+feedback_term_2G);
C13_2G = diag(xenon_absorbtion_term_2G_1);
C21_2G = diag(iodine_creation_2G);
C22_2G = -lambdaI;
C31_2G = diag(xenon_creation_2G+xenon_absorption_2G_2);
C32_2G = lambdaI;
C33_2G = diag(-lambdaX+xenon_absorption_2G_3);

C11_1G = Diffusion_term_1G+feedback_term_1G;
C13_1G = xenon_absorbtion_term_1G_1;
C21_1G = iodine_creation_1G;
C22_1G = -lambdaI;
C31_1G = xenon_creation_1G+xenon_absorption_1G_2;
C32_1G = lambdaI;
C33_1G = -lambdaX+xenon_absorption_1G_3;
C11_1G = diag(C11_1G(1:10,1:10));
C13_1G = diag(C13_1G(1:10,1:10));
C31_1G = diag(C31_1G(1:10,1:10));
C33_1G = diag(C33_1G(1:10,1:10));

% C11_1G = diag(C11_1G);
% C13_1G = diag(C13_1G);
% C31_1G = diag(C31_1G);
% C33_1G = diag(C33_1G);
Dif_k = 1./K_VALUE-1/K_VALUE(1);
Dif_k(1:4);

%%
FBC=feedback_term_1G(1,1)/feedback_term_2G(1,1);
XEDC_1 = xenon_absorbtion_term_1G_1(1,1)/xenon_absorbtion_term_2G_1(1,1);
XEDC_2 = xenon_absorption_1G_2(1,1)/xenon_absorption_2G_2(1,1);
XEDC_3 = xenon_absorption_1G_3(1,1)/xenon_absorption_2G_3(1,1);

%%
[ d_1G, ix_1G ] = min( abs( tot_time-3/3600 ) ); %Find closest value to three seconds to avoid prompt drop
tot_sol(ix_1G,4);
yyaxis right
plot(tot_time(ix_1G:end),tot_sol(ix_1G:end,([2 1 7]-1)*3+1),"LineWidth",2)
%ylim([-2e5,2e5])
ylim([-3e5,2e5])
ylabel("Amplitude 1G-HOM [cm^{-2}s^{-1}]","Fontsize", 22)
grid on
hold off
xlabel("Time [h]",'Fontsize', 22)
legend("Mode 2 (2G-HET)", "Mode 1 (2G-HET)", "Mode 3 (2G-HET)","Mode 2 (1G-HOM)","Mode 1 (1G-HOM)","Mode 3 (1G-HOM)","Location","southeast","FontSize",22)
ax2 = gca;
ax2.FontSize = 18;


normalised_feedbackterm_1G = feedback_term_1G./feedback_term_1G(1,1);
normalised_feedbackterm_2G = feedback_term_2G./feedback_term_2G(1,1);
normalised_xenon_1G_1 = xenon_absorbtion_term_1G_1./xenon_absorbtion_term_1G_1(1,1);
normalised_xenon_2G_1 = xenon_absorbtion_term_2G_1./xenon_absorbtion_term_2G_1(1,1);
normalised_xenon_1G_2 = xenon_absorption_1G_2./xenon_absorption_1G_2(1,1);
normalised_xenon_2G_2 = xenon_absorption_2G_2./xenon_absorption_2G_2(1,1);
normalised_xenon_1G_3 = xenon_absorption_1G_3./xenon_absorption_1G_3(1,1);
normalised_xenon_2G_3 = xenon_absorption_2G_3./xenon_absorption_2G_3(1,1);
normalised_iodine_2G = iodine_creation_2G./iodine_creation_2G(1,1);
normalised_xenon_2G = xenon_creation_2G./xenon_creation_2G(1,1);
%green_red_Heatmap(PHID_GAMMAI_PHI)
%green_red_Heatmap(PHID_GAMMAX_PHI)

%% Calculate axial offset in 2G-HET

AO = (sum(MOD1(:,:,18:34,1).*SIGF1(:,:,18:34)+MOD2(:,:,18:34,1).*SIGF2(:,:,18:34),'all')-sum(MOD1(:,:,1:17,1).*SIGF1(:,:,1:17)+MOD2(:,:,1:17,1).*SIGF2(:,:,1:17),'all'))/(sum(MOD1(:,:,:,1).*SIGF1(:,:,:)+MOD2(:,:,:,1).*SIGF2(:,:,:),'all'));


%Calculate radial offset in 2G-HET
RO = (sum(MOD1(:,17:32,:,1).*SIGF1(:,17:32,:)+MOD2(:,17:32,:,1).*SIGF2(:,17:32,:),'all')-sum(MOD1(:,1:16,:,1).*SIGF1(:,1:16,:)+MOD2(:,1:16,:,1).*SIGF2(:,1:16,:),'all'))/(sum(MOD1(:,:,:,1).*SIGF1(:,:,:)+MOD2(:,:,:,1).*SIGF2(:,:,:),'all'));




%% het radial
% fig = figure("Visible","off");
% hold off

% surf(MOD1(:,:,15,1))
% zlim([0,1])
% saveas(fig,"../shape_comparrison/mode_1_fast_rad.png")
% %set(gca,'visible','off')
% %% het axial adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1(16,16,:,1);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_1_fast_ax.png")
% %set(gca,'visible','off')
% %% het radial adjoint
% surf((MOD1_adj(:,:,15,1)))
% saveas(fig,"../shape_comparrison/mode_1_fast_rad_adj.png")
% zlim([0,1])
% grid on
% %set(gca,'visible','off')
% %shading interp
% %% het axial adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1_adj(16,16,:,1);
% plot(line,height,'LineWidth',2)
% saveas(fig,"../shape_comparrison/mode_1_fast_ax_adj.png")
% xlim([0,1])
% grid on
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf(MOD1(:,:,15,4))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_fast_rad.png")
% %set(gca,'visible','off')
% %% het axial 2 adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1(16,16,:,4);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_fast_ax.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf((MOD1_adj(:,:,10,4)))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_fast_rad_adj.png")
% %set(gca,'visible','off')
% %shading interp
% %% het axial 2 adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1_adj(16,16,:,4);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_fast_ax_adj.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf(MOD1(:,:,17,9))
% zlim([0,1])
% saveas(fig,"../shape_comparrison/mode_3_fast_rad.png")
% %set(gca,'visible','off')
% %% het axial 2 adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1(16,16,:,9);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_fast_ax.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf((MOD1_adj(:,:,17,9)))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_fast_rad_adj.png")
% %set(gca,'visible','off')
% %shading interp
% %% het axial 2 adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD1_adj(16,16,:,9);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_fast_ax_adj.png")
% %set(gca,'visible','off')
% 
% 
% %% het radial Thermal
% fig = figure("Visible","off");
% hold off
% grid on
% surf(MOD2(:,:,17,1))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_1_thermal_rad.png")
% %set(gca,'visible','off')
% %% het axial adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2(16,16,:,1);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_1_thermal_ax.png")
% %set(gca,'visible','off')
% %% het radial adjoint
% surf((MOD2_adj(:,:,17,1)))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_1_thermal_rad_adj.png")
% %set(gca,'visible','off')
% %shading interp
% %% het axial adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2_adj(16,16,:,1);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_1_thermal_ax_adj.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf(MOD2(:,:,17,4))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_thermal_rad.png")
% %set(gca,'visible','off')
% %% het axial 2 adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2(16,16,:,4);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_thermal_ax.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf((MOD2_adj(:,:,17,4)))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_thermal_rad_adj.png")
% %set(gca,'visible','off')
% %shading interp
% %% het axial 2 adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2_adj(16,16,:,4);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_2_thermal_ax_adj.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf(MOD2(:,:,17,9))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_thermal_rad.png")
% %set(gca,'visible','off')
% %% het axial 2 adjoint
% length_phi = size(MOD1);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2(16,16,:,9);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_thermal_ax.png")
% %set(gca,'visible','off')
% %% het radial 2 adjoint
% surf((MOD2_adj(:,:,17,9)))
% zlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_thermal_rad_adj.png")
% %set(gca,'visible','off')
% %shading interp
% %% het axial 2 adjoint
% length_phi = size(MOD1_adj);
% height = linspace(1,34,length_phi(3));
% line = zeros(1,length_phi(3));
% line(:) = MOD2_adj(16,16,:,9);
% plot(line,height,'LineWidth',2)
% xlim([0,1])
% grid on
% saveas(fig,"../shape_comparrison/mode_3_thermal_ax_adj.png")

%% Plot axial functions
hold off
close all
on = 1;
if on  
    [sizex_1G,sizey_1G,sizez_1G] = size(phi(:,:,:,1));
    h_1G = linspace(1,height,sizez_1G);
    X0_PHI1 = X0.*MOD2(:,:,:,1); % X0*Phi_thermal fundamental mode
    [~,max_idx] = max(X0_PHI1(:,:,:,1),[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(X0_PHI1(:,:,:,1)),max_idx);
    line_X0_PHI1(:) = X0_PHI1(x_max,y_max,:);
    
    X0_PHI2 = X0 .* MOD2(:,:,:,4);% X0*Phi_thermal fundamental mode
    [~,max_idx] = max(X0_PHI2(:,:,:,1),[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(X0_PHI2(:,:,:,1)),max_idx);
    line_X0_PHI2(:) = X0_PHI2(x_max,y_max,:);
    
    % Get axial profile of MODE 1 xenon 2G-HET
    PHID_X0_PHI1 = MOD1_adj(:,:,:,1).* X0_PHI1./(DV.*sum(G2_inner_product(MOD_adj(:,:,:,1),MOD(:,:,:,1),'vector','vector'),'all'));
    [~,max_idx] = max(PHID_X0_PHI1(:,:,:,1),[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(PHID_X0_PHI1(:,:,:,1)),max_idx);
    line_PHID_X0_PHI1(:) = PHID_X0_PHI1(x_max,y_max,:);
    
    % Get axial Profile of MODE 2 xenon 2G-HET
    PHID_X0_PHI2 = MOD1_adj(:,:,:,1).* X0_PHI2./(DV.*sum(G2_inner_product(MOD_adj(:,:,:,1),MOD(:,:,:,1),'vector','vector'),'all'));
    [~,max_idx] = max(PHID_X0_PHI2(:,:,:,1),[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(PHID_X0_PHI2(:,:,:,1)),max_idx);
    line_PHID_X0_PHI2(:) = PHID_X0_PHI2(x_max,y_max,:);
    
    %Get axial profile of MODE 3 xenon 2G-HET
    PHID_X0_PHI3 = MOD1(:,:,:,1).*X0.*MOD2(:,:,:,9)./(DV.*sum(G2_inner_product(MOD_adj(:,:,:,1),MOD(:,:,:,1),'vector','vector'),'all'));
    line_PHID_X0_PHI3 = get_axial(PHID_X0_PHI3);  
    
    %Get axial profile of phi_eq MODE 1 2G-HET
    PHID_eq_PHI1 = (MOD1_adj(:,:,:,1).*MOD_EQ_1_scaled.*MOD1(:,:,:,1)+MOD2_adj(:,:,:,1).*MOD_EQ_2_scaled.*MOD2(:,:,:,1))/PHID_F_PHI(1);
    [~,max_idx] = max(PHID_eq_PHI1(:,:,:,1),[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(PHID_eq_PHI1(:,:,:,1)),max_idx);
    line_PHID_eq_PHI1(:) = PHID_eq_PHI1(x_max,y_max,:);
    
    PHID_eq_PHI2 = (MOD1_adj(:,:,:,1).*MOD_EQ_1_scaled.*MOD1(:,:,:,4)+MOD2_adj(:,:,:,1).*MOD_EQ_2_scaled.*MOD2(:,:,:,4))./PHID_F_PHI(1);
    line_PHID_eq_PHI2 = get_axial(PHID_eq_PHI2);
    
    PHID_eq_PHI3 = (MOD1_adj(:,:,:,1).*MOD_EQ_1_scaled.*MOD1(:,:,:,9)+MOD2_adj(:,:,:,1).*MOD_EQ_2_scaled.*MOD2(:,:,:,9))/PHID_F_PHI(1);
    line_PHID_eq_PHI3 = get_axial(PHID_eq_PHI3);
    

    
    PHI_Xeq_PHI_1G_1(:,:,:) = phi(:,:,:,1).*Xeq.*phi(:,:,:,1);
    line_PHI_Xeq_PHI_1G_1 = get_axial(PHI_Xeq_PHI_1G_1);
    
    PHI_Xeq_PHI_1G_2(:,:,:) = phi(:,:,:,1).*Xeq.*phi(:,:,:,2);
    line_PHI_Xeq_PHI_1G_2 = get_axial(PHI_Xeq_PHI_1G_2);
    
    PHI_PHI_EQ_PHI_1G_1 = phi(:,:,:,1).*phieq.*phi(:,:,:,1);
    line_PHI_PHI_EQ_PHI_1G_1 = get_axial(PHI_PHI_EQ_PHI_1G_1);
    
    PHI_PHI_EQ_PHI_1G_2 = phi(:,:,:,1).*phieq.*phi(:,:,:,2);
    line_PHI_PHI_EQ_PHI_1G_2 = get_axial(PHI_PHI_EQ_PHI_1G_2);
    
    PHI_Xeq_PHI_1G_3(:,:,:) = phi(:,:,:,1).*Xeq.*phi(:,:,:,7);
    line_PHI_Xeq_PHI_1G_3 = get_axial(PHI_Xeq_PHI_1G_3);
    
    PHI_PHI_EQ_PHI_1G_3 = phi(:,:,:,1).*phieq.*phi(:,:,:,7);
    line_PHI_PHI_EQ_PHI_1G_3 = get_axial(PHI_PHI_EQ_PHI_1G_3);
    
    PHI_PHI_1G = phi(:,:,:,1).*phi(:,:,:,1);
    line_PHI_PHI_1G = get_axial(PHI_PHI_1G);
    %FIGURE 1
%     figure(1)
%     hold on
%     plot(line_Eq_PHI1,h)
%     plot(line_X0_PHI1,h)
%     grid on
%     ylabel height
%     legend("\Phi_{eq}(fast)\times\Phi(1,fast)","X0\times\Phi(1,thermal)")
%     hold off
%     %FIGURE 2
%     figure(2)
%     hold on
%     plot(line_PHID_eq_PHI1,h)
%     plot(line_PHID_X0_PHI1,h)
%     grid on
%     ylabel height
%     legend("\Phi^{T}(1)\times\Phi_{eq}(matrix)\times\Phi(1)","\Phi^{T}(1,fast)X0\times\Phi(1,thermal)")
%     hold off
%     %FIGURE 3
%     figure(3)
%     hold on
%     plot(line_PHID_eq_PHI2,h)
%     plot(line_PHID_X0_PHI2,h)
%     grid on
%     ylabel height
%     legend("\Phi^{T}(1)\times\Phi_{eq}(matrix)\times\Phi(2)","\Phi^{T}(1,fast)\times X0\times\Phi(2,thermal)")
    %FIGURE 4
    DZ_1G = height/sizez_1G;
    INT_XO = DZ*sum(line_PHID_X0_PHI1);
    INT_PHI_EQ = DZ*sum(line_PHID_eq_PHI1); 
    INT_Xeq_1G = DZ_1G*sum(line_PHI_Xeq_PHI_1G_1);
    INT_PHI_eq_1G = DZ_1G*sum(line_PHI_PHI_EQ_PHI_1G_1);
    INT_PHI_PHI1G = DZ_1G*sum(line_PHI_PHI_1G);
    figure(1)
    hold on
    %plot(h,line_PHID_eq_PHI2/INT_PHI_EQ)
    %plot(h,line_PHID_X0_PHI2/INT_XO)
    grid on
    
   % area graph compararison of MODE 1
    area(h,line_PHID_eq_PHI1/INT_PHI_EQ,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h,line_PHID_X0_PHI1/INT_XO,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    area(h_1G,line_PHI_PHI_EQ_PHI_1G_1/INT_PHI_eq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h_1G,line_PHI_Xeq_PHI_1G_1/INT_Xeq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    legend("2G-HET \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_2G(1))),"2G-HET X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_2G_2(1))),"1G-HOM \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_1G(1))),"1G-HOM X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_2G_2(1))),'Fontsize', 22)
    xlabel('Height (cm)','Fontsize', 18)
    ax2 = gca;
    ax2.FontSize = 18;
    hold off
   
    figure(2)
    hold on 
    grid on
    % area graph comparisson of MODE 2 
    area(h,line_PHID_eq_PHI2/INT_PHI_EQ,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h,line_PHID_X0_PHI2/INT_XO,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    area(h_1G,-line_PHI_PHI_EQ_PHI_1G_2/INT_PHI_eq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h_1G,-line_PHI_Xeq_PHI_1G_2/INT_Xeq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    legend("2G-HET \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_2G(4))),"2G-HET X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_2G_2(4))),"1G-HOM \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_1G(2))),"1G-HOM X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_2G_2(2))),'Fontsize',22)
    xlabel('Height (cm)','Fontsize', 22)
    ax2 = gca;
    ax2.FontSize = 18;
    hold off
    
    figure(3)
    hold on 
    grid on
    % area graph comparisson of MODE 3
    area(h,line_PHID_eq_PHI3/INT_PHI_EQ,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h,line_PHID_X0_PHI3/INT_XO,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    area(h_1G,-line_PHI_PHI_EQ_PHI_1G_3/INT_PHI_eq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','-','LineWidth',2)
    area(h_1G,-line_PHI_Xeq_PHI_1G_3/INT_Xeq_1G,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','--','LineWidth',2)
    legend("2G-HET \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_2G(9))),"2G-HET X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_2G_2(9))),"1G-HOM \Phi_{eq} integral = "+sprintf("%0.2f",abs(normalised_feedbackterm_1G(7))),"1G-HOM X_{eq} integral = "+sprintf("%0.2f",abs(normalised_xenon_1G_2(7))),'Fontsize', 22,'Location','north')
    ax2 = gca;
    ax2.FontSize = 18;
    xlabel('Height (cm)','Fontsize', 22)
    %display the difference between the top and bottom areas for both graphs respectively
    
    
    % compare the mode 1 feedback terms between the 1G-HOM and 2G-HET
    feedback_integrant_1G = phi(:,:,:,1).^2.*phieq;
    feedback_integrand_fast_2G = MOD1_adj(:,:,:,1).*MOD_EQ_1_scaled.*MOD1(:,:,:,1);
    feedback_integrand_thermal_2G = MOD2_adj(:,:,:,1).*MOD_EQ_2_scaled.*MOD2(:,:,:,1);
    feedback_integrand_2G = MOD1_adj(:,:,:,1).*MOD_EQ_1_scaled.*MOD1(:,:,:,1)+MOD2_adj(:,:,:,1).*MOD_EQ_2_scaled.*MOD2(:,:,:,1);
    line_feedback_integrant_1G = get_axial(feedback_integrant_1G);
    line_feedback_integrand_fast_2G = get_axial(feedback_integrand_fast_2G);
    line_feedback_integrand_thermal_2G = get_axial(feedback_integrand_thermal_2G);
    line_feedback_integrant_2G = get_axial(feedback_integrand_2G);
    figure(4)
    hold on
    grid on
    feedback_norm_constant_1G = v/(nodevol*sum(phi(:,:,:,1).^2,'all'));
    feedback_norm_const_2G = 1/(DV*sum(MOD1_adj(:,:,:,1).*1/v1.*MOD1(:,:,:,1)+MOD2_adj(:,:,:,1).*1/v2.*MOD2(:,:,:,1),'all'));
    area(h,line_feedback_integrant_2G/feedback_norm_const_2G,'FaceAlpha',0.3,'FaceColor',[0 0.4470 0.7410],'LineStyle','-')
    area(h,line_feedback_integrand_fast_2G/feedback_norm_const_2G,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','--')
    area(h,line_feedback_integrand_thermal_2G/feedback_norm_const_2G,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'LineStyle','-.')
    area(h_1G,line_feedback_integrant_1G/feedback_norm_constant_1G,'FaceAlpha',0.3,'FaceColor',[0.8500 0.3250 0.0980],'LineStyle','-')
    legend("2G-HET full","2G-HET fast","2G-HET thermal","1G-HOM")
    xlabel('Height (cm)','Fontsize', 12)
end

%% compare coefficients

Het_ratio_1G = feedback_term_1G(1,1)./xenon_absorption_1G_2(1,1)
Het_ratio_2G = feedback_term_2G(1,1)./xenon_absorption_2G_2(1,1)

Het_ratio_1G./Het_ratio_2G

%% create homogenous 2-G model to compare with
FLX_HOMO1 = sum(FLX1,'all'); %Homogenised Homogenise Fast flux
FLX_HOMO2 = sum(FLX2,'all'); %Homogenised Thermal flux

HOMO_SA1 = sum(ABS1.*FLX1,'all')./FLX_HOMO1; %Homogenised Absorbtion XS fast 
HOMO_SA2 = sum(ABS2.*FLX2,'all')./FLX_HOMO2; %Homogenised Absorbtion XS thermal

HOMO_NF1 = sum(NUFIS1.*FLX1,'all')./FLX_HOMO1; % Homogenised Nu*Fission XS fast
HOMO_NF2 = sum(NUFIS2.*FLX2,'all')./FLX_HOMO2; % Homogenised Nu*Fission XS thermal

HOMO_KF1 = sum(XS.KN.*NUFIS1.*FLX1,'all')./FLX_HOMO1; % Homogenised Kappa * Fission XS fast
HOMO_KF2 = sum(XS.KN.*NUFIS2.*FLX2,'all')./FLX_HOMO2; % Homogenised Kappa * Fission XS thermal

HOMO_D1 = sum(D1.*FLX1,'all')./FLX_HOMO1; 
HOMO_D2 = sum(D1.*FLX2,'all')./FLX_HOMO2;

HOMO_SR = sum(XS.SR1.*FLX1,'all')./FLX_HOMO1;

v1 = 2E9; % fast neutron velocity in cm/s 
v2 = 2.2E5; % thermal neutron velocity in cm/s
FB = 1.3e-18;
a = DZ*34/2;
Bg = pi/(2*a);
%% Test feedback term coefficient in a homoegnised 2G version
h = 1:DZ:2*a;
fphieq = PS*cos(Bg*h).^3;
fnorm = cos(Bg*h).^2;
fxeq = cos(Bg*h).^2*PS*(gammaI+gammaX).*(HOMO_NF1+HOMO_NF2*HOMO_SR/(HOMO_D2.*Bg^2+HOMO_SA2)).*cos(Bg*h)/(lambdaX+PS*sigmaX*HOMO_SR/(HOMO_D2*Bg^2+HOMO_SA2).*cos(Bg*h));

expr1 = FB*(1+(HOMO_NF2*HOMO_SR^2)/(keff*(Bg^2*HOMO_D2+HOMO_SA2)^3)*DZ*sum(fphieq,'all')); % numerator
expr2 = (1/v1+(HOMO_NF2*HOMO_SR)/(keff*(Bg^2*HOMO_D2+HOMO_SA2)^2*v2))*DZ*sum(fnorm,'all'); % denominator

coef_2G_feedback = expr1/expr2; % coefficient for the 2G version

FLX_CON = FLX1+FLX2; %Condensed Flux
FLX_HOMO = sum(FLX_CON,'all'); %Homogenised and condensed Flux

CON_v_1G = (1/v1*FLX1+1/v2*FLX2)./FLX_CON; % condense and homogenise the neutron speeds
FLX_CON(isnan(FLX_CON))=0;
CON_v_1G(isnan(CON_v_1G))=0;
HOMO_v_1G = sum(CON_v_1G.*FLX_CON,'all')./FLX_HOMO;

v_1G = 1/HOMO_v_1G;
coef_1G_feedback = FB*v_1G; % coefficient for the 1G version

coef_1G_feedback/coef_2G_feedback; % comparison of the two coefficients

%% Comparison for the xenon_eq term in a 2G homogenised test case


expr1 = sigmaX*HOMO_SR/(HOMO_D2*Bg^2+HOMO_SA2)*DZ*sum(fxeq,'all'); % feedback numerator
expr2 = (1+ (HOMO_SR*HOMO_NF2)/(keff*(HOMO_D2*Bg^2+HOMO_SA2)^2))*DZ*sum(fnorm,'all');% denominator

coef_2G_Xeq= expr1/expr2;

coef_1G_Xeq = sigmaX;

coef_1G_Xeq/coef_2G_Xeq;


%% Compare coefficient from feedback and xenon between 1G and 2G

coef_ratio_1G = coef_1G_feedback/coef_1G_Xeq;

coef_ratio_2G = coef_2G_feedback/coef_2G_Xeq;

coef_ratio_1G/coef_ratio_2G
