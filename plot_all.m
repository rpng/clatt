close all
SIMorEXP = 1;
start = n_init_steps+1;
incr = 1;

% % trajectory
for irun = 1 %floor(rand*nRuns+.5)+1%1:nRuns
    
    figure, hold on
    
    if ~isempty(xL_true)
        plot(xL_true(1,:,irun), xL_true(2,:,irun), 'm+','Linewidth',2,'Markersize',10)
    end
    
    if nT>0, plot(squeeze(xT_true(1,1,:,irun)), squeeze(xT_true(2,1,:,irun)), 'g--o','Linewidth',1), end
    if nT>1, plot(squeeze(xT_true(1,2,:,irun)), squeeze(xT_true(2,2,:,irun)), 'm--s','Linewidth',1), end
    if nT>2, plot(squeeze(xT_true(1,3,:,irun)), squeeze(xT_true(2,3,:,irun)), 'c-->','Linewidth',1), end
    if nT>3, plot(squeeze(xT_true(1,4,:,irun)), squeeze(xT_true(2,4,:,irun)), 'y--d','Linewidth',1), end
        
    if nR>0,plot(squeeze(xR_true(1,1,:,irun)), squeeze(xR_true(2,1,:,irun)), 'r-','Linewidth',2),end
    if nR>1,plot(squeeze(xR_true(1,2,:,irun)), squeeze(xR_true(2,2,:,irun)), 'k--','Linewidth',2),end
    if nR>2,plot(squeeze(xR_true(1,3,:,irun)), squeeze(xR_true(2,3,:,irun)), 'b-.','Linewidth',2),end
    if nR>3,plot(squeeze(xR_true(1,4,:,irun)), squeeze(xR_true(2,4,:,irun)), 'c-+','Linewidth',1),end
    if nR>4,plot(squeeze(xR_true(1,5,:,irun)), squeeze(xR_true(2,5,:,irun)), 'm-*','Linewidth',1),end
    if nR>5,plot(squeeze(xR_true(1,6,:,irun)), squeeze(xR_true(2,6,:,irun)), 'y-x','Linewidth',1),end
    
    if ~isempty(xL_true)
        hh=legend('Landmarks', 'Target', 'Robot');
    elseif nT==1 && nR==2
        hh=legend('Target 1',   'Robot 1','Robot 2');
    elseif nT==2 && nR==4
        hh=legend('Target 1', 'Target 2', 'Robot 1','Robot 2','Robot 3','Robot 4');
    elseif nT==3 && nR==6
        hh=legend('Target 1', 'Target 2','Target 3', 'Robot 1','Robot 2','Robot 3','Robot 4','Robot 5','Robot 6');
    elseif nT==4 && nR==4
        hh=legend('Target 1', 'Target 2','Target 3', 'Target 4', 'Robot 1','Robot 2','Robot 3','Robot 4');
    else
        %
    end
    
    xlabel('x (m)','FontWeight','bold','fontsize',12), ylabel('y (m)','FontWeight','bold','fontsize',12)
    
    if nR>0,plot(squeeze(xR_true(1,1,1,irun)), squeeze(xR_true(2,1,1,irun)), 'ro','Markersize',10,'Linewidth',2), end
    if nR>1,plot(squeeze(xR_true(1,2,1,irun)), squeeze(xR_true(2,2,1,irun)), 'ko','Markersize',10,'Linewidth',2),end
    if nR>2,plot(squeeze(xR_true(1,3,1,irun)), squeeze(xR_true(2,3,1,irun)), 'bo','Markersize',10,'Linewidth',2),end
    if nR>3,plot(squeeze(xR_true(1,4,1,irun)), squeeze(xR_true(2,4,1,irun)), 'co','Markersize',10,'Linewidth',2), end
    if nR>4,plot(squeeze(xR_true(1,5,1,irun)), squeeze(xR_true(2,5,1,irun)), 'mo','Markersize',10,'Linewidth',2), end
    
    if nT>0, plot(squeeze(xT_true(1,1,1,irun)), squeeze(xT_true(2,1,1,irun)), 'g>','Markersize',10,'Linewidth',2), end
    if nT>1, plot(squeeze(xT_true(1,2,1,irun)), squeeze(xT_true(2,2,1,irun)), 'm>','Markersize',10,'Linewidth',2),end
    if nT>2, plot(squeeze(xT_true(1,3,1,irun)), squeeze(xT_true(2,3,1,irun)), 'c>','Markersize',10,'Linewidth',2), end
    if nT>3, plot(squeeze(xT_true(1,4,1,irun)), squeeze(xT_true(2,4,1,irun)), 'y>','Markersize',10,'Linewidth',2), end
        
    axis equal
    
    if SIMorEXP
        print('-depsc2', ['traj_',int2str(irun)] )
    else
        print('-depsc2', 'traj_exp' )
    end
    %     close
end

% % Robot RMS
for ell = 1:nR
    figure
    subplot(2,1,1), hold on
    plot([start:incr:nSteps],squeeze(rmsRp_avg_isam2(ell,start:incr:end)),'g-.','Linewidth',2);
    plot([start:incr:nSteps],squeeze(rmsRp_avg_isam3(ell,start:incr:end)),'b-o','Linewidth',1);
    plot([start:incr:nSteps],squeeze(rmsRp_avg_bmap(ell,start:incr:end)),'r-','Linewidth',2);
    
    if SIMorEXP
        ylabel('Robot position RMSE (m)','FontWeight','bold','fontsize',12)
    else
        ylabel('Robot position error (m)','FontWeight','bold','fontsize',12)
    end
    
    hh = legend('Naive-UIS','UIS','MAP');
    
    subplot(2,1,2), hold on

    plot([start:incr:nSteps],squeeze(rmsRth_avg_isam2(ell,start:incr:end)),'g-.','Linewidth',2);
    plot([start:incr:nSteps],squeeze(rmsRth_avg_isam3(ell,start:incr:end)),'b-o','Linewidth',1);
    plot([start:incr:nSteps],squeeze(rmsRth_avg_bmap(ell,start:incr:end)),'r-','Linewidth',2);
        
    
    if SIMorEXP
        xlabel('Time (sec)','FontWeight','bold','fontsize',12),
        ylabel('Robot heading RMSE (rad)','FontWeight','bold','fontsize',12)
        
        if dim_target==6
            print('-depsc2',['robot_rms_bo_',int2str(ell)])
        else
            print('-depsc2',['robot_rms_',int2str(ell)])
        end
    else
        xlabel('Time (sec)','FontWeight','bold','fontsize',12),
        ylabel('Robot heading error (rad)','FontWeight','bold','fontsize',12)
        print('-depsc2',['robot_rms_exp_',int2str(ell)])
    end

end

RMS_Robot_Position = [ mean(rmsRp_avg_isam2(:,start:end),2), mean(rmsRp_avg_isam3(:,start:end),2), mean(rmsRp_avg_bmap(:,start:end),2) ]
RMS_Robot_Heading = [ mean(rmsRth_avg_isam2(:,start:end),2), mean(rmsRth_avg_isam3(:,start:end),2), mean(rmsRth_avg_bmap(:,start:end),2) ]


% % Target RMS
for ell = 1:nT
    figure
    subplot(2,1,1), hold on
    
    plot([start:incr:nSteps],squeeze(rmsT_avg_isam2(ell,start:incr:end)),'g-.','Linewidth',2);
    plot([start:incr:nSteps],squeeze(rmsT_avg_isam3(ell,start:incr:end)),'b-o','Linewidth',1);
    plot([start:incr:nSteps],squeeze(rmsT_avg_bmap(ell,start:incr:end)),'r-','Linewidth',2);
    
    ylim([0 5])
    
    if SIMorEXP
        ylabel('Target position RMSE (m)','FontWeight','bold','fontsize',12)
    else
        ylabel('Target position error (m)','FontWeight','bold','fontsize',12)
    end
    
    hh = legend('Naive-UIS','UIS','MAP');
    
    subplot(2,1,2), hold on

    plot([start:incr:nSteps],squeeze(rmsTv_avg_isam2(ell,start:incr:end)),'g-.','Linewidth',2);
    plot([start:incr:nSteps],squeeze(rmsTv_avg_isam3(ell,start:incr:end)),'b-o','Linewidth',1);
    plot([start:incr:nSteps],squeeze(rmsTv_avg_bmap(ell,start:incr:end)),'r-','Linewidth',2);
    
    ylim([0 .6])
    
    if SIMorEXP
        xlabel('Time (sec)','FontWeight','bold','fontsize',12),
        ylabel('Target velocity RMSE (m/sec)','FontWeight','bold','fontsize',12)
        if dim_target==6
            print('-depsc2',['target_rms_bo_',int2str(ell)])
        else
            print('-depsc2',['target_rms_',int2str(ell)])
        end
    else
        
        xlabel('Time (sec)','FontWeight','bold','fontsize',12),
        ylabel('Target velocity error (m/sec)','FontWeight','bold','fontsize',12)
        print('-depsc2',['target_rms_exp_',int2str(ell)])
    end
end

RMS_Target_Position = [ mean(rmsT_avg_isam2(:,start:end),2),mean(rmsT_avg_isam3(:,start:end),2), mean(rmsT_avg_bmap(:,start:end),2) ]
RMS_Target_Vel = [ mean(rmsTv_avg_isam2(:,start:end),2),mean(rmsTv_avg_isam3(:,start:end),2), mean(rmsTv_avg_bmap(:,start:end),2) ]


