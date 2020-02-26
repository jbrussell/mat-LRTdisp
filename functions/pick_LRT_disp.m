function [ tr ] = pick_LRT_disp( absR_Tv, mat, per_trace, phv_trace, phv_trace_std )
% Plot LRT panel and pick points on dispersion curve
%
% J. Russell
% github.com/jbrussell

tr = [];
ipk = zeros(10,1);
all_ks = [];
while 1
    figure(97); clf;
    set(gcf,'Position',[151           1        1050         704]);
    imagesc(mat.per_vec, mat.phv_vec, absR_Tv); hold on;
    errorbar(per_trace,phv_trace,phv_trace_std,'or','MarkerFaceColor',[1 1 1],'linewidth',1.5,'markersize',7);
    caxis([0 1]);
    xlim([min(mat.per_vec) max(mat.per_vec)]);
    ylim([mat.v_min mat.v_max]);
    ylabel('Velocity (km/s)'); xlabel('Period (s)');
    set(gca,'YDir','normal','FontSize',15,'linewidth',1.5,'TickDir','out');
    colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
    
    for ii = 1:length(tr)
        plot(tr(ii).per,tr(ii).phv,'bx');
        str = {};
        for jj = 1:length(tr(ii).per)
            str{jj} = num2str(tr(ii).mode);
        end
        text(tr(ii).per,tr(ii).phv*1.05,str,'color','w','fontsize',15,'HorizontalAlignment', 'center');
    end
    
    % ENTER USER COMMANDS
    disp('Place curser over data point and select with branch number [0-9]:');
    disp('"s" for box selection of several points');
    disp('"r" to deselect point');
    disp('"q" to quit');
    [x y bot] = ginput(1);
    
    if bot == 'q'
        disp('Quitting');
        break;
    end
    
    if bot == 's'
        while 1
            disp('Click lower left corner of box:');
            [xll,yll]=ginput(1);
            disp('Click upper right corner of box:');
            [xur,yur]=ginput(1);
            rectangle('Position',[xll yll xur-xll yur-yll],'EdgeColor','r','linewidth',1.5);
            disp('Continue? y/n')
            [~,~,bot]=ginput(1);
            if bot == 'y'
                break
            elseif bot == 'n'
                continue
            else
                disp('wrong option, reslect box');
                continue
            end
        end
        % Define mode branch
        disp('Input mode branch number 0-9 [0=Fundamental]')
        [~,~,bot]=ginput(1);
        mode = bot-48;
        disp(mode);
        
        % Index good points
        I_good = find(per_trace>xll & per_trace<xur &...
                      phv_trace>yll & phv_trace<yur) ;
        
        for ii = 1:length(I_good)
            k = I_good(ii);
            all_ks = [all_ks; k];
            ipk(mode+1) = ipk(mode+1)+1;
            tr(mode+1).mode = mode;
            tr(mode+1).per(ipk(mode+1)) = per_trace(k);
            tr(mode+1).phv(ipk(mode+1)) = phv_trace(k);
            tr(mode+1).phv_std(ipk(mode+1)) = phv_trace_std(k);   
        end
    end
    
    if bot == 'r'
        % Find closest point
        [k,dist] = dsearchn([per_trace', phv_trace'],[x, y]);
        if ismember(k,all_ks)
            disp('Removing point...');
        else
            continue
        end
        for itr = 1:length(tr)
            [~,I_remove] = find(tr(itr).per==per_trace(k) & tr(itr).phv==phv_trace(k));
            if isempty(I_remove)
                continue
            end
            tr(itr).per(I_remove) = [];
            tr(itr).phv(I_remove) = [];
            tr(itr).phv_std(I_remove) = [];
            if isempty(tr(itr).per)
                tr(itr).mode = [];
            end
        end
        continue
    end
    
    if bot>=48 && bot<=57 % bot is a number from 0-9 
        mode = bot-48;
        disp(mode);
        % Find closest point
        [k,dist] = dsearchn([per_trace', phv_trace'],[x, y]);
%         if ismember(k,all_ks)
%             disp('Point already selected. Try again...');
%             continue
%         end
        all_ks = [all_ks; k];
        ipk(mode+1) = ipk(mode+1)+1;
        tr(mode+1).mode = mode;
        tr(mode+1).per(ipk(mode+1)) = per_trace(k);
        tr(mode+1).phv(ipk(mode+1)) = phv_trace(k);
        tr(mode+1).phv_std(ipk(mode+1)) = phv_trace_std(k);
    end
end

% Sort final picks and remove duplicates
for itr = 1:length(tr)
    [~, isrt] = sort(tr(itr).per,'ascend');
    tr(itr).per = tr(itr).per(isrt);
    tr(itr).phv = tr(itr).phv(isrt);
    tr(itr).phv_std = tr(itr).phv_std(isrt);
    
    [~,iunique] = unique(tr(itr).per);
    tr(itr).per = tr(itr).per(iunique);
    tr(itr).phv = tr(itr).phv(iunique);
    tr(itr).phv_std = tr(itr).phv_std(iunique);
end


end

