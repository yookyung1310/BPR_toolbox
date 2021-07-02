function [fh, axg] = bprcorrSummary(edBC, varargin)
% It shows BPR correction result edBC.summary
% created at 2019.08.27 by Kyung Yoo
% last update at 2019.08.28 by Kyung Yoo
% It uses edBC whole data as input now. It will be better to get edBC.edBC.summary only as an input. Update it later.
nS = numel(edBC.summary);

if nargin > 1
    str_sub = varargin{1};
    isStrSub = true;
else
    isStrSub = false;
end
if nargin > 2
    figpath = varargin{2};
    isFigSave = true;
else
    isFigSave = false;
end
if nargin > 3
    nRow = varargin{3}(1); nCol = varargin{3}(1);
else
    [nRow, nCol] = minsquare(nS+1);
end

g_inv = @(x,lambda) exp(log(lambda*x+1)/lambda);

% check h and SF estimation
fh(1) = figure(20000); clf; 
for s = 1:nS
    axg{1}(s) = subplot(nRow,nCol,s); hold on;

    x = edBC.summary{s}.tVec_selblr; 
    y = edBC.summary{s}.selblr;
    obj1 = plot(x, y, 'k');    
       
    lambda = edBC.summary{s}.lambda;
    x2 = edBC.summary{s}.tVec; 
    y2 = g_inv(edBC.summary{s}.mu2_hat_org', lambda); 
    obj2 = plot(x2, y2);    
    
    offset1 = y2(1); % offset just to overlap        
    obj3 = plot(x2, edBC.summary{s}.h_org+offset1, 'r');
    
    obj4 = plot(x2, y2+edBC.summary{s}.h_org, 'c');
    
    axis tight
    xlim([-2 5]);    
    vline(0);    

    if isStrSub, title(str_sub{s}); end
end
axes1 = subplot(nRow,nCol,nRow*nCol);
legend(axes1, [obj1 obj2 obj3 obj4], 'data','mu_n_o_n_B_P_R hat','h hat','mu_n_o_n_B_P_R+h hat');
titleStr = sprintf('SF and h estimation result of BPR correct');  mh = mtit(titleStr);
mh.th.Position(2) = mh.th.Position(2)+0.03;
xh = xlabel(mh.ah, 'time from blink offset (sec)'); set(xh, 'Visible', 'On')
yh = ylabel(mh.ah, 'PD (mm)'); set(yh, 'Visible', 'On')
if isFigSave, savefigpng(figpath, 'SF and h estimation in BPR', titleStr); end

% the same but show 0~3 sec only
fh(2) = figure(20001); clf; 
for s = 1:nS
    axg{2}(s) = subplot(nRow,nCol,s); hold on;

    x = edBC.summary{s}.tVec_selblr; 
    y = edBC.summary{s}.selblr;
    obj1 = plot(x, y, 'k');    
       
    lambda = edBC.summary{s}.lambda;
    x2 = edBC.summary{s}.tVec; 
    y2 = g_inv(edBC.summary{s}.mu2_hat_org', lambda); 
    obj2 = plot(x2, y2);    
    
    offset1 = y2(1); % offset just to overlap        
    obj3 = plot(x2, edBC.summary{s}.h_org+offset1, 'r');
    
    obj4 = plot(x2, y2+edBC.summary{s}.h_org, 'c');
    
    axis tight
    xlim([0 3]); 
    vline(0);    
    
    if isStrSub, title(str_sub{s}); end   
end
axes1 = subplot(nRow,nCol,nRow*nCol);
legend(axes1, [obj1 obj2 obj3 obj4], 'data','mu_n_o_n_B_P_R hat','h hat','mu_n_o_n_B_P_R+h hat');
titleStr = sprintf('SF and h estimation result of BPR correct [0-3 sec]');  mh = mtit(titleStr);
mh.th.Position(2) = mh.th.Position(2)+0.03;
xh = xlabel(mh.ah, 'time from blink offset (sec)'); set(xh, 'Visible', 'On')
yh = ylabel(mh.ah, 'PD (mm)'); set(yh, 'Visible', 'On')
if isFigSave, savefigpng(figpath, 'SF and h estimation in BPR', titleStr); end

% all subjects in a figure (unified scale)
color_post = [0 158 115]/255; % bluish green
fh(3) = figure(20100); clf;
for s = 1:nS
    x = edBC.summary{s}.tVec;    
    
    axg{3}(s) = subplot(nRow,nCol,s); hold on;
    obj1   = shadedErrorBar(x, edBC.BPRaff{s}.irf_Avg, edBC.BPRaff{s}.irf_Sem, 'k', 1);
    obj2   = shadedErrorBar(x, edBC.BPRcor{s}.irf_Avg, edBC.BPRcor{s}.irf_Sem, color_post, 1);
    
    xlim([0 3]);
    
    if isStrSub, title(str_sub{s}); end
end
axes1 = subplot(nRow,nCol,nRow*nCol);
legend(axes1, [obj1.mainLine obj2.mainLine], 'data (pre)','BPR corrected (post)');
titleStr = sprintf('blink-locked pupil time course'); mh = mtit(titleStr);
mh.th.Position(2) = mh.th.Position(2)+0.03;
xh = xlabel(mh.ah, 'time from blink offset (sec)'); set(xh, 'Visible', 'On')
yh = ylabel(mh.ah, 'PD (mm)'); set(yh, 'Visible', 'On')
if isFigSave, savefigpng(figpath, sprintf('blink-locked'), titleStr); end

% all subjects in a figure (unified scale) (baseline corrected)
fh(4) = figure(20101); clf;
for s = 1:nS
    x = edBC.summary{s}.tVec;    
    
    axg{4}(s) = subplot(nRow,nCol,s); hold on;
    obj1   = shadedErrorBar(x, edBC.summary{s}.BPRaff.bl_irfAvg, edBC.summary{s}.BPRaff.bl_irfSem, 'k', 1);
    obj2   = shadedErrorBar(x, edBC.summary{s}.BPRcor.bl_irfAvg, edBC.summary{s}.BPRcor.bl_irfSem, color_post, 1);
    
    xlim([0 3]);   
    axis tight
    
    if isStrSub, title(str_sub{s}); end
end
linkaxes(axg{4},'y')


for s = 1:nS
   axes(axg{4}(s));
   ShowText(sprintf('%d',edBC.BPR{s}.nAllB));
end
axes1 = subplot(nRow,nCol,nRow*nCol);
% legend(axes1, [line1 line2], 'data (pre)','BPR corrected (post)'); % old
legend(axes1, [obj1.mainLine obj2.mainLine], 'data (pre)','BPR corrected (post)'); % new
titleStr = sprintf('blink-locked pupil time course (baselined), and # blinks'); mh = mtit(titleStr);
mh.th.Position(2) = mh.th.Position(2)+0.03;
xh = xlabel(mh.ah, 'time from blink offset (sec)'); set(xh, 'Visible', 'On')
yh = ylabel(mh.ah, 'PD (mm)'); set(yh, 'Visible', 'On')
if isFigSave, savefigpng(figpath, sprintf('blink-locked'), titleStr); end

end