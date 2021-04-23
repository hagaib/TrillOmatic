function xlim_binder( axes_array )
% binds the x - axis of all axes in axes_array to scroll together

for i=1:length(axes_array)
    addlistener(axes_array(i) , 'XLim' , 'PostSet' ,@(src , evnt)handleXLimChange(src , evnt, axes_array));
end



function handleXLimChange(src , evnt, axes_array)
% disp('xlim changed')
ax = evnt.AffectedObject;
for i=1:length(axes_array)
    if (~strcmp(ax.Tag , axes_array(i)))
        axes_array(i).XLim = ax.XLim;
    end
end

% function handleXLimChange(src , evnt, handles)
% % disp('xlim changed')
% ax = evnt.AffectedObject;
% if (strcmp(ax.Tag , 'axes_time'))
%     handles.axes_spect.XLim = ax.XLim;
% else
%     handles.axes_time.XLim = ax.XLim;
% end