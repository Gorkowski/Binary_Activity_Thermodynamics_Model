function [varargout]=progressbartext(progressfrac, rewrite_display)
%% Updated by Kyle Gorkowski [CMU_LAPTOP] on 2016-Aug-02  1:39 PM
%  a simpler version for command window display
%
%Description:
%
%
%   progressbar() provides an indication of the progress of some task using
% graphics and text. Calling progressbar repeatedly will update the figure and
% automatically estimate the amount of time remaining.
%   This implementation of progressbar is intended to be extremely simple to use
% while providing a high quality user experience.
%
% Features:
%   - Can add progressbar to existing m-files with a single line of code.
%   - Supports multiple bars in one figure to show progress of nested loops.
%   - Optional labels on bars.
%   - Figure closes automatically when task is complete.
%   - Only one figure can exist so old figures don't clutter the desktop.
%   - Remaining time estimate is accurate even if the figure gets closed.
%   - Minimal execution time. Won't slow down code.
%   - Randomized color. When a programmer gets bored...
%
% Example Function Calls For Single Bar Usage:
%   progressbar               % Initialize/reset
%   progressbar(0)            % Initialize/reset
%   progressbar('Label')      % Initialize/reset and label the bar
%   progressbar(0.5)          % Update
%   progressbar(1)            % Close
%

% Notes:
%   For best results, call progressbar with all zero (or all string) inputs
% before any processing. This sets the proper starting time reference to
% calculate time remaining.
%   Bar color is choosen randomly when the figure is created or reset. Clicking
% the bar will cause a random color change.
%
% Demos:
%     % Single bar
%     m = 500;
%     progressbartext % Init single bar
%     for i = 1:m
%       pause(0.01) % Do something important
%       progressbartext(i/m) % Update progress bar
%     end

if not(exist('rewrite_display'))
    % set default options
    rewrite_display='yes';
end

persistent progdata

% Init reset flag
resetflag = false;

% Set reset flag if first input is a string
if ischar(progressfrac)
    resetflag = true;
    disp(['Progress... ' progressfrac])
elseif progressfrac==0
    disp('Progress... ')
    resetflag = true;
end

if resetflag
    progdata.starttime = now;
    progdata.time=now;
    progdata.fractiondone=0;
    progdata.learningrate=0.9; % you can change this >0 and <1 addjusts how fast the moving slope is adjusted with longer/shor computation times
    progdata.slope=0;
    progdata.Ldisplay=0;
    progdata.min_seconds_between_calls=.5; % min time between calls to display and calc output
    progdata.previous_call_display=now;
end

% disp('123456789')
% fprintf(repmat('\b',1,10))
% disp('123456789')

% git time wall clock
current_time=now;
wirte_update='Timer Started';


% calc time since last call
time_since_last_call_sec=(current_time-progdata.previous_call_display).*86400;
good_to_calcNdisplay=time_since_last_call_sec>progdata.min_seconds_between_calls;

% Process inputs and update state of progdata
if good_to_calcNdisplay && not(resetflag) && progressfrac>0 && progressfrac<1
    
    % calc slope of time change
    
    if progdata.fractiondone==0 % first time
        current_slope =  (current_time-progdata.time)/(progressfrac-progdata.fractiondone);
        
    else % all times after with learning rate averaging
        slope_i = (current_time-progdata.time)/(progressfrac-progdata.fractiondone); % slope
        current_slope=progdata.slope.*progdata.learningrate+(1-progdata.learningrate).*slope_i; % learning rate averaging
    end
    
    time_left=(1-progressfrac).*current_slope;
    
    timestr=sec2timestr(time_left.*86400);
    
    wirte_update=['Percent Done ' num2str(progressfrac.*100) '% Time Left ' timestr];
    
    if progdata.Ldisplay==0
        disp(wirte_update)
    else
        if strcmpi(rewrite_display, 'yes')
            fprintf(repmat('\b',1,progdata.Ldisplay+1))
        else
%             disp('\n')
        end
        disp(wirte_update)
    end
    
    % save data for next
    progdata.slope=current_slope;
    progdata.fractiondone = progressfrac;
    progdata.time=current_time;
    progdata.Ldisplay=length(wirte_update);
    progdata.previous_call_display=current_time;
    
end



if not(resetflag) && progressfrac==1
    
    total_time=current_time-progdata.starttime;
    done_str=sec2timestr(total_time.*86400);
    
    if strcmpi(rewrite_display, 'yes')
        fprintf(repmat('\b',1,progdata.Ldisplay+1))
    end
        
    wirte_update=[ 'Finished at ' datestr(now,'mm/dd/yy ddd HH:MM PM') '  Total Time ' done_str ];
    disp(wirte_update)
end

% output options
output_call_number=nargout;

if output_call_number==1
    if good_to_calcNdisplay
        varargout{1}=wirte_update;
    else
        varargout{1}=[];
    end
end




% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

% Convert seconds to other units
w = floor(sec/604800); % Weeks
sec = sec - w*604800;
d = floor(sec/86400); % Days
sec = sec - d*86400;
h = floor(sec/3600); % Hours
sec = sec - h*3600;
m = floor(sec/60); % Minutes
sec = sec - m*60;
s = floor(sec); % Seconds

% Create time string
if w > 0
    if w > 9
        timestr = sprintf('%d week', w);
    else
        timestr = sprintf('%d week, %d day', w, d);
    end
elseif d > 0
    if d > 9
        timestr = sprintf('%d day', d);
    else
        timestr = sprintf('%d day, %d hr', d, h);
    end
elseif h > 0
    if h > 9
        timestr = sprintf('%d hr', h);
    else
        timestr = sprintf('%d hr, %d min', h, m);
    end
elseif m > 0
    if m > 9
        timestr = sprintf('%d min', m);
    else
        timestr = sprintf('%d min, %d sec', m, s);
    end
else
    timestr = sprintf('%d sec', s);
end
