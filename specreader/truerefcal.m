function truerefcal(varargin)
% TRUEREFCAL Called by trueref to calculate true reflectivity.
%
% Copyright 2004, Zhang Jiang

hFigTrueref = findall(0,'Tag','trueref_Fig');

% --- return if no data imported
if isappdata(hFigTrueref,'data')
   data = getappdata(hFigTrueref,'data');
else
    return;
end

% --- if true reflectivity exists plot and return
if ~isempty(data.tru)
    truePlot(hFigTrueref);
    return;
end
    

if isempty(data.ref)
    return;
end

if length(data.ref) == 1
    data.mrg{1} = data.ref{1};
elseif length(data.ref) > 1 & isempty(data.mrg{1})
    uiwait(msgbox('Merge curves first!','True Relfectivity Error','error','modal'));  
    return;
end

if length(data.pos) == 1
    data.mrg{2} = data.pos{1};
elseif length(data.pos) > 1 & isempty(data.mrg{2})
    uiwait(msgbox('Merge curves first!','True Relfectivity Error','error','modal'));
    return;
end

if length(data.neg) == 1
    data.mrg{3} = data.neg{1};
elseif length(data.neg) > 1 & isempty(data.mrg{3})
    uiwait(msgbox('Merge curves first!','True Relfectivity Error','error','modal'));
    return;
end

setappdata(hFigTrueref,'data',data);

% --- If there is no diffuse scattering, then conside reflectivity as true
if isempty(data.mrg{2}) & isempty(data.mrg{3})
    data.tru = data.mrg{1};
    setappdata(hFigTrueref,'data',data);
    truePlot(hFigTrueref);
    uiwait(msgbox({'No diffuse data imported.';'View current imported reflectivity as true reflectivity.'},'True Relfectivity Warning','warn','modal'));
    return;
end

% --- If diffuse exists, get user defined subtract factor
prompt = {'Enter FACTOR below (reflectivity - diffuse x FACTOR) :'};
dlgAnswer=inputdlg(prompt,'True Reflectivity',1,{'1'});
if isempty(dlgAnswer)
    return;
end
subtractFactor=str2num(dlgAnswer{1}); 
% if invalid factor, error and return;
if isempty(subtractFactor)
    uiwait(msgbox('Invalid factor !','True Relfectivity Error','error','modal'));    
    return;
end

% --- make the counting time for reflectivity and diffuse the same
tempdata = data.mrg;
for iDiff = 2:3
    if ~isempty(tempdata{iDiff})
        tempdata{iDiff}(:,2:3) = tempdata{iDiff}(:,2:3)*subtractFactor;
    end
end

% --- convert to relative error
for iTempdata = 1:length(tempdata)
    if ~isempty(tempdata{iTempdata})
        tempdata{iTempdata}(:,3) = tempdata{iTempdata}(:,3)./tempdata{iTempdata}(:,2);
    end
end

% --- if only positive or negative longitudinal diffuse scattering exists
if isempty(tempdata{2}) | isempty(tempdata{3})
    reflec = tempdata{1};
    if isempty(tempdata{2}) 
        diffuse = tempdata{3};
    else
        diffuse = tempdata{2};
    end
    if reflec(1,1) < diffuse(1,1) 
        diffuse  = [reflec(1,1),diffuse(1,2:3);diffuse];
    end
    if  reflec(end,1) > diffuse(end,1)
        diffuse = [diffuse;reflec(end,1),diffuse(end,2:3)];
    end
    interpDiff(:,1) = reflec(:,1);        % get angle of intercepted diffuse scattering
%    interpDiff(:,2:3) = interp1(diffuse(:,1),diffuse(:,2:3),reflec(:,1),'linear');
    interpDiff(:,2:3) = interp1(diffuse(:,1),diffuse(:,2:3),reflec(:,1),'linear','extrap');
    trueReflec(:,1) = reflec(:,1);                            % true reflectivity angle
    trueReflec(:,2) = reflec(:,2)-interpDiff(:,2);            % true reflectivity intensity
    absReflecError = reflec(:,3).*reflec(:,2);                % absolute reflectivity error
    absInterpDiffError = interpDiff(:,3).*interpDiff(:,2);    % absolute interpolated diffuse error
    trueReflec(:,3) = sqrt(absReflecError.^2+absInterpDiffError.^2)./trueReflec(:,2); % releative true reflectivity error
    data.tru = trueReflec;
end
% --- if both diffuse scattering data exist, calculate the average of the
% diffuse scattering before subtracting from reflectivity data
if ~isempty(tempdata{2}) & ~isempty(tempdata{3})
    reflec = tempdata{1};
    posDiff = tempdata{2};
    negDiff = tempdata{3};
    if reflec(1,1) < posDiff(1,1)
        posDiff  = [reflec(1,1),posDiff(1,2:3);posDiff];
    end
    if reflec(1,1) < negDiff(1,1)
        negDiff  = [reflec(1,1),negDiff(1,2:3);negDiff];
    end
    if  reflec(end,1) > posDiff(end,1)
        posDiff = [posDiff;reflec(end,1),posDiff(end,2:3)];
    end
    if  reflec(end,1) > negDiff(end,1)
        negDiff = [negDiff;reflec(end,1),negDiff(end,2:3)];
    end
    interpPosDiff(:,1) = reflec(:,1);       % set angle of intercepted positive diffuse scattering
    interpNegDiff(:,1) = reflec(:,1);       % set angle of intercepted negative diffuse scattering
    interpPosDiff(:,2:3) = interp1(posDiff(:,1),posDiff(:,2:3),reflec(:,1),'linear','extrap');    
    interpNegDiff(:,2:3) = interp1(negDiff(:,1),negDiff(:,2:3),reflec(:,1),'linear','extrap');            
    trueReflec(:,1) = reflec(:,1);                                                    % true reflectivity angle
    trueReflec(:,2) = reflec(:,2)-0.5*interpPosDiff(:,2)-0.5*interpNegDiff(:,2);      % true reflectivity intensity
    absReflecError = reflec(:,3).*reflec(:,2);                                        % absolute reflectivity error
    absInterpPosDiffError = interpPosDiff(:,3).*interpPosDiff(:,2);                   % absolute interpolated positive diffuse error    
    absInterpNegDiffError = interpNegDiff(:,3).*interpNegDiff(:,2);                   % absolute interpolated negtive diffuse error    
    trueReflec(:,3) = sqrt(absReflecError.^2+0.25*absInterpPosDiffError.^2+0.25*absInterpNegDiffError.^2)./trueReflec(:,2);
    data.tru = trueReflec;   
end
% --- Clear negative data
negPoint = find(data.tru(:,2) < 0);
data.tru(negPoint,:)=[];
% --- convert back to absolute error
data.tru(:,3) = data.tru(:,3).*data.tru(:,2);
% --- save data 
setappdata(hFigTrueref,'data',data);    
truePlot(hFigTrueref);


%===================================================
% --- plot function
%===================================================
function truePlot(hFigTrueref)
hAxes = findall(hFigTrueref,'Tag','trueref_Axes');
data = getappdata(hFigTrueref,'data');
cla(hAxes);
hLine = line('Parent',hAxes,...
    'XData',data.tru(:,1),...
    'YData',data.tru(:,2),...
    'Tag','trueReflectivity');
setappdata(hLine,'ydataError',data.tru(:,3));
set(hLine,...
    'Color','r',...
    'LineStyle','-',...
    'Marker','o',...
    'MarkerSize',3,...
    'MarkerFaceColor','m');
resettoolbar(hFigTrueref);
curvelegend(hFigTrueref);
