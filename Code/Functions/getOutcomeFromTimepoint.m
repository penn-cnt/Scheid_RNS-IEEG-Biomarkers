function [outcomes, outcome_date, all_outcomes, all_dates] = getOutcomeFromTimepoint(dt_ref, outcomeCellPair, targ_yr, options)
% USAGE: outs = getOutcomeFromTimepoint(t_ref, outcomeCellPair, t_targ)
%
% INPUTS:
% dt-ref: reference time point, datetime or datetime vector length N
% outcomeCellPair: cell or array of cells containing datetime, outcome pairs
% t-targ: scalar or vector of numbers (years), length T
% options: type = {'interpolate', 'locf'}

% OUTPUT:
% outcomes: N x 1 vector of outcomes for the N t_ref values and T t_targ
%           values. Outcomes are a percent reduction from seizure baseline
% outcome_date: N x 1 datetime array of outcome dates nearest to (dt_ref + targ_yr)
% all_outcomes: N x 1 cell array, each cell contains all outcomes. (This
%               output does not rely on dt_ref or targ_yr.)


arguments
    dt_ref
    outcomeCellPair
    targ_yr
    options.type = ''
end

nSets= length(dt_ref);
outcomes= nan(nSets, length(targ_yr));
outcome_date= NaT(nSets,length(targ_yr));

all_outcomes = cell(nSets,1);
all_dates = cell(nSets, 1); 

for i=1:nSets
    
    % Skip if data is missing
    if isempty(outcomeCellPair{i}), continue, end
    
    % Ensure dates are unique and sorted
    [pt_dates, i_sort] = sort(outcomeCellPair{i}{1});
    if length(pt_dates) == 1
        pt_outs = outcomeCellPair{i}{2};
    else
        pt_outs = outcomeCellPair{i}{2}(i_sort);
        [pt_dates, u1] = unique(pt_dates);
        pt_outs = pt_outs(u1);
    end
     
    % Get outcomes numbers
    if isa(pt_outs, 'double')
        outs = pt_outs;
        
    elseif isa(pt_outs, 'char')
        if strcmp(pt_outs,'No Change'), outs= 0;
        elseif strcmp(pt_outs,'Increased'), outs= -10;
        else error('unrecognized outcome'); end
        
    else % dealing with an array
        mtch = regexp(pt_outs,'\d*','Match');
        outs = zeros(length(mtch), 1); 
        for m=1:length(mtch)    
            if isempty(mtch) || isempty(mtch{m})
                if strcmp(pt_outs{m},'No Change'), outs(m)= 0;
                elseif strcmp(pt_outs{m},'Increased'), outs(m)= -10;
                else error('unrecognized outcome')
                end
            else
                outs(m) = mean(str2double(mtch{m}));
            end
        end
    end

    all_outcomes(i) = {outs};  
    all_dates(i) = {datetime(pt_dates, 'format', 'MM-dd-yyyy')}; 
    
    for i_targ= 1:length(targ_yr)
    % Get years b/w implant and measurement, and closest date to target 
        yrs= years(pt_dates-dt_ref(i));
        
        % Interpolate between datapoints
        if strcmp(options.type, 'interpolate')
            dateOut= dt_ref(i) + years(targ_yr(i_targ));
            interpOut = interp1(yrs, outs, targ_yr(i_targ));

            if isnan(interpOut) && targ_yr(i_targ) < min(yrs)
                [yrmin, i_min] = min(yrs);
                interpOut = outs(i_min); dateOut = dt_ref(i)+ years(yrmin);
            elseif isnan(interpOut) && targ_yr(i_targ) > max(yrs)
                [yrmax, i_max] = max(yrs);
                interpOut = outs(i_max); dateOut = dt_ref(i)+ years(yrmax);
            end    
            
            outcomes(i, i_targ)= round(interpOut);
            outcome_date(i, i_targ)= dateOut;
        
        % Report last observation carried forward to point
        elseif strcmp(options.type, 'locf')
            dateOut= dt_ref(i) + years(targ_yr(i_targ));
            
            i_locf = find(pt_dates <= dateOut, 1, 'last');
            if isempty(i_locf)
                warning('target time must be at least %0.2f years for set %d', min(years((pt_dates-dt_ref(i)))), i)
                continue
            end
            
            outcome_date(i, i_targ) = pt_dates(i_locf);
            outcomes(i, i_targ)= round(outs(i_locf));
            
        else
            [~,i_min]= min(abs(yrs-targ_yr(i_targ)));
            outcome_date(i, i_targ) = pt_dates(i_min);
            outcomes(i, i_targ)= round(outs(i_min));  
        end
    end

end


end