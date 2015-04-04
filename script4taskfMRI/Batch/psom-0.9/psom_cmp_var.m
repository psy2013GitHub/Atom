function flag_equal = psom_cmp_var(var1,var2)
%
% _________________________________________________________________________
% SUMMARY PSOM_CMP_VAR
%
% Compare two matlab variables, whatever their types may be
%
% SYNTAX:
% FLAG_EQUAL = PSOM_CMP_VAR(VAR1,VAR2)
%
% _________________________________________________________________________
% INPUTS:
%
% VAR1
%       (any type) a matlab variable
%
% VAR2
%       (any type) a matlab variable
%
% _________________________________________________________________________
% OUTPUTS:
%
% FLAG_EQUAL
%       (boolean) true if VAR1 and VAR2 are identical, false otherwise.
%
% _________________________________________________________________________
% COMMENTS:
%
% The actual matlab data types supported here are string (or char),
% numeric, logical,
% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008.
% Maintainer : pbellec@bic.mni.mcgill.ca
% See licensing information in the code.
% Keywords :

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%% Get type of variable 1

type_var1 = sub_type_var(var1);
type_var2 = sub_type_var(var2);

if ~strcmp(type_var1,type_var2)
    flag_equal = false;
    return
else
    switch type_var1

        case 'char'
            flag_equal = strcmp(var1,var2);

        case {'numeric','logical'}
            if min(size(var1) == size(var2)) == 0
                flag_equal = false;
                return
            else
                var1 = var1(:);
                var2 = var2(:);
                mask_nan = isnan(var1);
                flag_equal = false(size(var1));
                flag_equal(mask_nan) = isnan(var2(mask_nan));
                flag_equal(~mask_nan) = var1(~mask_nan) == var2(~mask_nan);
                flag_equal = min(flag_equal);                
            end

        case 'cell'

            if min(size(var1) == size(var2)) == 0
                flag_equal = false;
                return
            else
                var1 = var1(:);
                var2 = var2(:);
                flag_equal = true;
                for num_e = 1:length(var1)
                    if ~psom_cmp_var(var1{num_e},var2{num_e})
                        flag_equal = false;
                        return
                    end
                end
            end

        case 'struct'

            %% compare sizes
            if min(size(var1)==size(var2))==0
                flag_equal = false;
                return
            end
            
            var1 = var1(:);
            var2 = var2(:);
            
            if length(var1) > 1
                %% if there are more than one entry, loop over all entries
                flag_equal = true;
                for num_e = 1:length(var1)
                    flag_equal = flag_equal&&psom_cmp_var(var1(num_e),var2(num_e));
                end
                return
            else
                
                %% compare field names
                list_fields1 = fieldnames(var1);
                list_fields2 = fieldnames(var2);

                if ~psom_cmp_var(list_fields1,list_fields2)
                    flag_equal = false;
                    return
                end

                %% Compare the values of all fields 
                flag_equal = true;
                for num_e = 1:length(list_fields1)
                    if ~psom_cmp_var(var1.(list_fields1{num_e}),var2.(list_fields2{num_e}))
                        flag_equal = false;
                        return
                    end
                end

            end
    end
end

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%
function type_var = sub_type_var(var)

if ischar(var)
    type_var = 'char';
elseif isnumeric(var)
    type_var = 'numeric';
elseif islogical(var)
    type_var = 'logical';
elseif iscell(var)
    type_var = 'cell';
elseif isstruct(var)
    type_var = 'struct';
else
    var
    error('I could not determine the type of the preceeding variable')
end

