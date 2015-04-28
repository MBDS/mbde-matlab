function [ out ] = mbe_iff(cond, val_true, val_false)
%IFF Ternary conditinal operator, like "cond ? a:b" in C language. 
%   Usage: iff(condition, value_if_true, value_if_false)
% Though, not 100% efficient as in C, since here we force evaluating both sides of
% the condition.

	% -----------------------------------------------------------------------------
	% This file is part of MBDE-MATLAB.  See: https://github.com/MBDS/mbde-matlab
	% 
	%     MBDE-MATLAB is free software: you can redistribute it and/or modify
	%     it under the terms of the GNU General Public License as published by
	%     the Free Software Foundation, either version 3 of the License, or
	%     (at your option) any later version.
	% 
	%     MBDE-MATLAB is distributed in the hope that it will be useful,
	%     but WITHOUT ANY WARRANTY; without even the implied warranty of
	%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	%     GNU General Public License for more details.
	% 
	%     You should have received a copy of the GNU General Public License
	%     along with MBDE-MATLAB.  If not, see <http://www.gnu.org/licenses/>.
	% -----------------------------------------------------------------------------


    if (cond)
        out = val_true;
    else
        out = val_false;
    end
end

