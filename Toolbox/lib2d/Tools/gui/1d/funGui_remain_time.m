
function timestr = funGui_remain_time(step,nstep)
		% Returns remaing time string for an equivalent step process of
		% nstep at step. Useful for waitbar update or whatever.
		% tic must have been run once.
		% -------------------------------------------------
		% Denis Rouzaud, TOPO, EPFL, 2009
		% -------------------------------------------------
		ellapsed_time = toc;
		
		unit_proc_time = ellapsed_time / step;
		remain_proc = nstep-step;
		
		remain_time = unit_proc_time * remain_proc;
		
		remain_time_v = zeros(1,4); % d h m s
		
		sm   = 60;
		sh   = 60*sm;
		sday = 24*sh;
		
		remain_time_v(1) = floor( ( remain_time                            ) / sday );		% days
		remain_time_v(2) = floor( ( remain_time - (remain_time_v(1)* sday) ) / sh   );		% hours
		remain_time_v(3) = floor( ( remain_time - (remain_time_v(2)* sh  ) ) / sm   );		% min
		remain_time_v(4) = ceil ( ( remain_time - (remain_time_v(3)* sm  ) ) / 1    );		% sec
				
		%use only 2 values max! i.e. sec are not really intersting if more than 2 days! ;)
		ix = find(remain_time_v~=0, 1, 'first');
		istr = {' day(s)', ' h', ' min', ' sec'};
		
		timestr = [];
		if ix
			if ix < 4
				timestr = [sprintf('%3.0f', remain_time_v(ix)) char(istr(ix))];
			end
			ixe = min(4,ix+1);
			timestr = [timestr sprintf('%3.0f', remain_time_v(ixe)) char(istr(ixe))];
		end
	end
