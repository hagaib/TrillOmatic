function btime = longtrill_bulk_of_mass_weights(persegs , yin , energy , Nstds , returntypeflag)
%calculates start time and end time and indices of bulk of mass of longtrill.
%i.e. time which is most likely to contain most of the trill present in
%the recording

%Nstds = number of stds to take in as bulk of mass
%returntypeflag =   0 - return indices in yin.time
%                   1 - return time in seconds

%example: [istart , iend] = longtrill_bulk_of_mass(persegs , yin , 1.5 , 0);

times = mean(persegs,1);
mtimes = mean(times);
stdtimes = std(times);
tstart = mtimes-stdtimes*Nstds;
tend = mtimes+stdtimes*Nstds;
tstart = persegs(1,find(tstart<persegs(1,:),1,'first'));
tend = persegs(2, find(tend>persegs(2,:),1,'last') );
btime = [tstart ,tend];

disp(['long trill in interval:' , num2str(tstart) , ':' , num2str(tend)]);
if(returntypeflag == 0)
    istart = find(yin.time < tstart , 1 ,'last');
    iend = find(yin.time > tend , 1 , 'first');
    btime = [istart , iend];
end

