% makeDirectoryMPP(foldername)
% Makes the folder if is does not exist.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray, 2008 
% Distributed under the General Public License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function makeDirectory(fn)

fn=platformSpecificName(fn);
if isdir(fn)==0
    disp(['Creating directory ',fn]);
    mkdir(fn);
end