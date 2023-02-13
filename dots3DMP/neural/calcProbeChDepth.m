function real_depths = calcProbeChDepth(MDI_depth,ch_depths,probeID)

if nargin<3, probeID = 'DBC10034'; end

switch probeID
    case 'DBC10034'

        distance_to_first_contact = 750;
        distance_between_contacts = 65;
        numchs = 32;

%     case 'DBC'

end

real_depths = MDI_depth - (distance_to_first_contact + ((numchs-ch_depths) * distance_between_contacts));






        



