function real_depths = calcProbeChDepth(ch_depths,probe_info)

switch num2str(probe_info.probe_id)
    case '10034'
        distance_to_first_contact = 250; %?

    case '30227'
        distance_to_first_contact = 250;

    case '10225'
        distance_to_first_contact = 250;

    case '40479'
        distance_to_first_contact = 500;
        
    case '40482'
        distance_to_first_contact = 500; %?

    otherwise
        distance_to_first_contact = 500;

end


switch probe_info.probe_type
    case 'DBC-32A'
        distance_between_contacts = 65; 

        real_depths = probe_info.mdi_depth_um - ...
            (distance_to_first_contact + ((probe_info.probe_nch-ch_depths) * distance_between_contacts));
end






        



