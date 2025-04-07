function RMP_coils=BS_COMPASS_AUG_RMP_add_signs(RMP_coils,n,parity)
%BS_AUG_RMP_add_signs Determines the sign of a coil based on parity and mode
%   Used to define the sign of each coil in the RMP_coils-struct. This is
%   returned in the RMP_coils struct again

programmed=false;
switch n
    case 1
        %% n=1
        switch parity
            case 'even'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Bu
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,2,3,4}
                            RMP_coils.Bu(i).sign=1;
                        case {5,6,7,8}
                            RMP_coils.Bu(i).sign=-1;
                    end
                end
                % Bl
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {1,2,3,4}
                            RMP_coils.Bl(i).sign=1;
                        case {5,6,7,8}
                            RMP_coils.Bl(i).sign=-1;
                    end
                end
                % A
                for i=1:length(RMP_coils.Au)
                    switch i
                        case {1,2,3,4}
                            RMP_coils.Au(i).sign=1;
                        case {5,6,7,8}
                            RMP_coils.Au(i).sign=-1;
                    end
                end
                for i=1:length(RMP_coils.Al)
                    switch i
                        case {1,2,3,4}
                            RMP_coils.Al(i).sign=1;
                        case {5,6,7,8}
                            RMP_coils.Al(i).sign=-1;
                    end
                end
            case 'odd'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Bu
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,2}
                            RMP_coils.Bu(i).sign=1;
                        case {3,4}
                            RMP_coils.Bu(i).sign=-1;
                    end
                end
                % Bl
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {1,2}
                            RMP_coils.Bl(i).sign=-1;
                        case {3,4}
                            RMP_coils.Bl(i).sign=1;
                    end
                end
                % A
                for i=1:length(RMP_coils.Au)
                    switch i
                        case {1,2}
                            RMP_coils.Au(i).sign=-1;
                        case {3,4}
                            RMP_coils.Au(i).sign=1;
                    end
                end
                for i=1:length(RMP_coils.Al)
                    switch i
                        case {1,2}
                            RMP_coils.Al(i).sign=1;
                        case {3,4}
                            RMP_coils.Al(i).sign=-1;
                    end
                end
        end
    case 2
        %% n=2
        switch parity
            case 'even'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,2,5,6}
                            RMP_coils.Bu(i).sign=1;
                        case {3,4,7,8}
                            RMP_coils.Bu(i).sign=-1;
                    end
                end
                
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {1,2,5,6}
                            RMP_coils.Bl(i).sign=1;
                        case {3,4,7,8}
                            RMP_coils.Bl(i).sign=-1;
                    end
                end
                
                for i=1:length(RMP_coils.A)
                    switch i
                        case {1,2,5,6}
                            RMP_coils.A(i).sign=1;
                        case {3,4,7,8}
                            RMP_coils.A(i).sign=-1;
                    end
                end
            case 'odd'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Bu
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,3}
                            RMP_coils.Bu(i).sign=1;
                        case {2,4}
                            RMP_coils.Bu(i).sign=-1;
                    end
                end
                % Bl
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {1,3}
                            RMP_coils.Bl(i).sign=-1;
                        case {2,4}
                            RMP_coils.Bl(i).sign=1;
                    end
                end
                % A
                for i=1:length(RMP_coils.Au)
                    switch i
                        case {1,3}
                            RMP_coils.Au(i).sign=-1;
                        case {2,4}
                            RMP_coils.Au(i).sign=1;
                    end
                end
                for i=1:length(RMP_coils.Al)
                    switch i
                        case {1,3}
                            RMP_coils.Al(i).sign=1;
                        case {2,4}
                            RMP_coils.Al(i).sign=-1;
                    end
                end
            case '+0'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,5}
                            RMP_coils.Bu(i).sign=1;
                        case {3,7}
                            RMP_coils.Bu(i).sign=-1;
                        otherwise
                            RMP_coils.Bu(i).sign=0;
                    end
                end
                
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {1,5}
                            RMP_coils.Bl(i).sign=1;
                        case {3,7}
                            RMP_coils.Bl(i).sign=-1;
                        otherwise
                            RMP_coils.Bl(i).sign=0;
                    end
                end
                
                for i=1:length(RMP_coils.A)
                    RMP_coils.A(i).sign=0;
                end
            case '+180'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:length(RMP_coils.Bu)
                    switch i
                        case {1,5}
                            RMP_coils.Bu(i).sign=1;
                        case {3,7}
                            RMP_coils.Bu(i).sign=-1;
                        otherwise
                            RMP_coils.Bu(i).sign=0;
                    end
                end
                
                for i=1:length(RMP_coils.Bl)
                    switch i
                        case {3,7}
                            RMP_coils.Bl(i).sign=1;
                        case {1,5}
                            RMP_coils.Bl(i).sign=-1;
                        otherwise
                            RMP_coils.Bl(i).sign=0;
                    end
                end
                
                for i=1:length(RMP_coils.A)
                    RMP_coils.A(i).sign=0;
                end
            
        end
    
end

%% Control if fields have been assigned
if ~programmed && ~isfield(RMP_coils.Bu,'sign')
    error('The requested mode / parity is not programmed')
end
end
