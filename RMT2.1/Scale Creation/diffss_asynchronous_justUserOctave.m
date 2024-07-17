function dss = diffss_asynchronous_justUserOctave(ss,ODepd,OTime)
dss.otmin = ss.otmin;
dss.odmin = ss.odmin ;
dss.Od = ss.Od ;
dss.Ot = ss.Ot;
%dss.S = ss.S ;
dss.St = ss.St;
dss.Sd = ss.Sd;
dss.sigmat = ss.sigmat ;
dss.sigmad = ss.sigmad ;
% 
% for ODepd = 1 : size(ss.octave,2)
%     for OTime = 1 : size(ss.octave,1)
        [M, N, Time, Depd] = size(ss.octave{OTime, ODepd});
        RowVector = zeros(M, N, Depd-1); % depd difference
        ColumnVector = zeros(M, N, Time-1); % time difference
        DiaVectorDown = zeros(M, N, (Depd-1)); % DoGs from both dimensions
        for i = 1: (Depd-1)
            RowVector(:, :, i) = ss.octave{OTime, ODepd}( :, :, i+1, i)-ss.octave{OTime, ODepd}( : ,:, i, i);
        end
        
        for i = 1 : (Time-1)
            ColumnVector(:,:,i) = ss.octave{OTime, ODepd}( :, :, i, i+1)-ss.octave{OTime, ODepd}( : ,:, i, i);
        end
        
        for i = 1 : (Depd-1)
            % DiaVectorDown(:,:,i) = ss.octave{OTime, ODepd}( :, :, i+1, 1)-ss.octave{OTime, ODepd}( : ,:, i, 1);
            DiaVectorDown(:,:,i) = ss.octave{OTime, ODepd}( :, :, i+1, i+1)-ss.octave{OTime, ODepd}( : ,:, i, i);
        end
        
        dss.octave{OTime, ODepd}{1} = RowVector; % depd difference
        % dss.octave{OTime, ODepd}{1} = reshape(RowVector, M, N, Time*(Depd-1));
        dss.octave{OTime, ODepd}{2} = ColumnVector; % time difference
        % dss.octave{OTime, ODepd}{2} = reshape(ColumnVector, M, N, Time*(Depd-1));
        dss.octave{OTime, ODepd}{3} = DiaVectorDown; % both directions
        % dss.octave{OTime, ODepd}{3} = reshape(DiaVectorDown, M, N, (Time-1)*(Depd-1));
        clear RowVector ColumnVector DiaVector  vector row column dia
%     end
% end
