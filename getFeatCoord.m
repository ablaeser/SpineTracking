function [feat, Nfeat] = getFeatCoord( XMLpath, Ntype )
%FeatRead reads all the xml files listed in allComp structure and collects the coordinates of each feature.
%fprintf('Reading feature data...  ');
feat = cell(1,Ntype); Nfeat = zeros(1,Ntype); 
Xstruct = xml2struct( XMLpath ); 
try % xml2struct seems to sometimes organize its output differently?
    for w = 1:Ntype % flip( 1:(length(Xstruct.Children(4).Children)-3)/2 ) 
        if isfield( Xstruct.Children(4).Children(2*w+2), 'Name' )
            Nfeat(w) = (numel( Xstruct.Children(4).Children(2*w+2).Children )-3)/2; %#ok<*AGROW>
            feat{w} = zeros( Nfeat(w), 3 );
            if Nfeat(w) == 1
                    feat{w}(1,1) = str2double( Xstruct.Children(4).Children(2*w+2).Children(4).Children(2).Children.Data );
                    feat{w}(1,2) = str2double( Xstruct.Children(4).Children(2*w+2).Children(4).Children(4).Children.Data );
                    feat{w}(1,3) = str2double( Xstruct.Children(4).Children(2*w+2).Children(4).Children(6).Children.Data );
            else
                for q = 1:Nfeat(w)
                    feat{w}(q,1) = str2double( Xstruct.Children(4).Children(2*w+2).Children(2*q+2).Children(2).Children.Data );
                    feat{w}(q,2) = str2double( Xstruct.Children(4).Children(2*w+2).Children(2*q+2).Children(4).Children.Data );
                    feat{w}(q,3) = str2double( Xstruct.Children(4).Children(2*w+2).Children(2*q+2).Children(6).Children.Data );
                end
            end
            [~, sortInd] = sort( feat{w}(:,3), 'ascend' ); 
            feat{w} = feat{w}(sortInd,:); % sort from shallow to deep in stack
        end
        %fprintf('\n%s: Marker type problem (w = %d)', XMLpath, w);
    end
catch
    %
    for w = 1:numel( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type ) %1:Ntype % flip( find( ~cellfun( @isempty, Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type ) ) )  
        if isfield( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}, 'Marker' )
            Nfeat(w) = numel( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker );
            feat{w} = zeros( Nfeat(w), 3 );
            if Nfeat(w) == 1
                feat{w}(1,1) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker.MarkerX.Text );
                feat{w}(1,2) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker.MarkerY.Text );
                feat{w}(1,3) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker.MarkerZ.Text );
            else
                for q = 1:Nfeat(w)
                    feat{w}(q,1) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker{q}.MarkerX.Text );
                    feat{w}(q,2) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker{q}.MarkerY.Text );
                    feat{w}(q,3) = str2double( Xstruct.CellCounter_Marker_File.Marker_Data.Marker_Type{w}.Marker{q}.MarkerZ.Text );
                end
            end
            [~, sortInd] = sort( feat{w}(:,3), 'ascend' ); 
            feat{w} = feat{w}(sortInd,:); % sort from shallow to deep in stack
        else
            feat{w} = zeros(0,3);
        end
    end
end
feat = feat(1:Ntype); Nfeat = Nfeat(1:Ntype);
end