function um = feat2um( feat, conVec )
% convert feat (image coord system) to um (micron coord system) using conversion vector (x,y,z) resolution in microns

if isempty( feat )
    um = cell(1,7);
else
    um = cell( size(feat) );
    for w = reshape( find( ~cellfun( @isempty, feat ) ), 1, [] ) 
        um{w} = repmat( conVec, size(feat{w},1), 1).*feat{w};
    end
end

end

