function charPos = getCharPos( filename, charIn )
% getDashUndDot finds positions of selected characters. Used to dissect filenames.
% -----Inputs-------
% filename
% charIn = cell list of characters to find
% -----Outputs------
% charPos = positions of all instances of each character in the filename

for c = numel( charIn ):-1:1
    charPos{c} = strfind( filename, charIn{c} );
end

end

