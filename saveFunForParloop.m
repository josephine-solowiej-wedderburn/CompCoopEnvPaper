function saveFunForParloop( filename_as_string, var )

%note the file name should have the full path included and end .mat (or
%alternative)
%variable should NOT input as a string

var_name = inputname(2);
S.(var_name) = var;

    save(filename_as_string, '-struct', 'S')

end