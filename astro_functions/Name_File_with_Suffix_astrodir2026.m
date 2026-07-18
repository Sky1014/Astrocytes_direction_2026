function name_file_with_suffix = Name_File_with_Suffix_astrodir2026(filename)
    suffix = 1;
    while exist(filename, 'file')
        suffix = suffix + 1;
        [path, name, ext] = fileparts(filename);
        filename = strcat(path,'\', name, '_', num2str(suffix), ext);
    end
    name_file_with_suffix = filename;
end
