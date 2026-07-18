function ImgStack = ReadTiff_astrodir2026(filepath)

    Info=imfinfo(filepath);

    tif='tif';
    format=Info.Format;
    if  (strcmp(format ,tif)==0)
        disp('Input file is not a TIFF image.');
    end

    Slice=size(Info,1);
    Width=Info.Width;
    Height=Info.Height;

    ImgStack=zeros(Height,Width,Slice);

    for i=1:Slice
        ImgStack(:,:,i)=imread(filepath,i);
    end

end
