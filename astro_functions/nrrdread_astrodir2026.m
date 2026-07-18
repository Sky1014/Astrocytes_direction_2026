function [X, meta] = nrrdread_astrodir2026(filename)

fid = fopen(filename, 'rb');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

theLine = fgetl(fid);
assert(numel(theLine) >= 4, 'Bad signature in file.')
assert(isequal(theLine(1:4), 'NRRD'), 'Bad signature in file.')

meta = struct([]);

while (true)

  theLine = fgetl(fid);

  if (isempty(theLine) || feof(fid))
    break;
  end

  if (isequal(theLine(1), '#'))
      continue;
  end

  parsedLine = regexp(theLine, ':=?\s*', 'split','once');

  assert(numel(parsedLine) == 2, 'Parsing error')

  field = lower(parsedLine{1});
  value = parsedLine{2};

  field(isspace(field)) = '';
  meta(1).(field) = value;

end

datatype = getDatatype(meta.type);

assert(isfield(meta, 'sizes') && ...
       isfield(meta, 'dimension') && ...
       isfield(meta, 'encoding') && ...
       isfield(meta, 'endian'), ...
       'Missing required metadata fields.')

dims = sscanf(meta.sizes, '%d');
ndims = sscanf(meta.dimension, '%d');
assert(numel(dims) == ndims);

data = readData(fid, meta, datatype);
data = adjustEndian(data, meta);

X = reshape(data, dims');
X = permute(X, [2 1 3]);

function datatype = getDatatype(metaType)

switch (metaType)
 case {'signed char', 'int8', 'int8_t'}
  datatype = 'int8';

 case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
  datatype = 'uint8';

 case {'short', 'short int', 'signed short', 'signed short int', ...
       'int16', 'int16_t'}
  datatype = 'int16';

 case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
       'uint16_t'}
  datatype = 'uint16';

 case {'int', 'signed int', 'int32', 'int32_t'}
  datatype = 'int32';

 case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
  datatype = 'uint32';

 case {'longlong', 'long long', 'long long int', 'signed long long', ...
       'signed long long int', 'int64', 'int64_t'}
  datatype = 'int64';

 case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
       'uint64', 'uint64_t'}
  datatype = 'uint64';

 case {'float'}
  datatype = 'single';

 case {'double'}
  datatype = 'double';

 otherwise
  assert(false, 'Unknown datatype')
end

function data = readData(fidIn, meta, datatype)

switch (meta.encoding)
 case {'raw'}

  data = fread(fidIn, inf, [datatype '=>' datatype]);

 case {'gzip', 'gz'}

  tmpBase = tempname();
  tmpFile = [tmpBase '.gz'];
  fidTmp = fopen(tmpFile, 'wb');
  assert(fidTmp > 3, 'Could not open temporary file for GZIP decompression')

  tmp = fread(fidIn, inf, 'uint8=>uint8');
  fwrite(fidTmp, tmp, 'uint8');
  fclose(fidTmp);

  gunzip(tmpFile)

  fidTmp = fopen(tmpBase, 'rb');
  cleaner = onCleanup(@() fclose(fidTmp));

  meta.encoding = 'raw';
  data = readData(fidTmp, meta, datatype);

 case {'txt', 'text', 'ascii'}

  data = fscanf(fidIn, '%f');
  data = cast(data, datatype);

 otherwise
  assert(false, 'Unsupported encoding')
end

function data = adjustEndian(data, meta)

[~,~,endian] = computer();

needToSwap = (isequal(endian, 'B') && isequal(lower(meta.endian), 'little')) || ...
             (isequal(endian, 'L') && isequal(lower(meta.endian), 'big'));

if (needToSwap)
    data = swapbytes(data);
end
