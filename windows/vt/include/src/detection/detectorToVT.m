function detectorToVT(detector, filename, multiplyThresholdsBy255)

n = detector.clf.nClf;
res=struct('dims',detector.modelDsPad/detector.pPyramid.pChns.shrink,'nWeak',n,...
  'hs',reshape(detector.clf.hs,[],n),'thrs',reshape(detector.clf.thrs,[],n),...
  'fIds',reshape(detector.clf.fIds,[],n),'scoreThr',min(detector.clf.thrs0));
res.thrs=res.thrs(1:3,:); res.fIds=res.fIds(1:3,:);

% undo scaling done inside of chnsCompute, when vtchannels() byte images
% are converted to float. 
if (multiplyThresholdsBy255)
  res.thrs = res.thrs*255;
end

% "ound thresholds
% NOTE: we need to ceil() rather than round() here. This ensures that the 
% ">=" tests (when evaluating trees using integers for thresholds and inputs) 
% have the same outcome as they would have when using floating point math. 
% Example: 1 >= 1.2 <=> 1 >= ceil(1.2)
res.thrs = ceil(res.thrs);

% scale tree outputs and round
outscalebits = 7;
outscale = 2^outscalebits;
res.hs = round(res.hs * outscale);
assert(min(res.hs(:)) >= -outscale);
assert(max(res.hs(:)) <= outscale);

% scale score threshold and round
res.scoreThr = round(res.scoreThr * outscale);

h = res.dims(1,1); w = res.dims(1,2);
for j=1:res.nWeak
  for k=1:3
    res.fIds(k,j) = convertId(res.fIds(k,j),w,h);
  end
end

chns = chnsCompute(uint8(rand(64, 64, 3)),detector.pPyramid.pChns);
numChannels = 0;
for i = 1:chns.nTypes
  numChannels = numChannels + size(chns.data{i},3);
end

fid = fopen(filename, 'wt');
fprintf(fid, '%d %d %d %d %d %d\n', w, h, numChannels, res.nWeak, ...
res.scoreThr, outscalebits);
for j=1:res.nWeak
  for k=1:3
    fprintf(fid, '%d %d ', res.fIds(k,j), res.thrs(k,j));
  end
  fprintf(fid, '%d ', res.hs(:,j));
  fprintf(fid, '\n');
end
fclose(fid);

end

function newId = convertId(id, w, h)

windowSize = w * h;
imgOffset = id - mod(id, windowSize);
rem = mod(id, windowSize);
origCol = idivide(rem, h, 'floor');
origRow = mod(rem, h);
newRem = origRow * w + origCol;
newId = imgOffset + newRem;

end