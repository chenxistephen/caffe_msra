vtr = getenv('VisionToolsRoot');

im = imread([vtr '\tests\detectortest\images\finput.png']);
if (size(im,3) == 3)
    im = rgb2gray(im);
end

% need to fix this: channels in matlab have small roundoff errors 
% compared to c++ version..
channels = uint8(255*vtimex('load', ...
    [vtr '\tests\detectortest\images\fchannels.vti']));

figure(1)
colormap jet
subplot(3,3,1);
imagesc(im, [0 255]);
axis image
for i=1:min(size(channels,3),9)
    subplot(3,3,i + 1);
    imagesc(channels(:,:,i), [0 255]);
    axis image
end

figure(2);
vtdetector('init', [vtr '\tests\detectortest\images\fdetector.txt']);
heatmap = vtdetector('run', channels); 
heatmap = heatmap{1};
subplot(2,2,1);
imagesc(im, [0 255]);
axis image
subplot(2,2,2);
imagesc(heatmap);
axis image
hmgt = imread([vtr '\tests\detectortest\images\fheatmap.png']);

vtchannels('init');
pyr = vtchannels('run', im);
for i=1:min(length(pyr),9)
    figure(3);
    subplot(3,3,i);
    imagesc(pyr{i}(:,:,2), [0 100]);
    figure(4);
    subplot(3,3,i);
    tmp = vtdetector('run', pyr{i});
    imagesc(tmp{1}, [120 150]);
    axis image
end

