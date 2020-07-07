


load('./CA_files/decoder.mat','decoder')
load('./CA_files/mask.mat','mask')
decoder(decoder == -1) = -0.7059;

backgroundRate = 1000;
sourceStrength = 100;
xCoord = 8;
yCoord = 1;

image_array = zeros(16,16);
image_array(xCoord,yCoord) = image_array(xCoord,yCoord) + sourceStrength;
image_array = image_array + poissrnd(backgroundRate/256, [16,16]);

detector_array = conv2(mask, image_array, 'same');

decoded_image = conv2(detector_array, decoder);

figure(1);
contourf(rot90(decoded_image, 3));
colorbar();
