clc
clear

fc=857.5e6;
%% Input an image file and convert to binary stream
load('Picture_all.mat');
index = 1;
fData = Picture_all(index).data;      % Read image data from file
scale = 0.02;                      % Image scaling factor
origSize = size(fData);            % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Resize image
imsize = size(fData);              % Store new image size
binData = dec2bin(fData(:),8);     % Convert to 8 bit unsigned binary
trData = reshape((binData-'0').',1,[]).'; % Create binary stream
Global_Parameters;
[eNodeBOutput,txGrid] = lteRMCDLTool(rmc3,trData);
eNodeBOutput(687:823) = eNodeBOutput(687:823)*3;%%加功率
eNodeBOutput(10287:10423) = eNodeBOutput(10287:10423)*3;
sync = ones(19200,1)*0.01;
sync(687:823) = sync(687:823)*3;%%加功率
sync(10287:10423) = sync(10287:10423)*3;

chn=2;
%% Connect to Radio
radio = comm.SDRuTransmitter('Platform','X310','IPAddress', '192.168.40.2');
radio.MasterClockRate = 200e6;
radio.InterpolationFactor = 80;
if chn==2
    radio.ChannelMapping = [1 2]; 
    tx_sign = [eNodeBOutput sync];
else 
    radio.ChannelMapping = 1; 
    tx_sign = eNodeBOutput;
end
radio.CenterFrequency = fc;
radio.Gain = 0;
radio.UnderrunOutputPort = true;
radio

%% Send Signal over One Antennas

% Scale signal to make maximum magnitude equal to 1
tx_sign = tx_sign/max(abs(tx_sign(:)));
disp('Starting transmission');

currentTime = 0;
while currentTime<1600                    
bufferUnderflow = step(radio,tx_sign);
if bufferUnderflow~=0
    warning('sdru:examples:DroppedSamples','Dropped samples')
end
currentTime = currentTime+0.01;
end
release(radio);
disp('Transmission finished')
displayEndOfDemoMessage(mfilename)


