%%%% Multi-Frequency representation
clear;
clc;
path = 'D:\Motif_Results\Datasets\BirdSong\Features_RMT\';
done = 1;
for tsID = 1 : 154
    if done == 0
        load([path,num2str(tsID),'\feature_',num2str(tsID),'.mat']);
        
        [nanID,~] = find(isnan(data));
        
        data(nanID,:)=[];
        [M,N]= size(data);
        windowsize= 8;
        noOfWindows=M-8+1;
        
        temp = data(:,7);
        sampling_frequency=256;
        %samplesInOneWindow=sampling_frequency*lengthOfWindow;
        %stride= samplesInOneWindow/2;
        %noOfWindows=(totlength) ./ (samplesInOneWindow);% 216 million/2048=105,944
        %noOfWindows=floor(noOfWindows);
        tempSeries = [];
        for i =1 : 13
            yfreq = [];
            for j = 1:noOfWindows
                
                firstWindow= temp( j:j+ windowsize-1);
                %% forst way
                xdft = fft(firstWindow);
                psdx = (1/sampling_frequency*N)* abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                
                %% Matlab way
                %psdx = psd(firstWindow);
                psdx = spectrum(firstWindow);
                
                
                
                %   k = psd(temp( j:j+ windowsize-1));
                %   y{i,j+1} = abs(fftshift(fft(temp( j:j+ windowsize-1)))) ;
                y{i,j+1} = psdx;
                % N is the length of power spectrum=2048
                N = length(y{i,j+1});
                % x is the discrete frequency range -1024/2048 to 1024/ 2048, values are
                % -0.5 to +0.5
                x{i,j+1} =( -N/2 : N/2-1)./N;
                % sampling frequency is 256Hz,rate at which we sample the EEG data
                f = 256;
                % discretizing the frequency space into discrete integer multiples of f
                % x is -128 to +128 Hz
                x{i,j+1} = x{i,j+1} .*f;
                
                %% fignal in DB/HZ
                % y{i,j+1}= 10*log10(psdx)
                
                % grouping of frequency
                % and averaging energy components based on grouped frequency power
                % spectrum
                xf{i,j+1}= [0:1:19];
                yf{i,j+1} =  interp1(x{i,j+1}, y{i,j+1}, xf{i,j+1});
                
                % yf is the energy corresponding to 0 to 19 Hz
                %interpl function interpolates to find y, the values of the underlying function y at the points
                %in the vector or array xf
                
                y1=mean(yf{i,j+1}(1:5)) ;
                yfreq(1,j)=y1;
                yf1{i, j+1}{1,1} = y1;
                
                y2=mean(yf{i,j+1}(6:10)) ;
                yfreq(2,j)=y2;
                yf1{i, j+1}{1,2} = y2;
                
                y3=mean(yf{i,j+1}(11:15)) ;
                yfreq(3,j)=y3;
                yf1{i, j+1}{1,3} = y3;
                
                y4=mean(yf{i,j+1}(16:20))   ;
                yfreq(4,j)=y4;
                yf1{i, j+1}{1,4} = y4;
                %         y5=mean(yf{i,j+1}(17:20)) ;
                %         yfreq(5,j)=y5;
                %         yf1{i, j+1}{1,5} = y5;
                
            end
            tempSeries = [tempSeries;yfreq];
        end
        csvwrite(['D:\Motif_Results\Datasets\BirdSong\ExpandedData\',num2str(tsID),'.csv'],tempSeries);
    end
    data= csvread(['D:\Motif_Results\Datasets\BirdSong\ExpandedData\',num2str(tsID),'.csv']);
    save(['D:\Motif_Results\Datasets\BirdSong\ExpandedData\',num2str(tsID),'.mat'],'data');
end