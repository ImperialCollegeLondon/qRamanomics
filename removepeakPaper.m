%% This function is used to remove the cosmic rays within the input Raman spectrum
function [x,y]=removepeak(x,y)
    dy=diff(y)./diff(x); % Since the cosmic-ray peak is very sharp, we use the gradient to decide whether is it a cosmic ray
    while max(dy)>8 % The threshold of gradient is set to 8, any peak with a gradient larger than 8 is a cosmic-ray peak
       [~,index]=max(dy);
       y((index-5):(index+5))=min(y); % If the gradient is larger than the threshold, the nearby area will be be smoothed, and thus the cosmic ray is removed
       dy((index-10):(index+10))=min(dy);
    end
end
