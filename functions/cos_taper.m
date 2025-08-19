% cos_taper.m
% Applies a 5% cosine taper
% usage:
% tapered = cos_taper(data);
function tapered = cos_taper(data)
    M=floor(((length(data))*5/100)/2+0.5);

    tapered=zeros(1,length(data));
    
    for j=1:length(data)
        if j<=M+1
            tapered(j)=data(j) * (0.5 * ( 1-cos(j*pi/(M+1))));
        elseif (j<length(data) - M-1)
            tapered(j) = data(j);
        elseif j<=length(data)
            tapered(j) = data(j) * (0.5 * (1-cos((length(data)-j)*pi/(M+1))));
        end
    end
    
                
    return
