fid = fopen([pathname filename '.hdr']);

for k=1:23
    tline=fgets(fid);
    if k==10 % Profile Interval, s
        disp(tline)
        ProfileInterval=cell2mat(textscan(tline,'%*s %*s %f %*s'));
       
    elseif k==12 % Cell size, m
        disp(tline)
        cellsize=cell2mat(textscan(tline,'%*s %*s %f %*s'))/100;
        
    elseif k==16 % Blanking distance, m
        disp(tline)
        blankdist=cell2mat(textscan(tline,'%*s %*s %f %*s'));
    end
end

fclose(fid);
