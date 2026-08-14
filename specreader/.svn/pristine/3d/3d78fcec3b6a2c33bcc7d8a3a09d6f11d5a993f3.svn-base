function y = litcount(filename, literal)
% Search for number of string matches per line.  

fid = fopen(filename, 'rt');
y = 0;
while feof(fid) == 0
   tline = fgetl(fid);
   matches = findstr(tline, literal);
   num = length(matches);
   if num > 0
      y = y + num;
%      fprintf(1,'%d:%s\n',num,tline);
   end
end
fclose(fid); 