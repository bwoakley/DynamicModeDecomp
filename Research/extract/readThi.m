st='../Cases/Chmo/Thi';
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,251,2);