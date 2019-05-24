function[im]=retrieve_im(im0,IM0,IM)

[r,c]= size(im0);
NIM  = length(IM);
im = zeros(r,NIM);

% only one im vector is provided.  Not wise though
if c==1
    im = repmat(im0,1,length(IM));
    return;
end

for i=1:NIM
   if IM(i)<0 
       ind = IM0==IM(i);
       im(:,i)=im0(:,ind);
   else
       disc = abs(IM0-IM(i));
       [~,ind] = min(disc);
       im(:,i)=im0(:,ind);
   end
end
